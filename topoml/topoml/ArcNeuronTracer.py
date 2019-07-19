# Standard library imports
import json
from networkx.readwrite import json_graph
import os

# Third party imports
import numpy as np
from sklearn import linear_model, model_selection
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
import networkx as nx
import scipy.stats


# Local application imports
from topoml.ml.utils import gaussian_fit
from topoml.NeuronTracer import NeuronTracer
from topoml.topology.utils import (
    get_pixel_values_from_arcs,
    is_ridge_arc,
    is_valley_arc,
)
from topoml.image.feature import (
    mean_filter,
    variance_filter,
    median_filter,
    minimum_filter,
    maximum_filter,
    gaussian_blur_filter,
    difference_of_gaussians_filter,
    sobel_filter,
    laplacian_filter,
    neighbor_filter,
)
from topoml.ui.colors import blue, light_blue, red, light_red, green, light_green, purple, light_purple, orange, light_orange
from topoml.ui.vis import show_image


class ArcNeuronTracer(NeuronTracer):
    def __init__(
            self, fname=None, blur_sigma=2, persistence=0,
            valley=True, ridge=True,
            construct_msc = True
    ):
        self.primal = None
        self.dual = None
        self.use_valley_arcs = valley
        self.use_ridge_arcs = ridge
        self.images = {}
        super(ArcNeuronTracer, self).__init__(fname, blur_sigma,
                                              persistence, construct_msc)
        if not construct_msc:
            self.in_arcs = None
            self.out_arcs = None
            self.msc = None
            self.arcs = None
        self.arcs = []
        self.positive_arcs = None
        self.negative_arcs = None
        if construct_msc:
            for arc in self.msc.arcs:
                if (
                        not self.use_ridge_arcs and is_ridge_arc(arc, self.msc)
                ) or (
                    not self.use_valley_arcs and is_valley_arc(arc, self.msc)
                ):
                    continue
                self.arcs.append(arc)
                    

    def _add_default_features(self):
        # Functions to apply to the pixels of an arc
        self.features["length"] = lambda pixels: len(pixels)
        self.features["min"] = lambda pixels: np.min(pixels)
        self.features["max"] = lambda pixels: np.max(pixels)
        self.features["median"] = lambda pixels: np.median(pixels)
        self.features["mode"] = lambda pixels: scipy.stats.mode(np.round(pixels, 2))[0][0]
        self.features["mean"] = lambda pixels: gaussian_fit(pixels)[0]
        self.features["std"] = lambda pixels: gaussian_fit(pixels)[1]
        self.features["var"] = lambda pixels: gaussian_fit(pixels)[2]
        self.features["skew"] = lambda pixels: scipy.stats.skew(pixels)
        self.features["kurtosis"] = lambda pixels: scipy.stats.kurtosis(pixels)
        self.features["range"] = lambda pixels: np.max(pixels) - np.min(pixels)

        # The input for these is fundamentally different, so for now
        # we will key off the fact that the name will always start with
        # "neighbor" in order to pass the right pixels to them.
        # Alternatively, we can send all of the pixels, and in the
        # methods above, just operate on the first row which could be
        # guaranteeed to be the pixels for the arc being operated on.
        self.features["neighbor_degree"] = lambda connected_pixels: len(
            connected_pixels
        )
        self.features["neighbor_min"] = lambda connected_pixels: np.min(
            [np.min(pixels) for pixels in connected_pixels]
        )
        self.features["neighbor_max"] = lambda connected_pixels: np.max(
            [np.max(pixels) for pixels in connected_pixels]
        )
        self.features["neighbor_mean"] = lambda connected_pixels: np.mean(
            [np.mean(pixels) for pixels in connected_pixels]
        )
        self.features["neighbor_std"] = lambda connected_pixels: np.std(
            [np.mean(pixels) for pixels in connected_pixels]
        )

        # Pixel values to use for the aforementioned functions:
        self.images["identity"] = self.image
        self.images["sobel"] = sobel_filter(self.image)
        self.images["laplacian"] = laplacian_filter(self.image)
        for i in range(1, 5):
            pow1 = 2 ** (i - 1)
            pow2 = 2 ** i
            self.images["mean_{}".format(i)] = mean_filter(self.image, i)
            self.images["variance_{}".format(i)] = variance_filter(
                self.image, i
            )
            self.images["median_{}".format(i)] = median_filter(self.image, i)
            self.images["min_{}".format(i)] = minimum_filter(self.image, i)
            self.images["max_{}".format(i)] = maximum_filter(self.image, i)
            self.images["gauss_{}".format(pow2)] = gaussian_blur_filter(
                self.image, pow2
            )
            self.images[
                "delta_gauss_{}_{}".format(pow1, pow2)
            ] = difference_of_gaussians_filter(self.image, pow1, pow2)
        for i, neighbor_image in enumerate(neighbor_filter(self.image, 3)):
            self.images["shift_{}".format(i)] = neighbor_image

    def compile_features(self, selection=None, return_labels=False):
        #G = self.construct_dual_graph()
        G = self.construct_dual_graph()
        arc_map = {}
        arc_pixel_map = {}
        for i, arc in enumerate(self.arcs):
            index = tuple(arc.node_ids) + (len(arc.line),)
            arc_map[index] = i
        #G = self.construct_dual_graph()
        arc_features = []
        feature_names = []
        for arc in self.arcs:
            index = tuple(arc.node_ids) + (len(arc.line),)
            i = arc_map[index]

            arc_feature_row = []

            for image_name, image in self.images.items():
                if i not in arc_pixel_map:
                    arc_pixel_map[i] = get_pixel_values_from_arcs([arc], image)
                arc_pixels = arc_pixel_map[i]

                connected_pixels = []
                for j in G.neighbors(i):
                    arc = self.arcs[j]
                    if j not in arc_pixel_map:
                        arc_pixel_map[j] = get_pixel_values_from_arcs(
                            [arc], self.image
                        )
                    j_pixels = arc_pixel_map[j]
                    connected_pixels.append(j_pixels)

                for function_name, foo in self.features.items():
                    if function_name.startswith("neighbor_"):
                        arc_feature_row.append(foo(connected_pixels))
                    else:
                        arc_feature_row.append(foo(arc_pixels))
                    if len(arc_features) == 0:
                        feature_names.append(image_name + "_" + function_name)
            for node_id in arc.node_ids:
                if self.msc.nodes[node_id].index == 1:
                    saddle_value = self.msc.nodes[node_id].value
                else:
                    maximum_value = self.msc.nodes[node_id].value
            arc_feature_row.append(maximum_value - saddle_value)
            if len(arc_features) == 0:
                feature_names.append("persistence")
            arc_features.append(arc_feature_row)

        arc_features = np.array(arc_features)
        mu = np.mean(arc_features, axis=0)
        std = np.std(arc_features, axis=0)
        arc_features = (arc_features - mu) / std

        if selection is not None:
            
            selected_arc_indices = []
            for i, arc in enumerate(selection):
                index = tuple(arc.node_ids) + (len(arc.line),)
                selected_arc_indices.append(arc_map[index])
            if return_labels:
                return arc_features[selected_arc_indices, :], feature_names
            else:
                return arc_features[selected_arc_indices, :]

        if return_labels:
            return arc_features, feature_names
        else:
            return arc_features

    def _add_default_models(self):
        logit = linear_model.LogisticRegression(
            penalty="l2",
            dual=False,
            tol=0.0001,
            C=1.0,
            fit_intercept=True,
            intercept_scaling=1,
            class_weight=None,
            random_state=None,
            solver="liblinear",
            max_iter=100,
            verbose=0,
            warm_start=False,
            n_jobs=1,
        )
        param_grid = {"penalty": ["l1", "l2"], "C": [1.0, 10]}
        self.models["logit"] = model_selection.GridSearchCV(logit, param_grid)

    def _extract_selection(self, in_arcs, out_arcs, out_pixels):
        self.positive_arcs = in_arcs
        self.negative_arcs = out_arcs
        return in_arcs, out_arcs

    def show_classifier_result(self, classifier):
        compiled_data = self.compile_features()
        arc_classes = classifier.predict(compiled_data)

        predicted_classification = np.zeros(self.image.shape)
        pruned_arcs = None
        kept_arcs = None
        for arc, label in zip(self.arcs, arc_classes):
            points = np.array(np.round(arc.line), dtype=int)
            if label == 1:
                predicted_classification[points[:, 1], points[:, 0]] = label
                if kept_arcs is None:
                    kept_arcs = points
                else:
                    kept_arcs = np.vstack((kept_arcs, points))
            else:
                if pruned_arcs is None:
                    pruned_arcs = points
                else:
                    pruned_arcs = np.vstack((pruned_arcs, points))

        show_image(predicted_classification)
        plt.show(block=True)
        show_image(self.image)
        if pruned_arcs is not None:
            plt.scatter(
                pruned_arcs[:, 0],
                pruned_arcs[:, 1],
                facecolor=blue,
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.03,
            )
        if kept_arcs is not None:
            plt.scatter(
                kept_arcs[:, 0],
                kept_arcs[:, 1],
                facecolor=light_green,
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.03,
            )
        plt.show(block=True)

    def show_classified_graph(self, graph_prefix = '', G = None, train_view = False):
        #compiled_data = self.compile_features()
        #arc_classes = classifier.predict(compiled_data)

        predicted_classification = np.zeros(self.image.shape)

        if graph_prefix:
            predicted_path =  os.path.join(cwd,'data','json_graphs',graph_prefix+"-G.json")
            G = json_graph.node_link_graph(json.load(open(predicted_path)))
        pruned_arcs = None
        kept_arcs = None
        kept_bg = None
        switched_to_arc = None
        switched_to_bg = None
        to_arc_count = 0
        to_bg_count = 0
        same_count= 0
        arc_class=None
        bg_class=None
        val_set = None
        for arc, node in zip(self.arcs, G.nodes_iter()):
            points = np.array(np.round(arc.line), dtype=int)
            if train_view:
                if 'train' in G.node[node]:
                    if G.node[node]['train']:
                        if arc_class is None:
                            arc_class = points
                        else:
                            arc_class = np.vstack((arc_class, points))
                    if G.node[node]['val']:
                        if bg_class is None:
                            bg_class = points
                        else:
                            bg_class = np.vstack((bg_class, points))
            elif "prediction" in G.node[node]:
                if G.node[node]['label'] == G.node[node]['prediction']:
                    if arc_class is None:
                        arc_class = points
                    else:
                        arc_class = np.vstack((arc_class, points))
                if G.node[node]['label'] != G.node[node]['prediction']:
                    if bg_class is None:
                        bg_class = points
                    else:
                        bg_class = np.vstack((bg_class, points))
                """
                # remained arc
                if G.node[node]['label'][1] == 1 and G.node[node]["prediction"][1] == 1: #'prediction
                    predicted_classification[points[:, 1], points[:, 0]] = 1
                    if kept_arcs is None:
                        kept_arcs = points
                    else:
                        kept_arcs = np.vstack((kept_arcs, points))
                if G.node[node]['label'][0] ==  1 and G.node[node]["prediction"][0] == 1: #'prediction
                    predicted_classification[points[:, 1], points[:, 0]] = 0
                    if kept_bg is None:
                        kept_bg = points
                    else:
                        kept_bg = np.vstack((kept_bg, points))
                # predicted arc when originaly background
                if G.node[node]['label'][1] == 0 and G.node[node]["prediction"][1] == 1: #'prediction
                    predicted_classification[points[:, 1], points[:, 0]] = 1.0
                    if switched_to_arc is None:
                        switched_to_arc = points
                    else:
                        switched_to_arc = np.vstack((switched_to_arc, points))
                # predicted background when originaly arc
                if G.node[node]['label'][1] == 1 and G.node[node]["prediction"][1] == 0: #'prediction
                    predicted_classification[points[:, 1], points[:, 0]] = 0
                    if switched_to_bg is None:
                        switched_to_bg = points
                    else:
                        switched_to_bg = np.vstack((switched_to_bg, points))
                """
            elif 'train' in G.node[node] and G.node[node]['train']=='fish':
                if val_set is None:
                    val_set = points
                else:
                    val_set = np.vstack((val_set, points))
            else:
                #print(G.node[node])
                if pruned_arcs is None:
                    pruned_arcs = points
                else:
                    pruned_arcs = np.vstack((pruned_arcs, points))

        #show_image(predicted_classification)
        plt.subplot(211)
        im_light = (self.image*1.3)
        show_image(im_light/np.max(im_light))#((self.image+3)/10.)/(np.max(self.image+3)/10.))
        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        if pruned_arcs is not None:
            plt.scatter(
                pruned_arcs[:, 0],
                pruned_arcs[:, 1],
                facecolor=colors['dodgerblue'],#colors['aqua'], #if kept_bg[:, 0]%2==0 else colors['dodgerblue'],
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.3,
            )
        if arc_class is not None:
            plt.scatter(
                arc_class[:, 0],
                arc_class[:, 1],
                facecolor=orange,
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.1,
            )
        if bg_class is not None:
            plt.scatter(
                bg_class[:, 0],
                bg_class[:, 1],
                facecolor=colors['darkviolet'], #if kept_bg[:, 0]%2==0 else colors['mediumorchid'],
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.1,
            )
        """
        if kept_arcs is not None:
            plt.scatter(
                kept_arcs[:, 0],
                kept_arcs[:, 1],
                facecolor=orange,
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.1,
            )
        if kept_bg is not None:
            plt.scatter(
                kept_bg[:, 0],
                kept_bg[:, 1],
                facecolor=colors['darkviolet'], #if kept_bg[:, 0]%2==0 else colors['mediumorchid'],
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.1,
            )
        if switched_to_arc is not None:
            plt.scatter(
                switched_to_arc[:, 0],
                switched_to_arc[:, 1],
                facecolor=colors['forestgreen'], #if switched_to_arc[:, 0]%2==0 else colors['lime'],
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.1,
            )
        """
        
        if val_set is not None:
            plt.scatter(
                val_set[:, 0],
                val_set[:, 1],
                facecolor = colors['crimson'], #if switched_to_bg[:, 0]%2==0 else colors['magenta'],
                edgecolor="none",
                s=1,
                marker=",",
                zorder=2,
                alpha=0.1,
            )
        
        
        unseen = mpatches.Patch(color=colors['dodgerblue'], label='unseen/pruned')
        arc_arc = mpatches.Patch(color=orange, label='correct')
        #to_bg = mpatches.Patch(color=colors['forestgreen'], label='bg to arc')
        #to_arc= mpatches.Patch(color=colors['crimson'], label='arc to bg')
        bg_bg= mpatches.Patch(color=colors['darkviolet'], label='incorrect')
        vals = mpatches.Patch(color=colors['crimson'], label='val. set')
        handles=[unseen,arc_arc,bg_bg, vals]
        plt.legend(handles=[unseen,arc_arc,bg_bg, vals], bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=3, borderaxespad=0.) #bbox_to_anchor=(0.5,-0.1))#, mode="expand", borderaxespad=0.)  bbox_to_anchor=(0., 1.02, 1., .102),# mode="expand",
        plt.show(block=True)

    def construct_primal_graph(self):
        """ Constructs the primal graph where Critical Points become
            nodes and Morse-Smale arcs are represented as edges.
        """
        if self.primal is None:
            self.primal = nx.Graph()
            for arc in self.arcs:
                self.primal.add_edge(arc.node_ids[0], arc.node_ids[1])
        return self.primal.copy()

    def construct_dual_graph(self):
        """ Constructs a dual graph where Morse-Smale Arcs become nodes
            and edges represent when two Arcs share a critical point or
            otherwise intersect.
        """
        if self.dual is None:
            node_map = {}
            G = nx.Graph()
            current_arc_id = 0
            for arc in self.arcs:
                for node_id in arc.node_ids:
                    if node_id not in node_map:
                        node_map[node_id] = []
                    node_map[node_id].append(current_arc_id)
                G.add_node(
                    current_arc_id, index=arc.node_ids, size=len(arc.line)
                )
                current_arc_id += 1
                
            for arc_id, arc in list(G.nodes_iter(data=True)):
                
                for node_id in arc["index"]:
                    for connected_arc_id in node_map[node_id]:
                        if arc_id != connected_arc_id:
                            G.add_edge(arc_id, connected_arc_id)
            self.dual = G
        return self.dual.copy()

    def create_graphsage_input(self, in_arcs, out_arcs, gt_in_arcs, gt_out_arcs,label_scheme = 'train', val_count=20, val_in_arcs=None, val_out_arcs=None):
        """ Creates the necessary input data structures for GraphSage's
            GNN implementation.

        Keyword arguments:
            in_arcs -- a set of arcs representing "foreground" selection
            background_indices -- a set of arcs representing
                                  "background" selection

        Returns:
            tuple: (string -- A networkx-specified json representation
                              describing the input graph. Nodes have
                              'val' and 'test' attributes specifying if
                              they are a part of the validation and test
                              sets, respectively.
                    string -- A json-stored dictionary mapping the graph
                              node ids to consecutive integers.
                    string -- A json-stored dictionary mapping the graph
                              node ids to classes.
                    ndarray -- An array of node features; ordering given
                               by second item. Can be omitted and only
                               identity features will be used.)
        """
        if not gt_in_arcs:
            gt_in_arcs = in_arcs
            gt_out_arcs = out_arcs
        
        selected_indices = set()
        background_indices = set()
        gt_pos_indices = set()
        gt_neg_indices = set()
        val_pos_indices = set()
        val_neg_indices = set()

        for arc in in_arcs:
            selected_indices.add(arc)
        for arc in out_arcs:
            background_indices.add(arc)
        cvt_indices = in_arcs.union( out_arcs )

        for arc in gt_in_arcs:
            gt_pos_indices.add(arc)
        for arc in gt_out_arcs:
            gt_neg_indices.add(arc)

        for arc in val_in_arcs:
            val_pos_indices.add(arc)
        for arc in val_out_arcs:
            val_neg_indices.add(arc)
        val_set = val_pos_indices.union( val_neg_indices )
            
            
        compiled_data = self.compile_features()

        node_map = {}
        G = self.construct_dual_graph()

        current_arc_id = 0
        node_ids = {}
        node_labels = {}

        i_val = 0 ##!!!! need to set size for val set
        
        for arc, features in zip(self.arcs, compiled_data):
            index = tuple(arc.node_ids) + (len(arc.line),)
            for node_id in arc.node_ids:
                #
                # !! need to fix here to accomadate only ascending ridges
                #
                if node_id not in node_map:
                    node_map[node_id] = []
                node_map[node_id].append(current_arc_id)

                label = [
                    int(index in gt_neg_indices),
                    int(index in gt_pos_indices),
                ]

                node = G.node[current_arc_id]
                node["index"] = arc.node_ids
                node["size"] = len(arc.line)
                node["features"] = features.tolist()
                
                # labeled nodes assigned as train, test, or val
                if bool(np.sum(label)):
                    modified = 0
                    if index in cvt_indices:# and  i_val < val_count:
                        modified=1
                        node["train"] = True
                        node["test"] = False
                        node["val"] = False
                    if index in val_set:
                        modified=1
                        node["train"] = False
                        node["test"] = False
                        node["val"] = True
                    elif not bool(modified):
                        node["test"] = True
                        node["val"] = False
                        node["train"] = False
                    
                """Label all non-selected arcs as test"""
                #if not  bool(np.sum(label)):
                    #node["test"] = True
                    
                if bool(np.sum(label)):
                    node["label"] = label
                    node["prediction"] = None
                #G.node[current_arc_id] = node
                
                #current_arc_id += 1
                node_ids[current_arc_id] = current_arc_id
                node_labels[current_arc_id] = label
            current_arc_id += 1
        
        for arc_id, arc in list(G.nodes_iter(data=True)):#G.nodes.items():
            for node_id in arc["index"]:
                for connected_arc_id in node_map[node_id]:
                    G.add_edge(arc_id, connected_arc_id)

        data1 = json_graph.node_link_data(G) 
        s1 = json.dumps(data1) # input graph
        s2 = json.dumps(node_ids) #  dict: nodes to ints
        s3 = json.dumps(node_labels) #  dict: node_id to class

        return  (data1, node_ids, node_labels, compiled_data)#(s1, s2, s3, compiled_data)

