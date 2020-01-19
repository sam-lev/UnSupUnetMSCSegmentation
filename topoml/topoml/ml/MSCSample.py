#3rd party imports
import numpy as np
import networkx as nx
import scipy.stats
import json
from networkx.readwrite import json_graph

#local imports
from topoml.ui.ArcSelector import ArcSelector
from topoml.topology.utils import (
    get_pixel_values_from_arcs,
    is_ridge_arc,
    is_valley_arc,
)
from topoml.ml.utils import gaussian_fit
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

class MSCSample():
    def __init__(self, msc=None, geomsc=None, labeled_segmentation=None, ridge=True, valley=False):

        self.msc = msc
        self.geomsc = geomsc
        self.arcs = None
        self.ridge = ridge
        self.valley= valley
        self.arc_map = None
        self.arc_accuracies = None
        self.labeled_segmentation = None
        # msc after labeling arcs (updating arc.label_accuracy
        self.labeled_msc = None
        self.image = {}
        self.features = {}

    def msc_arc_accuracy(self, arc=None, msc=None, geomsc=None, labeled_segmentation=None, invert=True):
        if labeled_segmentation is not None:
            self.labeled_segmentation = labeled_segmentation
        arc_accuracy = 0
        for point in arc.line:
            x = point[0]
            y = point[1]
            if invert:
                x = point[1]
                y = point[0]
            if self.labeled_segmentation[int(x),int(y)] > 0:
                arc_accuracy += 1.
        label_accuracy = arc_accuracy/float(len(arc.line))
        if label_accuracy == 0.:
            label_accuracy = 1e-6
        return label_accuracy



    def map_labeling(self, image=None, msc=None, geomsc=None, labeled_segmentation=None, invert=True):
        if msc is not None:
            self.msc = msc
        if geomsc is not None:
            self.geomsc = geomsc
            self.msc = geomsc
            self.arc_accuracies = self.geomsc.arc_accuracies
        if labeled_segmentation is not None:
            self.labeled_segmentation = labeled_segmentation

        # This ensures smaller arcs take precedence in the event of
        # pixel overlap
        sorted_arcs = sorted(self.msc.arcs, key=lambda arc: len(arc.line))
        arc_points = []
        self.arcs = []
        for arc in sorted_arcs:
            if (not self.ridge and is_ridge_arc(arc, self.msc)) or (
                    not self.valley and is_valley_arc(arc, self.msc)
            ):
                continue
            index = tuple(arc.node_ids) + (len(arc.line),)#make_arc_id(arc)
            arc_points.extend(arc.line)
            self.arc_map.extend([index] * len(arc.line))
            self.arcs.append(arc)

        for arc in self.msc.arcs:
            arc.label_accuracy = self.msc_arc_accuracy(arc=arc
                                                       , labeled_segmentation=self.labeled_segmentation,
                                                       invert=invert)
        return self.labeled_msc


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

    def set_default_features(self, images=None):
        if images is not None:
            self.images = images
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

    def compile_features(self, selection=None, return_labels=False, images=None):
        if images is not None:
            self.images = images

        if not self.features:
            self.set_default_features()

        G = self.construct_dual_graph()
        arc_map = {}
        arc_pixel_map = {}
        for i, arc in enumerate(self.arcs):
            index = tuple(arc.node_ids) + (len(arc.line),)
            arc_map[index] = i
        # G = self.construct_dual_graph()
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

    def partition_msc_edges_train_test_val(self, number_samples, hops, msc=None):
        if msc is not None:
            self.msc = msc

        # collect/compute features before partition
        compiled_data = self.compile_features()

        self.positive_arcs = set()
        self.negative_arcs = set()
        self.validation_set = {}
        self.training_set = {}
        self.test_set = {}

        def fill_set(list):
            new_set = set()
            for s in list:
                new_set.add(s)
            return new_set

        for arc in self.msc.arcs:
            if arc.label_accuracy == 1.:
                self.positive_arcs.add(arc)
            else:
                self.negative_arcs.add(arc)
            # could add class for high/mid accuracy arcs

        # use cvt sampling to obtain validation/test edges
        edge_selector = ArcSelector(msc=self.msc, valley=False)
        in_arcs, out_arcs = edge_selector.sample_selection(count=number_samples, rings=hops, seed=123)
        val_in_arcs, val_out_arcs = edge_selector.sample_selection(count=1, rings=int(hops / 2), seed=111)
        test_in_arcs, test_out_arcs = edge_selector.sample_selection(count=1, rings=int(hops / 2), seed=323)

        # set validation and test test arcs from cvt sampling
        self.validation_set["positive"] = fill_set(val_in_arcs)
        self.validation_set["negative"] = fill_set(val_out_arcs)
        self.test_set["positive"] = fill_set(test_in_arcs)
        self.test_set["negative"] = fill_set(test_out_arcs)

        val_and_test = self.validation_set["positive"].union(self.validation_set["negative"]).union(self.test_set["positive"]).union(self.test_set["negative"])

        node_map = {}
        G = self.construct_dual_graph()

        current_arc_id = 0
        node_ids = {}
        node_labels = {}

        i_val = 0  ##!!!! need to set size for val set

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
                    int(index in self.negative_arcs),
                    int(index in self.positive_arcs),
                ]

                node = G.node[current_arc_id]
                node["index"] = arc.node_ids
                node["size"] = len(arc.line)
                node["features"] = features.tolist()

                # labeled nodes assigned as train, test, or val
                if bool(np.sum(label)):
                    modified = 0
                    if index not in val_and_test:  # and  i_val < val_count:
                        modified = 1
                        node["train"] = True
                        node["test"] = False
                        node["val"] = False
                    if index in self.validation_set["positive"].union(self.validation_set["negative"]):
                        modified = 1
                        node["train"] = False
                        node["test"] = False
                        node["val"] = True
                    elif index in self.test_set["positive"].union(self.test_set["negative"]):
                        node["test"] = True
                        node["val"] = False
                        node["train"] = False

                """Label all non-selected arcs as test"""
                # if not  bool(np.sum(label)):
                # node["test"] = True

                if bool(np.sum(label)):
                    node["label"] = label #arc.label_accuracy
                    node["label_accuracy"] = arc.label_accuracy
                    node["prediction"] = None
                # G.node[current_arc_id] = node

                # current_arc_id += 1
                node_ids[current_arc_id] = current_arc_id
                node_labels[current_arc_id] = label
            current_arc_id += 1

        for arc_id, arc in list(G.nodes_iter(data=True)):  # G.nodes.items():
            for node_id in arc["index"]:
                for connected_arc_id in node_map[node_id]:
                    G.add_edge(arc_id, connected_arc_id)

        data1 = json_graph.node_link_data(G)
        s1 = json.dumps(data1)  # input graph
        s2 = json.dumps(node_ids)  # dict: nodes to ints
        s3 = json.dumps(node_labels)  # dict: node_id to class

    def create_graphsage_input(self, in_arcs, out_arcs, gt_in_arcs, gt_out_arcs, label_scheme='train', val_count=20,
                               val_in_arcs=None, val_out_arcs=None):
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
        cvt_indices = in_arcs.union(out_arcs)

        for arc in gt_in_arcs:
            gt_pos_indices.add(arc)
        for arc in gt_out_arcs:
            gt_neg_indices.add(arc)

        for arc in val_in_arcs:
            val_pos_indices.add(arc)
        for arc in val_out_arcs:
            val_neg_indices.add(arc)
        val_set = val_pos_indices.union(val_neg_indices)

        compiled_data = self.compile_features()

        node_map = {}
        G = self.construct_dual_graph()

        current_arc_id = 0
        node_ids = {}
        node_labels = {}

        i_val = 0  ##!!!! need to set size for val set

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
                    if index in cvt_indices:  # and  i_val < val_count:
                        modified = 1
                        node["train"] = True
                        node["test"] = False
                        node["val"] = False
                    if index in val_set:
                        modified = 1
                        node["train"] = False
                        node["test"] = False
                        node["val"] = True
                    elif not bool(modified):
                        node["test"] = True
                        node["val"] = False
                        node["train"] = False

                """Label all non-selected arcs as test"""
                # if not  bool(np.sum(label)):
                # node["test"] = True

                if bool(np.sum(label)):
                    node["label"] = label
                    node["prediction"] = None
                # G.node[current_arc_id] = node

                # current_arc_id += 1
                node_ids[current_arc_id] = current_arc_id
                node_labels[current_arc_id] = label
            current_arc_id += 1

        for arc_id, arc in list(G.nodes_iter(data=True)):  # G.nodes.items():
            for node_id in arc["index"]:
                for connected_arc_id in node_map[node_id]:
                    G.add_edge(arc_id, connected_arc_id)

        data1 = json_graph.node_link_data(G)
        s1 = json.dumps(data1)  # input graph
        s2 = json.dumps(node_ids)  # dict: nodes to ints
        s3 = json.dumps(node_labels)  # dict: node_id to class

        return (data1, node_ids, node_labels, compiled_data)  # (s1, s2, s3, compiled_data)
