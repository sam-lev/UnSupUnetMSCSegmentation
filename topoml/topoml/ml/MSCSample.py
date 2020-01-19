#3rd party imports
import numpy as np
import networkx as nx

#local imports
from topoml.topology.utils import (
    get_pixel_values_from_arcs,
    is_ridge_arc,
    is_valley_arc,
)

class SupervisedMSC:
    def __init__(self, msc=None, geomsc=None, labeled_segmentation=None, ridge=True, valley=False):
        self.msc = msc
        self.geomsc = geomsc
        self.arcs = None
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
            self.labeled_msc = msc
        if geomsc is not None:
            self.geomsc = geomsc
            self.msc = geomsc
            self.arc_accuracies = self.geomsc.arc_accuracies
            self.labeled_msc = self.geomsc
        if labeled_segmentation is not None:
            self.labeled_segmentation = labeled_segmentation

        # This ensures smaller arcs take precedence in the event of
        # pixel overlap
        sorted_arcs = sorted(self.labeled_msc.arcs, key=lambda arc: len(arc.line))
        arc_points = []
        self.arcs = []
        """for arc in sorted_arcs:
            if (not self.use_ridge_arcs and is_ridge_arc(arc, self.labeled_msc)) or (
                    not self.use_valley_arcs and is_valley_arc(arc, self.labeled_msc)
            ):
                continue
            index = tuple(arc.node_ids) + (len(arc.line),)#make_arc_id(arc)
            arc_points.extend(arc.line)
            self.arc_map.extend([index] * len(arc.line))
            self.arcs.append(arc)"""

        for arc in self.labeled_msc.arcs:
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

    def compile_features(self, selection=None, return_labels=False, images = {}):
        # G = self.construct_dual_graph()
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
