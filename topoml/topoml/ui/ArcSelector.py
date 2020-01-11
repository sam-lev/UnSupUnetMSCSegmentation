# Standard library imports
import sys
import os

# Third party imports
import matplotlib.pyplot as plt
import matplotlib.colors as pltcolor
from skimage import io
import imageio
from PIL import Image

import numpy as np
import scipy
from skimage import morphology
import samply

# Local application imports
from topoml.ui.colors import red, blue, green, orange, purple
orange = red
from topoml.image.utils import (
    bounding_box,
    box_intersection,
    make_dilated_arc_image,
    make_arc_mask,
    make_mc_arc_mask,
)
from topoml.topology.utils import is_ridge_arc, is_valley_arc


def make_arc_id(a):
    return tuple(a.node_ids) + (len(a.line),)

class ArcSelector(object):
    def __init__(
            self, image, msc, selection_radius=10, valley=True, ridge=True, invert=False, kdtree=False
    ):
        # Needed for kdTree on large point sets:
        sys.setrecursionlimit(10000)

        self.image = image
        self.msc = msc
        self.msc.arcs = msc.arcs

        self.invert = invert
        
        self.selection_radius = selection_radius
        self.selection_shape = morphology.disk(self.selection_radius)
        self.fat_mask = make_dilated_arc_image(
            self.image, self.msc, self.selection_radius, self.invert
        )

        self.in_color = orange
        self.out_color = purple
        self.in_arcs = set()
        self.out_arcs = set()
        self.out_pixels = set()
        self.arc_map = []
        self.arc_drawings = {}
        self.use_valley_arcs = valley
        self.use_ridge_arcs = ridge

        # This ensures smaller arcs take precedence in the event of
        # pixel overlap
        sorted_arcs = sorted(self.msc.arcs, key=lambda arc: len(arc.line))
        arc_points = []
        self.arcs = []
        for arc in sorted_arcs:
            if (not self.use_ridge_arcs and is_ridge_arc(arc, self.msc)) or (
                not self.use_valley_arcs and is_valley_arc(arc, self.msc)
            ):
                continue
            index = make_arc_id(arc)
            arc_points.extend(arc.line)
            self.arc_map.extend([index] * len(arc.line))
            self.arcs.append(arc)
        # only needed for selection ui to choose neighboring arcs
        # can cause error with sparse MSC
        self.kdtree = None
        if kdtree:
            self.kdtree = scipy.spatial.KDTree(arc_points, leafsize=1000)

    def get_closest_arc_index(self, point):
        distance, index = self.kdtree.query(point)
        return self.arc_map[index]

    def launch_ui(self, xlims=None, ylims=None):
        plt.ion()
        self.fig = plt.figure()
        plt.imshow(self.image, cmap=plt.cm.Greys_r, zorder=1)

        for arc in self.msc.arcs:
            if (not self.use_ridge_arcs and is_ridge_arc(arc, self.msc)) or (
                not self.use_valley_arcs and is_valley_arc(arc, self.msc)
            ):
                continue

            arc_index = make_arc_id(arc)
            points = np.array(arc.line)
            self.arc_drawings[arc_index] = plt.scatter(
                points[:, 0],
                points[:, 1],
                facecolor="none",
                edgecolor="none",
                s=2,
                marker=",",
                zorder=3,
            )
            if arc_index in self.in_arcs:
                self.arc_drawings[arc_index].set_facecolor(self.in_color)
                self.arc_drawings[arc_index].set_visible(True)
            elif arc_index in self.out_arcs:
                self.arc_drawings[arc_index].set_facecolor(self.out_color)
                self.arc_drawings[arc_index].set_visible(True)
            else:
                self.arc_drawings[arc_index].set_visible(False)

        if self.use_ridge_arcs:
            arc_mask = make_mc_arc_mask(self.image, self.msc, False)
            plt.imshow(
                arc_mask,
                cmap="Oranges",
                vmin=0,
                vmax=4,
                interpolation="none",
                alpha=0.8,
                zorder=2,
            )
        if self.use_valley_arcs:
            arc_mask = make_mc_arc_mask(self.image, self.msc, True)
            plt.imshow(
                arc_mask,
                cmap="Blues",
                vmin=0,
                vmax=4,
                interpolation="none",
                alpha=0.8,
                zorder=2,
            )

        extrema_points = [[], [], []]
        for node in self.msc.nodes.values():
            x, y = node.xy
            extrema_points[node.index].append([x, y])

        for i, color in enumerate([blue, green, red]):
            xy = np.array(extrema_points[i])
            plt.scatter(
                xy[:, 0],
                xy[:, 1],
                facecolor=color,
                edgecolor="none",
                s=1,
                marker=",",
                zorder=4,
            )

        plt.gca().set_xlim(0, self.image.shape[1])
        plt.gca().set_ylim(self.image.shape[0], 0)
        if xlims is not None:
            plt.gca().set_xlim(xlims[0], xlims[1])
        if ylims is not None:
            plt.gca().set_ylim(ylims[1], ylims[0])

        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        plt.show(block=True)

        in_arcs = []
        out_arcs = []
        for arc in self.msc.arcs:
            index = make_arc_id(arc)
            if index in self.in_arcs:
                in_arcs.append(arc)
            elif index in self.out_arcs:
                out_arcs.append(arc)

        return (in_arcs, out_arcs, np.array(list(self.out_pixels)))

    def write_image(self, filename):
        self.fig = plt.figure()
        plt.imshow(self.image, cmap=plt.cm.Greys_r, zorder=1)

        for arc in self.msc.arcs:
            if (not self.use_ridge_arcs and is_ridge_arc(arc, self.msc)) or (
                not self.use_valley_arcs and is_valley_arc(arc, self.msc)
            ):
                continue

            arc_index = make_arc_id(arc)
            if arc_index in self.in_arcs:
                points = np.array(arc.line)
                self.arc_drawings[arc_index] = plt.scatter(
                    points[:, 0],
                    points[:, 1],
                    facecolor= self.in_color,
                    edgecolor="none",
                    s=2,
                    marker=",",
                    zorder=3,
                )
                self.arc_drawings[arc_index].set_visible(True)
            elif arc_index in self.out_arcs:
                points = np.array(arc.line)
                self.arc_drawings[arc_index] = plt.scatter(
                    points[:, 0],
                    points[:, 1],
                    facecolor=self.out_color,
                    edgecolor="none",
                    s=2,
                    marker=",",
                    zorder=3,
                )
                self.arc_drawings[arc_index].set_visible(True)

        if self.use_ridge_arcs:
            arc_mask = make_mc_arc_mask(self.image, self.msc, False)
            plt.imshow(
                arc_mask,
                cmap="Oranges",
                vmin=0,
                vmax=4,
                interpolation="none",
                alpha=0.8,
                zorder=2,
            )
        if self.use_valley_arcs:
            arc_mask = make_mc_arc_mask(self.image, self.msc, True)
            plt.imshow(
                arc_mask,
                cmap="Blues",
                vmin=0,
                vmax=4,
                interpolation="none",
                alpha=0.8,
                zorder=2,
            )

        extrema_points = [[], [], []]
        for node in self.msc.nodes.values():
            x, y = node.xy
            extrema_points[node.index].append([x, y])

        for i, color in enumerate([blue, green, red]):
            xy = np.array(extrema_points[i])
            plt.scatter(
                xy[:, 0],
                xy[:, 1],
                facecolor=color,
                edgecolor="none",
                s=1,
                marker=",",
                zorder=4,
            )

        plt.gca().set_xlim(0, self.image.shape[1])
        plt.gca().set_ylim(self.image.shape[0], 0)
        plt.savefig(filename)

        in_arcs = []
        out_arcs = []
        for arc in self.msc.arcs:
            index = make_arc_id(arc)
            if index in self.in_arcs:
                in_arcs.append(arc)
            elif index in self.out_arcs:
                out_arcs.append(arc)

        return (in_arcs, out_arcs, np.array(list(self.out_pixels)))
    
    def make_to_scale_image(image_in, outputname, size=(1, 1), dpi=80):
        fig = plt.figure()
        fig.set_size_inches(size)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        plt.set_cmap('binary')
        ax.imshow(image_in, aspect='equal')
        plt.savefig(outputname, dpi=dpi)

    def draw_binary_segmentation(self, filename, msc = None, invert=False, reshape_out = True, dpi = True):
        if not msc:
            self.msc = msc
        
        black_box = np.zeros((self.image.shape[0],self.image.shape[1])) if not invert else np.zeros((self.image.shape[1],self.image.shape[0]))
        #print(self.image.shape+(,,3))
        fig = plt.imshow(black_box,cmap='gray', alpha=None,interpolation = 'nearest') #plt.figure() #in
        plt.axis('off')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        for arc in self.msc.arcs:
            if (not self.use_ridge_arcs and is_ridge_arc(arc, self.msc)) or (
                not self.use_valley_arcs and is_valley_arc(arc, self.msc)
            ):
                continue

            arc_index = make_arc_id(arc)
            if arc_index in self.in_arcs:
                points = np.array(arc.line)
                arc_drawings[arc_index] = plt.scatter(
                    points[:, 0],
                    points[:, 1],
                    facecolor= 'white',#self.in_color,
                    alpha=None,
                    edgecolor="none",
                    s=5,
                    marker=",",
                    zorder=3,
                )
                arc_drawings[arc_index].set_visible(True)
            elif arc_index in self.out_arcs:
                points = np.array(arc.line)
                arc_drawings[arc_index] = plt.scatter(
                    points[:, 0],
                    points[:, 1],
                    facecolor= 'white',#self.out_color,
                    edgecolor="none",
                    alpha=None,
                    s=5,
                    marker=",",
                    zorder=3,
                )
                arc_drawings[arc_index].set_visible(True)

        if self.use_ridge_arcs:
            arc_mask = make_mc_arc_mask(self.image, self.msc, invert)#False)
            plt.imshow(
                arc_mask,
                cmap="binary",
                vmin=0,
                vmax=0,
                norm=plt.Normalize(vmin=0, vmax=1),
                interpolation="nearest",
                alpha=None,
                zorder=2,
            )
        if self.use_valley_arcs:
            arc_mask = make_mc_arc_mask(self.image, self.msc, invert)#True)
            plt.imshow(
                arc_mask,
                cmap="binary",
                vmin=0,
                vmax=0,
                interpolation="nearest",
                alpha=None,
                zorder=2,
            )
        """
        extrema_points = [[], [], []]
        for node in self.msc.nodes.values():
            x, y = node.xy
            extrema_points[node.index].append([x, y])

        for i, color in enumerate([blue, green, red]):
            xy = np.array(extrema_points[i])
            plt.scatter(
                xy[:, 0],
                xy[:, 1],
                facecolor=color,
                edgecolor="none",
                s=1,
                marker=",",
                zorder=4,
            )
        """

        #plt.gca().set_axis_off()
        #plt.gca().set_xlim(0, self.image.shape[0])
        #plt.gca().set_ylim(self.image.shape[1], 0)
        dpi_ = None
        if dpi:
            if self.image.shape[0] >= 600:
                dpi_ = 600
            else:
                dpi_ = 156
            if isinstance(dpi, int):
                dpi_ = dpi
        
        plt.savefig(filename, bbox_inches='tight', pad_inches = 0.0, dpi = dpi_, transparent = False, cmap=None)

        img = imageio.imread(filename)
        
        if reshape_out:
            if self.image.shape[0] >= 600:
                img = Image.fromarray(img).resize((img.shape[1]+9, img.shape[1]-91+5))
            else:
                img = Image.fromarray(img).resize((img.shape[1]+1, img.shape[0]+1))
        else:
            img = Image.fromarray(img).resize((img.shape[1]+1, img.shape[0]+1)) #if not invert else Image.fromarray(img).resize((img.shape[0], img.shape[1]))
        img = np.array(img)[:,:,0]
        Img = Image.fromarray(img)
        Img.save(filename)

        plt.close()
        
        in_arcs = []
        out_arcs = []
        for arc in self.msc.arcs:
            index = make_arc_id(arc)
            if index in self.in_arcs:
                in_arcs.append(arc)
            elif index in self.out_arcs:
                out_arcs.append(arc)

        return (in_arcs, out_arcs, np.array(list(self.out_pixels)))

    def toggle_arc(self, x, y, in_class=True):
        pt = np.array([x, y])
        min_index = self.get_closest_arc_index(pt)
        if in_class:
            if min_index in self.in_arcs:
                self.in_arcs.remove(min_index)
            elif min_index in self.out_arcs:
                self.out_arcs.remove(min_index)
            else:
                self.in_arcs.add(min_index)
        else:
            if min_index in self.out_arcs:
                self.out_arcs.remove(min_index)
            elif min_index in self.in_arcs:
                self.in_arcs.remove(min_index)
            else:
                self.out_arcs.add(min_index)
        return min_index

    def highlight_pixels(self, x, y):
        start_y = y - self.selection_radius
        start_x = x - self.selection_radius
        for i in range(0, 2 * self.selection_radius + 1):
            xi = start_x + i
            if xi > 0 and xi < self.fat_mask.shape[1]:
                for j in range(0, 2 * self.selection_radius + 1):
                    yj = start_y + j
                    if yj > 0 and yj < self.fat_mask.shape[0]:
                        if (
                            self.selection_shape[j, i]
                            and not self.fat_mask[yj, xi]
                        ):
                            self.out_pixels.add((xi, yj))

    def on_click(self, event):
        if event is None or None in [event.xdata, event.ydata, event.button]:
            return

        if self.fat_mask[int(event.ydata), int(event.xdata)]:
            selected_index = self.toggle_arc(
                event.xdata, event.ydata, event.button == 1
            )
            self.arc_drawings[selected_index].set_visible(
                not self.arc_drawings[selected_index].get_visible()
            )
            if event.button == 1:
                self.arc_drawings[selected_index].set_facecolor(self.in_color)
            else:
                self.arc_drawings[selected_index].set_facecolor(self.out_color)
        else:
            self.highlight_pixels(int(event.xdata), int(event.ydata))

        if len(self.out_pixels):
            bg_points = np.array(list(self.out_pixels))
            plt.scatter(
                bg_points[:, 0],
                bg_points[:, 1],
                edgecolor="none",
                facecolor=self.out_color,
                s=1,
                marker=",",
                zorder=4,
            )
        self.fig.canvas.draw()

    def save_arcs(self, filename="arcs.csv", mode="a"):
        f = open(filename, mode)
        for index in self.in_arcs:
            f.write("{},{},{},{}\n".format(1, *index))
        for index in self.out_arcs:
            f.write("{},{},{},{}\n".format(0, *index))
        f.close()

    def load_arcs(self, filename="arcs.csv"):
        if os.path.exists(filename):
            f = open(filename, "r")
            for line in f:
                tokens = list(map(int, line.strip().split(",")))
                if tokens[0]:
                    self.in_arcs.add(tuple(tokens[1:]))
                else:
                    self.out_arcs.add(tuple(tokens[1:]))
        return (self.in_arcs, self.out_arcs)

    def cull_selection(self, window):
        for a in self.msc.arcs:
            arc_index = make_arc_id(a)
            for arc_set in [self.in_arcs, self.out_arcs]:
                if arc_index in arc_set:
                    if not box_intersection(window, bounding_box(a.line)):
                        arc_set.remove(arc_index)
        return self.in_arcs, self.out_arcs

    def construct_node_map(self):
        node_map = {}
        current_arc_index = 0
        for arc in self.arcs:
            for node_id in arc.node_ids:
                if node_id not in node_map:
                    node_map[node_id] = []
                node_map[node_id].append(current_arc_index)
            current_arc_index += 1
        return node_map

    def sample_selection(self, count=20, rings=1, seed=0):
        np.random.seed(seed)
        X = samply.hypercube.cvt(count, 2)
        X[:, 0] *= self.image.shape[1]
        X[:, 1] *= self.image.shape[0]

        node_map = self.construct_node_map()
        
        seed_arc_keys = list()

        for x in X:
            arc_key = self.get_closest_arc_index(x)
            seed_arc_keys.append(arc_key)
        ring = 0
        ring_index = 0
        ring_count = len(seed_arc_keys)

        while ring <= rings:
            next_ring = seed_arc_keys[ring_index:(ring_index+ring_count)]

            ring_count = 0
            for arc_key in next_ring:
                for node_id in arc_key[:2]:
                    for arc_index in node_map[node_id]:
                        neighbor_arc_key = make_arc_id(self.arcs[arc_index])
                        if neighbor_arc_key not in seed_arc_keys:
                            seed_arc_keys.append(neighbor_arc_key)
                            ring_count += 1
                ring_index += 1
            ring += 1

        seed_arc_keys = set(seed_arc_keys)
        
        """
        for arc_set in [self.in_arcs, self.out_arcs]:
            for arc_key in self.arc_map:
                if arc_key in arc_set and arc_key not in seed_arc_keys:
                    arc_set.remove(arc_key)
        """
        return (self.in_arcs.intersection(seed_arc_keys), self.out_arcs.intersection(seed_arc_keys))
