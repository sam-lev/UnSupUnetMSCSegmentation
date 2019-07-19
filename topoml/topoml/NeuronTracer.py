# Standard library imports
from abc import ABC, abstractmethod

# Third party imports
import numpy as np
from skimage import io

# Local application imports
from topoml.topology.utils import build_msc
from topoml.ui.ArcSelector import ArcSelector
from topoml.image.utils import blur_and_save


class NeuronTracer(ABC):
    def __init__(self, fname=None, blur_sigma=2, persistence=0, construct_msc = True):
        self.compiled_features = None
        self.feature_names = None

        

        if construct_msc:
            self.raw_image = io.imread(fname, as_gray=True)
            self.image, fname_raw = blur_and_save(
                self.raw_image, fname.rsplit(".", 1)[0], blur_sigma
            )
        else:
            self.raw_image = None
            self.image = None
            

        self.features = {}
        self.models = {}

        self.foreground_arcs = set()
        self.background_arcs = set()

        if construct_msc:
            self.msc = build_msc(
                fname_raw, self.image.shape[1], self.image.shape[0], persistence
            )

            self._add_default_features()
            self._add_default_models()
        else:
            self.msc = None

        self.positive_arcs = []
        self.negative_arcs = []

    @abstractmethod
    def _add_default_features(self):
        pass

    @abstractmethod
    def _add_default_models(self):
        pass

    @abstractmethod
    def compile_features(self, selection, return_labels=False):
        pass

    @abstractmethod
    def _extract_selection(self, in_arcs, out_arcs, out_pixels):
        pass

    def _build_models(self, in_selection, out_selection):
        # Models can accept one of two inputs: either a list of in/out
        # pixels stacked with different features, or a collection of
        # in/out arcs stacked with arc-specific features, each model
        # will get both inputs and will have to parse out what they
        # want
        in_data = self.compile_features(in_selection)
        out_data = self.compile_features(out_selection)

        features = np.vstack((in_data, out_data))
        labels = np.hstack((np.ones(len(in_data)), np.zeros(len(out_data))))
        for model in self.models.values():
            model.fit(features, labels)

    def prune_msc(self, classifier_name):
        return self.models[classifier_name].predict(self.image, self.msc)

    def add_feature(self, name, feature):
        self.features[name] = feature
        # We have to remove any cached data, since we are modifying its
        # structure
        self.compiled_features = None
        self.feature_names = None

    def add_classifier(self, name, model):
        self.models[name] = model

    def do_selection(self, selector=None, xlims=None, ylims=None, filename=None):
        if selector is None:
            selector = ArcSelector(self.raw_image, self.msc)

        if filename is not None:
            in_arcs, out_arcs, out_pixels = selector.write_image(filename)
        else:
            in_arcs, out_arcs, out_pixels = selector.launch_ui(xlims, ylims)
        self.positive_arcs = in_arcs
        self.negative_arcs = out_arcs

        in_selection, out_selection = self._extract_selection(
            in_arcs, out_arcs, out_pixels
        )
        self._build_models(in_selection, out_selection)
        return (in_arcs, out_arcs, out_pixels)
