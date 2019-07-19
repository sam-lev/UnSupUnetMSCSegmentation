# Standard library imports

# Third party imports
import numpy as np
from sklearn import linear_model, model_selection
from sklearn.preprocessing import scale
from topoml.ui.vis import show_pixel_classifier_result

# Local application imports
from topoml.topology.utils import get_points_from_arcs
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
from topoml.NeuronTracer import NeuronTracer


class PixelNeuronTracer(NeuronTracer):
    def _add_default_features(self):
        self.features["identity"] = self.image
        self.features["sobel"] = sobel_filter(self.image)
        self.features["laplacian"] = laplacian_filter(self.image)
        for i in range(1, 5):
            pow1 = 2 ** (i - 1)
            pow2 = 2 ** i
            self.features["mean_{}".format(i)] = mean_filter(self.image, i)
            self.features["variance_{}".format(i)] = variance_filter(
                self.image, i
            )
            self.features["median_{}".format(i)] = median_filter(self.image, i)
            self.features["min_{}".format(i)] = minimum_filter(self.image, i)
            self.features["max_{}".format(i)] = maximum_filter(self.image, i)
            self.features["gauss_{}".format(pow2)] = gaussian_blur_filter(
                self.image, pow2
            )
            self.features[
                "delta_gauss_{}_{}".format(pow1, pow2)
            ] = difference_of_gaussians_filter(self.image, pow1, pow2)
        for i, neighbor_image in enumerate(neighbor_filter(self.image, 3)):
            self.features["shift_{}".format(i)] = neighbor_image

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

    def compile_features(self, selection=None, return_labels=False):
        if self.compiled_features is None:
            compiled_data = np.dstack(self.features.values())
            feature_names = list(self.features.keys())
            mu = np.mean(compiled_data, axis=0)
            std = np.std(compiled_data, axis=0)
            compiled_data = (compiled_data - mu) / std
            self.compiled_features = compiled_data
            self.feature_names = feature_names

        if selection is None:
            if return_labels:
                return self.compiled_features, self.feature_names
            else:
                return self.compiled_features
        if return_labels:
            return (
                self.compiled_features[selection[:, 1], selection[:, 0], :],
                self.feature_names,
            )
        else:
            return self.compiled_features[selection[:, 1], selection[:, 0], :]

    def _extract_selection(self, in_arcs, out_arcs, out_points):
        in_points = get_points_from_arcs(in_arcs)
        arc_out_points = get_points_from_arcs(out_arcs)
        if len(out_points) and len(arc_out_points):
            out_points = np.vstack((arc_out_points, out_points))
        elif len(arc_out_points):
            out_points = arc_out_points
        return in_points, out_points

    def show_classifier_result(self, classifier):
        compiled_data = self.compile_features()
        show_pixel_classifier_result(compiled_data, classifier)
