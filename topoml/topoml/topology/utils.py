# Standard library imports
import subprocess
import os

# Third party imports
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

# Local application imports
from topoml.image.utils import blur_and_save
from topoml.ml.utils import gaussian_fit
from topoml.topology.msc import MSC
from topoml.ui.vis import show_image
from topoml.ui.colors import blue


cached_msc_path = None


def find_test_2d_viewer(starting_path="/"):
    global cached_msc_path
    if cached_msc_path is None:
        name = "test_2d_viewer"
        for root, dirs, files in os.walk(starting_path):
            if name in files:
                cached_msc_path = os.path.join(root, name)
                break
    if cached_msc_path is None:
        raise FileNotFoundError()
    return cached_msc_path


def filtered_msc_training_set(filtered_images, kernel_size, persistence=8.5):
    pos_pixels = []
    neg_pixels = []
    gauss_msc_pos = []
    gauss_msc_neg = []

    for i in range(kernel_size):
        # Compute the morse smale complex of all blurred images and
        # filter arcs of interest based on mu from gaussian fit to pixel
        # histogram for each image
        filtered_smooth, fname = blur_and_save(
            filtered_images[:, :, i], "filtered_msc_training_set_{}".format(i)
        )
        full_msc = build_msc(
            fname,
            filtered_smooth.shape[1],
            filtered_smooth.shape[0],
            persistence,
        )
        mu, sigma, var = gaussian_fit(filtered_images[:, :, i])
        filtered_msc = msc_filtration(filtered_images[:, :, i], full_msc, mu)

        # Cover filtered MSC with closed balls to select negative space
        # as background and arcs as positive samples
        pos_samples, neg_samples = collect_training_from_filtered(
            filtered_images[:, :, i], filtered_msc
        )

        # a = np.array(neg_samples)

        # get pixel values under arc for positives
        for arc_pos in pos_samples:
            for point_pos in arc_pos:
                pos_pixels.append(
                    [filtered_images[int(point_pos[0]), int(point_pos[1]), i]]
                )

        # get negative pixel values
        for bg_group in neg_samples:
            for point_neg in bg_group:
                neg_pixels.append(
                    [filtered_images[int(point_neg[0]), int(point_neg[1]), i]]
                )

        # fit gaussian to negative pixel values in background and take
        # mean, standard deviation, variance
        for arc_neg in neg_samples:
            filtered_negatives = []
            for point_neg in arc_neg:
                filtered_negatives.append(
                    [filtered_images[int(point_neg[0]), int(point_neg[1]), i]]
                )
            # fit gaussian to negative sample pixels in background
            mu_neg, sigma_neg, var_neg = gaussian_fit(
                np.array(filtered_negatives)
            )
            gauss_msc_neg.append([mu_neg, sigma_neg, var_neg])

        # fit gaussian to positive pixel values under arcs of MSC and
        # take mean, standard deviation, and variance
        for arc in pos_samples:
            filtered_positives = []
            for point in arc:
                filtered_positives.append(
                    [filtered_images[int(point[0]), int(point[1]), i]]
                )
            # fit gaussian to positive pixel values under arc and take
            # mean, standard deviation, and variance
            # NOTE! can filter based on arc length, currently doing
            # three bc three degrees of freedom in gaussian fit.
            if len(filtered_positives) >= 3:
                mu_pos, sigma_pos, var_pos = gaussian_fit(
                    np.array(filtered_positives)
                )
                gauss_msc_pos.append([mu_pos, sigma_pos, var_pos])

    train_x = np.vstack((pos_pixels, neg_pixels))
    train_y = np.hstack(
        (
            np.ones(np.array(pos_pixels).shape[0]),
            np.zeros(np.array(neg_pixels).shape[0]),
        )
    )

    gaussian_msc_arc_X = np.vstack((gauss_msc_pos, gauss_msc_neg))
    gaussian_msc_arc_Y = np.hstack(
        (
            np.ones(np.array(gauss_msc_pos).shape[0]),
            np.zeros(np.array(gauss_msc_neg).shape[0]),
        )
    )

    return train_x, train_y, gaussian_msc_arc_X, gaussian_msc_arc_Y


def build_msc(raw_filename, width, height, persistence=0):
    """Calls Attila's C++ code using subprocess

       TODO:
       Dan: write a proper swig or python interface to Attila's code

    Keyword arguments:
        raw_filename -- a string representing the anme of the image file
        in raw format
        width -- The width of the image
        height --  The height of the image
        persistence -- a floatint point value specifying the persistence
        to use

    Returns:
        MSC: An msc object
    """
    # TODO: this works for me specifically because of where I installed things,
    # the default is to start from root which could take a potentially long time.
    # Hard-code your own path here for now if your libraries are not installed
    # side-by-side. 
    starting_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..",'..'
    )
    proc = subprocess.Popen(
        [   'sh',
            find_test_2d_viewer(starting_dir),
            raw_filename,
            str(width),
            str(height),
            "1",
            "1",
            str(persistence),
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = proc.communicate()
    # print(stdout.decode("utf-8"))
    # print(stderr.decode("utf-8"))
    msc = MSC()
    msc.read_from_file(raw_filename)
    return msc


def msc_filtration(image, msc, mu):
    # store arcs containing pixel values above mean of gaussian fit to
    # image histogram
    filtered_arcs = []
    for a in msc.arcs:
        points = np.array(a.line)
        arc_points = []
        for point in points:
            if image[int(point[1]), int(point[0])] > mu:
                arc_points.append(tuple((point[1], point[0])))
        if len(arc_points):
            filtered_arcs.append(arc_points)

    # make id dictionary for each arc
    filtered_min = []
    filtered_max = []
    filtered_saddle = []
    for node in msc.nodes.values():
        if node.xy in filtered_arcs:
            if node.index == 0:
                filtered_min.append(node.xy)
            if node.index == 1:
                filtered_saddle.append(node.xy)
            if node.index == 2:
                filtered_max.append(node.xy)

    # return filtered msc as [arcs, minima, saddle, maxima]
    return filtered_arcs, filtered_min, filtered_saddle, filtered_max


def classify_msc(image, classifier, persistence=8.5, show_pruned=False):
    predicted_classification = np.copy(image) / 1.1
    # np.zeros(image.shape)

    # Compute the morse smale complex of all blurred images and filter
    # arcs of interest based on mu from gaussian fit to pixel histogram
    # for each image

    blurred_image, fname = blur_and_save(image, "classify_msc")
    msc = build_msc(
        fname, blurred_image.shape[1], blurred_image.shape[0], persistence
    )
    arc_set_pixels = []
    # gauss_arcs = []
    pruned_cord = []
    for arc in msc.arcs:
        arc_pixels = []
        gauss_arc = []
        points = np.array(arc.line)
        # arc_points.append(tuple((point[1], point[0])))
        for point in points:
            arc_pixels.append(image[int(point[1]), int(point[0])])
        gauss_arc.append(gaussian_fit(np.array(arc_pixels)))
        arc_set_pixels.append(arc_pixels)
        # predict given obsereved arc attributes
        arc_classification = classifier.predict(gauss_arc)

        for point in points:
            if not show_pruned:
                predicted_classification[int(point[1]), int(point[0])] = (
                    arc_classification
                    if arc_classification == 1
                    else predicted_classification[int(point[1]), int(point[0])]
                )
            else:
                predicted_classification[int(point[1]), int(point[0])] = (
                    arc_classification
                    if arc_classification == 1
                    else predicted_classification[int(point[1]), int(point[0])]
                )
                if arc_classification != 1:
                    pruned_cord.append([int(point[1]), int(point[0])])

    show_image(predicted_classification)
    if show_pruned:
        pruned = np.vstack(pruned_cord)
        plt.scatter(
            pruned[:, 1],
            pruned[:, 0],
            facecolor=blue,
            edgecolor="none",
            s=1,
            marker=",",
            zorder=2,
            alpha=0.03,
        )


def filtered_msc_training_kernel(image, persistence=8.5):
    pos_pixels = []
    neg_pixels = []
    gauss_msc_pos = []
    gauss_msc_neg = []

    # Compute the morse smale complex of all blurred images and filter
    # arcs of interest based on mu from gaussian fit to pixel histogram
    # for each image
    blurred_image, fname = blur_and_save(image, "filtered_msc_training_kernel")
    full_msc = build_msc(
        fname, blurred_image.shape[1], blurred_image.shape[0], persistence
    )
    # mu, sigma, var = gaussian_fit(blurred_image)
    mu, _, _ = gaussian_fit(blurred_image)
    filtered_msc = msc_filtration(blurred_image, full_msc, mu)

    # Cover filtered MSC with closed balls to select negative space as
    # background and arcs as positive samples
    pos_samples, neg_samples = collect_training_from_filtered(
        blurred_image, filtered_msc, plot=True
    )

    # a = np.array(neg_samples)

    # get pixel values under arc for positives
    for arc_pos in pos_samples:
        for point_pos in arc_pos:
            pos_pixels.append([image[int(point_pos[0]), int(point_pos[1])]])

    # get negative pixel values
    for bg_group in neg_samples:
        for point_neg in bg_group:
            neg_pixels.append([image[int(point_neg[0]), int(point_neg[1])]])

    # fit gaussian to negative pixel values in background and take mean,
    # deviation, variation
    for arc_neg in neg_samples:
        filtered_negatives = []
        for point_neg in arc_neg:
            filtered_negatives.append(
                [image[int(point_neg[0]), int(point_neg[1])]]
            )
        # fit gaussian to negative sample pixels in background
        mu_neg, sigma_neg, var_neg = gaussian_fit(np.array(filtered_negatives))
        gauss_msc_neg.append([mu_neg, sigma_neg, var_neg])

    # fit gaussian to positive pixel values under arcs of MSC and take
    # mean, standard deviation, and variances
    for arc in pos_samples:
        filtered_positives = []
        for point in arc:
            filtered_positives.append([image[int(point[0]), int(point[1])]])
        # fit gaussian to positive pixel values under arc and take mean,
        # standard deviation, and variance
        # NOTE! can filter based on arc length, currently doing three bc
        # three degrees of freedom in gaussian fit.
        if len(filtered_positives) >= 3:
            mu_pos, sigma_pos, var_pos = gaussian_fit(
                np.array(filtered_positives)
            )
            gauss_msc_pos.append([mu_pos, sigma_pos, var_pos])

    train_x = np.vstack((pos_pixels, neg_pixels))
    train_y = np.hstack(
        (
            np.ones(np.array(pos_pixels).shape[0]),
            np.zeros(np.array(neg_pixels).shape[0]),
        )
    )

    gaussian_msc_arc_X = np.vstack((gauss_msc_pos, gauss_msc_neg))
    gaussian_msc_arc_Y = np.hstack(
        (
            np.ones(np.array(gauss_msc_pos).shape[0]),
            np.zeros(np.array(gauss_msc_neg).shape[0]),
        )
    )

    return train_x, train_y, gaussian_msc_arc_X, gaussian_msc_arc_Y


def collect_training_from_filtered(image, filtered_msc, width=10):
    # Collect Positive Samples
    fat_mask_image = make_dilated_arc_image(image, filtered_msc, width)

    # Collect Negative Samples
    background = np.array(zip(*np.where(fat_mask_image == 0)))
    ncluster = 16
    k_means = KMeans(init="k-means++", n_clusters=ncluster, n_init=10)
    k_means.fit(background)

    background_groups = []
    for k in range(ncluster):
        sub_group = []
        group_k_indx = zip(*np.where(k_means.labels_ == k))
        for loc_k in group_k_indx:
            sub_group.append(background[loc_k])
        background_groups.append(sub_group)

    return filtered_msc[0], background_groups


def mapped_msc_filtered_training_set(
    original_image, filtered_images, msc, filter_msc=False, bg_pixels=False
):  # , image_set = filtered_images):
    pos_pixels = []
    neg_pixels = []
    gauss_msc_pos = []
    gauss_msc_neg = []

    # compute morse smale complex of original, blurred image, to project
    # onto filtered images
    original_image = filtered_images[:, :, 0]

    # If the msc of original image is to be filtered
    if filter_msc:
        # mu, sigma, var = gaussian_fit(original_image)
        mu, _, _ = gaussian_fit(original_image)
        filtered_msc = msc_filtration(original_image, msc, mu)
        msc = filtered_msc
    else:
        # no filtration just sorted msc
        mu = 0
        msc = msc_filtration(original_image, msc, mu)

    filtered_images = np.dstack((filtered_images, original_image))
    kernel_size = filtered_images.shape[2]
    # print("kernel size: ", kernel_size)
    for i in range(kernel_size):
        # Cover filtered MSC with closed balls to select negative space
        # as background and arcs as positive samples
        # Using msc obtained from original image
        pos_samples, neg_samples = collect_training_from_filtered(
            filtered_images[:, :, i], msc
        )

        # get pixel values under arc for positives
        for arc_pos in pos_samples:
            for point_pos in arc_pos:
                pos_pixels.append(
                    [filtered_images[int(point_pos[0]), int(point_pos[1]), i]]
                )

        # get negative pixel values
        if bg_pixels:
            for bg_group in neg_samples:
                for point_neg in bg_group:
                    neg_pixels.append(
                        [
                            filtered_images[
                                int(point_neg[0]), int(point_neg[1]), i
                            ]
                        ]
                    )

            # fit gaussian to negative pixel values in background and
            # take mean, standard deviation, variance
            for arc_neg in neg_samples:
                filtered_negatives = []
                for point_neg in arc_neg:
                    filtered_negatives.append(
                        [
                            filtered_images[
                                int(point_neg[0]), int(point_neg[1]), i
                            ]
                        ]
                    )
                # fit gaussian to negative sample pixels in background
                mu_neg, sigma_neg, var_neg = gaussian_fit(
                    np.array(filtered_negatives)
                )
                # max_pix = np.max(np.array(filtered_negatives))
                # min_pix = np.min(np.array(filtered_negatives))
                gauss_msc_neg.append([mu_neg, sigma_neg, var_neg])

        # fit gaussian to positive pixel values under arcs of MSC and
        # take mean, standard deviation, and variance
        for arc in pos_samples:
            filtered_positives = []
            for point in arc:
                filtered_positives.append(
                    [filtered_images[int(point[0]), int(point[1]), i]]
                )
            # fit gaussian to positive pixel values under arc and take
            # mean, standard deviation, variance
            # NOTE! can filter based on arc length, currently doing
            # three bc three degrees of freedom in gaussian fit.
            if len(filtered_positives) >= 3:
                mu_pos, sigma_pos, var_pos = gaussian_fit(
                    np.array(filtered_positives)
                )
                # max_pix = np.max(np.array(filtered_positives))
                # min_pix = np.min(np.array(filtered_positives))
                gauss_msc_pos.append([mu_pos, sigma_pos, var_pos])

    train_x = np.vstack((pos_pixels, neg_pixels))
    train_y = np.hstack(
        (
            np.ones(np.array(pos_pixels).shape[0]),
            np.zeros(np.array(neg_pixels).shape[0]),
        )
    )

    gaussian_msc_arc_X = np.vstack((gauss_msc_pos, gauss_msc_neg))
    gaussian_msc_arc_Y = np.hstack(
        (
            np.ones(np.array(gauss_msc_pos).shape[0]),
            np.zeros(np.array(gauss_msc_neg).shape[0]),
        )
    )

    return train_x, train_y, gaussian_msc_arc_X, gaussian_msc_arc_Y


def get_points_from_arcs(arcs):
    points = tuple()
    for arc in arcs:
        points += tuple(np.array(np.round(arc.line), dtype=int))
    points = np.vstack(points)
    return points


def get_pixel_values_from_arcs(arcs, image):
    points = get_points_from_arcs(arcs)
    return image[points[:, 1], points[:, 0]].flatten()


def is_ridge_arc(arc, msc):
    return 0 not in [
        msc.nodes[arc.node_ids[0]].index,
        msc.nodes[arc.node_ids[1]].index,
    ]


def is_valley_arc(arc, msc):
    return 2 not in [
        msc.nodes[arc.node_ids[0]].index,
        msc.nodes[arc.node_ids[1]].index,
    ]

