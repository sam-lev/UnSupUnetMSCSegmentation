# Standard library imports

# Third party imports
import numpy as np
from skimage import morphology, exposure
from skimage.io import imsave

# Local application imports
from topoml.image.feature import gaussian_blur_filter


def make_arc_image(image, msc,labeled_accuracy=1, invert=False):
    arc_mask_image = np.zeros(image.shape)
    print("SHAPEEEEE")
    print(image.shape)
    x = 0 if not invert else 1
    y = 1 if not invert else 0
    for a in msc.arcs:
        points = np.array(a.line)
        for point in points:
            arc_mask_image[int(point[x]), int(point[y])] = a.label_accuracy if a.label_accuracy is not None else 1
    return arc_mask_image


def make_mc_arc_image(image, msc, labeled_accuracy=1, invert=False):
    arc_mask_image = np.zeros(image.shape)

    mask_index = 2 if invert else 0

    x = 0 if invert else 1
    y = 1 if invert else 0

    for a in msc.arcs:
        if mask_index not in [
            msc.nodes[a.node_ids[0]].index,
            msc.nodes[a.node_ids[1]].index,
        ]:
            points = np.array(a.line)
            for point in points:
                arc_mask_image[int(point[x]), int(point[y])] = a.label_accuracy if a.label_accuracy is not None else 1
    return arc_mask_image


def make_dilated_arc_image(image, msc, width,labeled_accuracy=1, invert=True):
    return morphology.dilation(
        make_arc_image(image, msc,labeled_accuracy=labeled_accuracy, invert=invert), selem=morphology.disk(width)
    )


def make_arc_mask(image, msc, labeled_accuracy=1, invert=False):
    arc_mask_image = make_arc_image(image, msc, labeled_accuracy=labeled_accuracy,invert=invert)
    return np.ma.masked_where(arc_mask_image == 0, arc_mask_image)


def make_mc_arc_mask(image, msc,labeled_accuracy=1, invert=False):
    arc_mask_image = make_mc_arc_image(image, msc,labeled_accuracy=labeled_accuracy, invert=invert)
    return np.ma.masked_where(arc_mask_image == 0, arc_mask_image)


def blur_and_save(original_image, fname_base, blur_sigma=2, grey_scale=True):
    blurred_image = gaussian_blur_filter(original_image, sigma=blur_sigma, as_grey=grey_scale).astype(
        "float32"
    )
    fname_raw = fname_base + "_smoothed.raw"
    blurred_image.tofile(fname_raw)
    return blurred_image, fname_raw

def augment_channels(original_image, fname_base, channels = [0,1]):
    import copy
    import cv2
    augmented_image = copy.deepcopy(original_image)
    
    #[ 0 = blue, 1 = green, 2 = red ]
    for c in channels:
        augmented_image[:,:,c] = 0#cv2.equalizeHist(augmented_image[:,:,c])
    fname_components = fname_base.split('.')
    fname_aug = fname_components[0]+"_aug."+fname_components[1]
    imsave(fname_aug, augmented_image)
    return augmented_image, fname_aug

def scale_intensity(original_image, fname_base, scale_range=(0,255)):
    scaled_image = exposure.rescale_intensity(original_image)#, in_range=(0, 255))
    fname_raw = fname_base + "_scaled.raw"
    scaled_image.tofile(fname_raw)
    return scaled_image, fname_raw  


def bounding_box(points):
    points = np.array(points)
    xmin = np.min(points[:, 0])
    xmax = np.max(points[:, 0])
    ymin = np.min(points[:, 1])
    ymax = np.max(points[:, 1])

    return (xmin, xmax, ymin, ymax)


def range_overlap(a_min, a_max, b_min, b_max):
    return (a_min <= b_max) and (b_min <= a_max)


def box_intersection(b1, b2):
    return range_overlap(b1[0], b1[1], b2[0], b2[1]) and range_overlap(
        b1[2], b1[3], b2[2], b2[3]
    )
