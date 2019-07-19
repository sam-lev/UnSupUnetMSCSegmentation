# Standard library imports

# Third party imports
import numpy as np

# Local application imports


def gaussian_fit(image, plot=False):
    n, bins_ = np.histogram(image.flatten())
    mids = 0.5 * (bins_[1:] + bins_[:-1])
    mu = np.average(mids, weights=n)
    var = np.average((mids - mu) ** 2, weights=n)
    sigma = np.sqrt(var)
    # right_inflection = mu + sigma
    return mu, sigma, var  # , right_inflection
