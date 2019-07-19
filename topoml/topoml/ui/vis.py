# Standard library imports
import collections

# Third party imports
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.stats import norm
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

# Local application imports
from topoml.ui.colors import orange, purple, blue, bold_palette


def show_image(image):
    plt.figure(num=None, dpi=200)
    plt.imshow(image, cmap=plt.cm.Greys_r)
    plt.gca().axis("off")


def show_training_samples(image, in_samples, out_samples):
    show_image(image)
    plt.scatter(
        in_samples[:, 0], in_samples[:, 1], color=orange, s=1, marker=","
    )
    plt.scatter(
        out_samples[:, 0], out_samples[:, 1], color=purple, s=1, marker=","
    )
    plt.show()


def show_histogram(image):
    (mu, sigma) = norm.fit(image.flatten())
    n, bins, patches = plt.hist(
        image.flatten(), 60, normed=1, facecolor="blue", alpha=0.75
    )

    # add a 'best fit' line
    y = mlab.normpdf(bins, mu, sigma)
    # l = plt.plot(bins, y, "r--", linewidth=2)

    # inflection points at x = mean +/- standard deviation
    # inflection_2 = mu + sigma
    x = mu + sigma
    plt.scatter(x, y=norm.pdf(x), s=3, marker=",", alpha=0)
    # plot
    plt.xlabel("pixel value")
    plt.ylabel("probability")
    plt.title(r" mu=%.3f, \sigma=%.3f$" % (mu, sigma))
    plt.grid(True)

    # Peak value probability density function at mean
    # gauss_peak = norm.pdf(mu)
    print("inflection point = ", mu + sigma, ",", norm.pdf(mu + sigma))
    print("peak = ", mu, ",", norm.pdf(mu))

    """ sanity check """
    """
    n_, bins_ = np.histogram(image.flatten())
    mids = 0.5*(bins_[1:] + bins_[:-1])
    mu_ = np.average(mids, weights=n_)
    var = np.average((mids - mu_)**2, weights=n_)
    sigma_ = np.sqrt(var)
    print('mu', mu_, 'var',var,'sigma',sigma_)
    """
    plt.show(block=True)


def show_filtered_msc(image, filtered_msc):
    filtered_arcs = filtered_msc[0]
    filtered_min = filtered_msc[1]
    filtered_saddle = filtered_msc[2]
    filtered_max = filtered_msc[3]

    show_image(image)

    plt.scatter(
        np.array(filtered_arcs)[:, 1],
        np.array(filtered_arcs)[:, 0],
        facecolor="orange",
        edgecolor="none",
        s=1,
        marker=",",
        alpha=0.5,
    )
    plt.scatter(
        np.array(filtered_min)[:, 1],
        np.array(filtered_min)[:, 0],
        facecolor=bold_palette[0],
        edgecolor="none",
        marker=",",
        s=2,
    )
    plt.scatter(
        np.array(filtered_saddle)[:, 1],
        np.array(filtered_saddle)[:, 0],
        facecolor=bold_palette[1],
        edgecolor="none",
        marker=",",
        s=2,
    )
    plt.scatter(
        np.array(filtered_max)[:, 1],
        np.array(filtered_max)[:, 0],
        facecolor=bold_palette[2],
        edgecolor="none",
        marker=",",
        s=2,
    )

    plt.gca().set_xlim(0, image.shape[1])
    plt.gca().set_ylim(image.shape[0], 0)
    plt.show(block=True)


def show_pixel_classifier_result(compiled_data, classifier):
    flattened_data = compiled_data.reshape(-1, compiled_data.shape[-1])
    predicted_labels = classifier.predict(flattened_data)
    predicted_labels = predicted_labels.reshape(
        compiled_data.shape[0], compiled_data.shape[1]
    )
    show_image(predicted_labels)
    plt.show(block=True)


def show_spectral_graph_embedding(graph):
    plt.figure()
    cmap = [purple, orange, blue]
    nodes = graph.nodes()
    colors = [cmap[graph.node[n]['label']] for n in nodes]
    nx.draw_spectral(graph, node_size=10, nodelist=nodes, node_color=colors)
    plt.show()


def show_spring_graph_layout(graph):
    plt.figure()
    Gcc = sorted(
        nx.connected_component_subgraphs(graph), key=len, reverse=True
    )[0]

    cmap = [purple, orange, blue]
    nodes = graph.nodes()
    colors = [cmap[graph.node[n]['label']] for n in nodes]
    pos = nx.spring_layout(Gcc)
    plt.axis("off")
    nx.draw_networkx_edges(Gcc, pos, alpha=0.4)
    nx.draw_networkx_nodes(Gcc, pos, node_size=10, nodelist=nodes, node_color=colors)
    plt.show()


def show_graph_degree_histogram(graph):
    plt.figure()
    degree_sequence = sorted([d for n, d in graph.degree()], reverse=True)
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    plt.bar(deg, cnt, width=0.80, color="b")
    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    plt.gca().set_xticks([d + 0.4 for d in deg])
    plt.gca().set_xticklabels(deg)
    plt.show()


def show_pca(X, labels=None):
    X = np.array(X)
    if labels is None:
        labels = -1*np.ones(X.shape[0], dtype=int)
    colors = np.array([purple, orange, blue])
    plt.figure()
    pca = PCA(n_components=X.shape[1])
    X_reduced = pca.fit_transform(X)
    plt.scatter(X_reduced[:, 0], X_reduced[:, 1], c=colors[labels], s=10)
    plt.show()
    return pca

def show_lda(X, labels):
    X = np.array(X)
    plt.figure()
    lda = LinearDiscriminantAnalysis(solver='svd', shrinkage=None, priors=None, n_components=None, store_covariance=False, tol=0.0001)
    X_reduced = lda.fit_transform(X, labels)
    colors = np.array([purple, orange, blue])
    plt.scatter(X_reduced[:, 0], labels, c=colors[labels], s=10)
    plt.show()
    return lda
