import warnings

import fastcluster
import matplotlib.pyplot as plt
import numpy as np
import polo
import seaborn as sns

FIGURE_FOLDER = 'figures'
DATA_FOLDER = 'data'


cluster_id_to_name = {24: 'Rods', 25: 'Cones', 26: 'Bipolar cells (group1)',
                      27: 'Bipolar cells (group2)',
                      33: 'Bipolar cells (group3)', 34: 'Muller glia'}
cluster_name_to_ids = {'Horizontal cells': 1, 'Retinal ganglion cells': 2,
                       'Amacrine cells': range(3, 24), "Rods": 24,
                       'Cones': 25, 'Bipolar cells': range(26, 34),
                       'Muller glia': 34, 'Astrocytes': 35,
                       'Fibroblasts': 36, 'Vascular endothelium': 37,
                       'Pericytes': 38, 'Microglia': 39}


cluster_ids = np.array([24, 25, 26, 27, 33, 34])

colors = sns.color_palette(palette='Set2', n_colors=len(cluster_ids))
id_to_color = dict(zip(cluster_ids, colors))
cluster_names_to_color = dict((cluster_id_to_name[i], id_to_color[i]) for i in cluster_ids)


def heatmap(data, ticklabels=False, **kwargs):
    """Visualize a 2d matrix, by default omitting zero-values and column/row labels
    
    Parameters
    ----------
    data : pandas.DataFrame or numpy.array
        Matrix of values to plot visually using colors
    ticklabels : bool
        If True, show both x- and y-ticklabels. If False, hide both
    
    Returns
    -------
    ax : matplotlib Axes
        Axes object with the heatmap
    """
    if 'mask' not in kwargs:
        mask = data == 0
    if not ticklabels:
        kwargs['xticklabels'] = []
        kwargs['yticklabels'] = []
    
    sns.heatmap(data, mask=mask, **kwargs)
    
def heatmaps(datas, axsize=3, **kwargs):
    """Visualize many heatmaps of matrices
    
    Parameters
    ----------
    datas : dict
        A mapping of the name of the matrix as a string to the matrix itself
    axsixe : int
        How big to make each individual axes
    """
    
    ncols = len(datas)
    nrows = 1

    width = ncols * axsize * 1.25
    height = nrows * axsize

    fig, axes = plt.subplots(ncols=ncols, figsize=(width, height))
    
    for ax, (key, data) in zip(axes, datas.items()):
        heatmap(data, ax=ax, **kwargs)
        ax.set(title=key)

# TODO
class MatrixGrid(sns.axisgrid.Grid):
    """Generalized plotter for many matrices"""
    
    def __init__(self, datas, size=3):
        pass


def optimal_linkage(data, rows=True, method='ward', metric='euclidean'):
    if not rows:
        data = data.T

    distance = fastcluster.pdist(data, metric=metric)
    linkage = fastcluster.linkage(distance, method=method)
    optimal_linkage = polo.optimal_leaf_ordering(linkage, distance)
    return optimal_linkage


def clustermap(data, xticklabels=[], yticklabels=[], method='ward',
               metric='euclidean', **kwargs):

    row_linkage = optimal_linkage(data, method=method, metric=metric)
    col_linkage = optimal_linkage(data, rows=False, method=method,
                                  metric=metric)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sns.clustermap(data, row_linkage=row_linkage,
                              xticklabels=xticklabels, yticklabels=yticklabels,
                              col_linkage=col_linkage, **kwargs)