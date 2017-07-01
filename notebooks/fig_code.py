import warnings

import fastcluster
from ipywidgets import interact, IntSlider
import macosko2015
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polo
from scipy.spatial import distance
import seaborn as sns



FIGURE_FOLDER = 'figures'
DATA_FOLDER = 'data'

expression, cell_metadata, gene_metadata = macosko2015.load_big_clusters()
cluster_ids_unique = cell_metadata['cluster_id'].unique()

cluster_n_to_name = {24: 'Rods', 25: 'Cones',
                      26: 'Bipolar cells\n(group1)',
                      27: 'Bipolar cells\n(group2)',
                      33: 'Bipolar cells\n(group3)',
                      34: 'Muller glia'}

cluster_id_to_name = dict(('cluster_{}'.format(str(i).zfill(2)), name)
                          for i, name in cluster_n_to_name.items())

colors = sns.color_palette(palette='Set2', n_colors=len(cluster_ids_unique))
id_to_color = dict(zip(cluster_ids_unique, colors))

color_labels = [id_to_color[i] for i in cell_metadata['cluster_id']]
cluster_names_to_color = dict((cluster_id_to_name[i], id_to_color[i])
                              for i in cluster_ids_unique)
cluster_names_to_color = pd.Series(cluster_names_to_color)




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
    """Wraps seaborn clustermap with optimal leaf ordering and no ticklabels"""

    row_linkage = optimal_linkage(data, method=method, metric=metric)
    col_linkage = optimal_linkage(data, rows=False, method=method,
                                  metric=metric)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sns.clustermap(data, row_linkage=row_linkage,
                              xticklabels=xticklabels, yticklabels=yticklabels,
                              col_linkage=col_linkage, **kwargs)


def plot_dropout(percent_gene_dropout=50,
                 correlation='pearson', distance_metric='euclidean',
                 linkage_method='ward'):

    title = '{}% dropout, {} correl, {} distance, {} linkage'.format(
        percent_gene_dropout, correlation, distance_metric, linkage_method)

    threshold = percent_gene_dropout / 100.
    print('threshold', threshold)

    # Use a boolean (0 and 1) mask to randomly choose genes to stay (equal to
    # 1) or leave (equal to 0)
    mask = np.random.uniform(size=expression.shape) > threshold
    # print('mask.shape', mask.shape)

    # '~' means 'not' so it will flip 0s to 1s and make it easier to sum
    removed = pd.DataFrame(~mask)
    print("Mean genes removed per cell:", removed.sum(axis=1).mean())
    print("Median genes removed per cell:", removed.sum(axis=1).median())
    data = expression.multiply(mask)
    # print('data.shape', data.shape)
    # print(data.head())

    correlated = data.T.corr(method=correlation)

    g = clustermap(correlated,
                   col_colors=color_labels,
                   row_colors=color_labels,
                   metric=distance_metric,
                   method=linkage_method,
                   figsize=(4, 4))
    g.fig.suptitle(title)
    # plt.show()


def plot_dropout_interactive():
    interact(plot_dropout,
             percent_gene_dropout=IntSlider(value=0, min=0, max=90, step=10),
             correlation=['pearson', 'spearman'],
             distance_metric=['euclidean', "cityblock"],
             linkage_method=['ward', 'average', 'single', "complete"])


def plot_color_legend():
    """Show how colors map back to the labeled cell types"""
    sns.palplot(cluster_names_to_color)
    ax = plt.gca()
    xticks = np.arange(0, cluster_names_to_color.shape[0])
    ax.set(xticklabels=cluster_names_to_color.index, xticks=xticks)
    ax.grid(False)
    # plt.show()

### --- Anscombe's quartet --- ###
# Convienence function that reads the data
# Equivalent to:
# anscombe = pd.read_csv("https://github.com/mwaskom/seaborn-data/raw/master/anscombe.csv")
# But is much easier to write
anscombe = sns.load_dataset('anscombe')


def explore_anscombe(summary):

    statistical = ('mean', 'var', 'std')
    grouped = anscombe.groupby('dataset')

    col = None

    if summary in statistical:
        summarized = getattr(grouped, summary)()
        tidy = summarized.unstack().reset_index()
        tidy = tidy.rename(columns={'level_0': 'variable', 0: summary})
        col = 'variable'
    else:
        if summary.endswith('correlation'):
            method = summary.split()[0].lower()
            summarized = grouped.apply(
                lambda df: df['x'].corr(df['y'], method=method))
        elif summary.endswith('distance'):
            metric = getattr(distance, summary.split()[0].lower())
            summarized = grouped.apply(lambda df: metric(df['x'], df['y']))
        tidy = summarized.reset_index()
        tidy = tidy.rename(columns={'index': 'variable', 0: summary})
    print(summarized.T)

    g = sns.factorplot(data=tidy, col=col, x='dataset',
                       y=summary, kind='bar', size=3, zorder=-1)
    for ax in g.axes.flat:
        # add a white grid on top
        ax.grid(axis='y', color='white', zorder=2)
    g.fig.suptitle(f'summary: {summary}')
    # plt.show()


def interact_anscombe():
    interact(explore_anscombe,
             summary=['mean', 'var', 'std','Pearson correlation',
                      'Spearman correlation', 'Euclidean distance',
                      'Cityblock distance'])