import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA, FastICA, NMF
import math
from matplotlib.colors import rgb2hex

from bokeh.plotting import figure, show, output_file

from bokeh.models import HoverTool

from bokeh.models import ColumnDataSource

from sklearn.datasets import load_digits
from ipywidgets import interact

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

from sklearn.decomposition import FastICA, PCA
import seaborn as sns

sns.set(style='white', context='notebook')
import ipywidgets
from ipywidgets import interact, IntSlider

### Manifold learning dimensionality reduction: MDS and t-SNE


# Author: Jake Vanderplas -- <vanderplas@astro.washington.edu>


import scipy
from time import time
import warnings

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter

from sklearn import manifold, datasets

# Next line to silence pyflakes. This import is needed.
Axes3D

import fig_code


digits = load_digits()

def plot_image_components(x, coefficients=None, mean=0, components=None,
                          imshape=(8, 8), n_components=6, fontsize=12):
    if coefficients is None:
        coefficients = x

    if components is None:
        components = np.eye(len(coefficients), len(x))

    mean = np.zeros_like(x) + mean

    fig = plt.figure(figsize=(1.2 * (5 + n_components), 1.2 * 2))
    g = plt.GridSpec(2, 5 + n_components, hspace=0.3)

    def show(i, j, x, title=None):
        ax = fig.add_subplot(g[i, j], xticks=[], yticks=[])
        ax.imshow(x.reshape(imshape), interpolation='nearest')
        if title:
            ax.set_title(title, fontsize=fontsize)

    show(slice(2), slice(2), x, "True")

    approx = mean.copy()
    show(0, 2, np.zeros_like(x) + mean, r'$\mu$')
    show(1, 2, approx, r'$1 \cdot \mu$')

    for i in range(0, n_components):
        approx = approx + coefficients[i] * components[i]
        show(0, i + 3, components[i], r'$c_{0}$'.format(i + 1))
        show(1, i + 3, approx,
             r"${0:.2f} \cdot c_{1}$".format(coefficients[i], i + 1))
        plt.gca().text(0, 1.05, '$+$', ha='right', va='bottom',
                       transform=plt.gca().transAxes, fontsize=fontsize)

    show(slice(2), slice(-2, None), approx, "Approx")


def plot_pca_interactive(data, n_components=6):


    pca = PCA(n_components=n_components)
    Xproj = pca.fit_transform(data)

    def show_decomp(i=0):
        plot_image_components(data[i], Xproj[i],
                              pca.mean_, pca.components_)

    interact(show_decomp, i=(0, data.shape[0] - 1));


def explore_smushers():
    import macosko2015
    six_clusters, six_clusters_cells, six_clusters_genes = \
        macosko2015.load_big_clusters()

    hue = 'cluster_n_celltype'
    groupby = six_clusters_cells[hue]
    palette = fig_code.cluster_names_to_color

    csv = macosko2015.BASE_URL + 'differential_clusters_lowrank.csv'
    lowrank = pd.read_csv(csv, index_col=0)

    lowrank_big_clusters = lowrank.loc[six_clusters.index,
                                       six_clusters.columns]

    algorithms = {'PCA': PCA, 'ICA': FastICA, 'NMF': NMF}

    def _smush(n_components):

        smushers = {}
        smusheds = {}

        summaries = {}
        datas = {'big clusters (lowrank)': lowrank_big_clusters,
                 'big clusters': six_clusters}
        for data_name, data in datas.items():
            # print('data_name', data_name)
            for algo_name, algorithm in algorithms.items():
                # print('\talgo_name', algo_name)
                smusher = algorithm(n_components=n_components,
                                    random_state=2017)

                if data.min().min() < 0:
                    nonnegative = data - data.min().min()
                else:
                    nonnegative = data

                if 'digits' in data_name:
                    #             groupby = digits.target
                    index = digits.target
                elif 'amacrine' in data_name:
                    #             groupby = amacrine_cells['cluster_n']
                    index = data.index
                else:
                    #             groupby = cell_metadata['cluster_n_celltype']
                    index = data.index

                smushed = pd.DataFrame(smusher.fit_transform(nonnegative),
                                       index=index)
                smushed.columns.name = f'{algo_name} components'

                median = smushed.groupby(groupby).median()
                mean = smushed.groupby(groupby).mean()
                prefix = f'{data_name} {algo_name}'
                smusheds[prefix] = smushed
                smushers[prefix] = smusher

                summaries[f'{prefix}: mean'] = mean
                summaries[f'{prefix}: median'] = median
        return smusheds, smushers, summaries
    smusheds_n10, smushers_n10, summaries_n10 = _smush(n_components=10)

    sns.set(style='whitegrid', context='notebook')

    def plot_smushed(plot_type, algorithm, statistic, lowrank, n_components):

        if n_components != 10:
            smusheds_nX, smushers_nX, summaries_nX = _smush(n_components)
            smusheds, smushers, summaries = \
                smusheds_nX, smushers_nX, summaries_nX
        else:
            smusheds, smushers, summaries = \
                smusheds_n10, smushers_n10, summaries_n10

        key = 'big clusters'
        if lowrank:
            key += ' (lowrank)'
        key += f' {algorithm}'

        if plot_type == 'heatmap':
            # statistics = 'mean',  # 'median'
            key_summary = f'{key}: {statistic}'
            summary = summaries[key_summary]
            fig, ax = plt.subplots()
            sns.heatmap(summary)
            ax.set(title=key_summary)

        if plot_type == 'pairplot':
            smushed = smusheds[key]
            smushed_clusters = smushed.join(groupby)
            sns.pairplot(smushed_clusters, hue=hue, palette=palette)
            fig = plt.gcf()
            fig.suptitle(key)

    interact(plot_smushed, plot_type=['heatmap', 'pairplot'],
             algorithm=algorithms.keys(),
             statistic=['mean', 'median'], lowrank=False,
             n_components=IntSlider(value=10, min=2, max=10))



class MatrixGrid(sns.axisgrid.Grid):
    """Generalized plotter for many matrices"""

    def __init__(self, datas, size=3, sharex=True, sharey=False,
                 col_wrap=None, aspect=1, order=None):
        self.datas = datas

        n = len(datas) if order is None else len(order)

        if col_wrap is None:
            nrows = 1
            ncols = n
        else:
            ncols = col_wrap
            nrows = math.floor(n / col_wrap)

        print('nrows, ncols', nrows, ncols)

        figwidth = aspect * size * ncols
        figheight = size * nrows
        figsize = figwidth, figheight
        print('figsize', figsize)
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                 figsize=figsize, sharex=sharex, sharey=sharey)
        self.fig = fig
        self.axes = axes
        self.order = self.datas.keys() if order is None else order

    def map(self, function, *args, **kwargs):
        for ax, name in zip(self.axes.flat, self.order):
            #             print('data name', name)
            plt.sca(ax)
            data = self.datas[name]
            function(data, *args, **kwargs)
            ax.set_title(name)


def plot_bokeh_smushed(smushed_metadata, hue, prefix):

    categories = smushed_metadata[hue].unique()
    palette = map(rgb2hex, sns.color_palette('husl', n_colors=len(categories)))

    colormap = dict(zip(categories, palette))
    colors = [colormap[x] for x in groupby]
    smushed_metadata['color'] = colors

    hover = HoverTool(tooltips=[
        ("Cell Barcode", "@barcode"),
        ("Cluster #", '@cluster_n'),
        ('Celltype', '@celltype'),
        ('x, y', '($x, $y)')
    ])
    TOOLS = "crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"

    # p = figure(plot_height=600, plot_width=700, title="", toolbar_location=None, )
    p = figure(title=prefix, tools=[TOOLS, hover])

    xlabel = 0
    ylabel = 1
    # p.xaxis.axis_label = xlabel
    # p.yaxis.axis_label = ylabel

    source = ColumnDataSource(smushed_metadata)

    p.circle('0', '1',
             source=source,
             color='color',
             fill_alpha=0.2, size=10)
    show(p)





def load_pancreas():
    url = 'https://media.githubusercontent.com/media/olgabot/segerstolpe_palasantza2016/emily_V2/data/01_subset_metadata/rpkm_filter_cells_and_genes.txt'
    pancreas = pd.read_table(url, index_col=0).T.dropna(how='all', axis=1)
    url = 'https://media.githubusercontent.com/media/olgabot/segerstolpe_palasantza2016/emily_V2/data/01_subset_metadata/metadata_subset_792_cells.txt'
    pancreas_cells = pd.read_table(url, index_col=0)
    return pancreas, pancreas_cells


def compare_mds_tsne(algorithm='MDS', stddev=0, #metric='euclidean',
                     tsne_init='pca', random_state=0):
    cmap = plt.cm.viridis
    n_points = 1000
    # random_state = 0
    X, sample_order = datasets.samples_generator.make_s_curve(n_points,
                                                              random_state=0)

    if stddev > 0:
        X = X + np.random.normal(size=np.product(X.shape), scale=0.1).reshape(
            X.shape)

    fig = plt.figure(figsize=(8, 4))
    # plt.suptitle("Manifold Learning with %i points" % (n_points))

    # Plot original data
    ax = fig.add_subplot(121, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=sample_order, cmap=cmap)
    ax.view_init(4, -72)

    # Add noise if necessary
    n_components = 2
    mds_kws = dict(n_components=n_components, random_state=random_state)
    tsne_kws = dict(init=tsne_init)
    tsne_kws.update(mds_kws)

    # Perform MDS and plot it
    if algorithm == "MDS":
        smusher = manifold.MDS(max_iter=100, n_init=1, **mds_kws)
    # Perform t-SNE and plot it
    if algorithm == 'TSNE':
        smusher = manifold.TSNE(**tsne_kws)
    if algorithm == 'PCA':
        smusher = PCA(n_components=2)
    if algorithm == 'ICA':
        smusher = FastICA(n_components=2)
    if algorithm == 'NMF':
        smusher = NMF(n_components=2)
        X -= X.min()

    Y = smusher.fit_transform(X)

    ax = fig.add_subplot(1, 2, 2)
    plt.scatter(Y[:, 0], Y[:, 1], c=sample_order, cmap=cmap, linewidth=0.5,
                edgecolor='grey')
    plt.title(f"{algorithm}")
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')



def explore_manifold():
    ipywidgets.interact(
        compare_mds_tsne,
        algorithm=['TSNE', 'MDS', 'PCA', 'ICA', 'NMF'],
        # metric=ipywidgets.Dropdown(options=['euclidean', 'cityblock'],
        #                            value='euclidean', description="Distance Metric"),
        tsne_init=['random', 'pca'],
        n_points=ipywidgets.IntSlider(value=1000, min=50, max=2000, step=50),
        random_state=ipywidgets.IntSlider(value=0, min=0, max=10),
        stddev=ipywidgets.FloatSlider(value=0, min=0, max=1,
                                      description='Random Noise',
                                      step=0.1))