# Interactive jupyter widgets - use IntSlider directly for more control
from ipywidgets import IntSlider, interact

# Convert RGB colors to hex colors for portability
from matplotlib.colors import rgb2hex
import matplotlib.pyplot as plt

# Visualize networks
import networkx

# Numerical python
import numpy as np

# Pandas for dataframes
import pandas as pd

# K-nearest neighbors cell clustering from Dana Pe'er's lab
import phenograph

# Make color palettes
import seaborn as sns

# Bokeh - interactive plotting in the browser
from bokeh.plotting import figure, show, output_file
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.models.widgets import Panel, Tabs

import macosko2015

# Bokeh tools in sidebar
TOOLS = "crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset," \
        "tap,save,box_select,poly_select,lasso_select"
X_COL = 'xs'
Y_COL = 'ys'



def correlation_to_distance(correlations):
    """Convert -1 to 1 similarity measure to a metric distance"""
    return np.sqrt(2*(1-correlations))


def plot_graph(nodes_source, edges_source, legend_col, color_col, title,
               tab=False):
    """Draw a network in 2d space using Bokeh"""

    # names=['cell'] matches p1.circle(... name='cell' ...) so that the
    # hover tooltip only shows up on the cells, not on the edges
    hover1 = HoverTool(names=['cell'], tooltips=[
        ("Cell Barcode", "@barcode"),
        ("Group", f"@{legend_col}"),
    ])

    p1 = figure(tools=[TOOLS, hover1], plot_height=500, plot_width=750,
                title=title)
    edges = p1.multi_line(X_COL, Y_COL, line_width=1.5,
                          alpha='alphas', color='black',
                          source=edges_source)
    nodes = p1.circle(X_COL, Y_COL, size=10,
                      fill_alpha=0.5, hover_alpha=1,
                      hover_fill_color=color_col,
                      muted_alpha=0.2,
                      source=nodes_source, legend=legend_col,
                      color=color_col, name='cell')
    p1.legend.click_policy = "mute"
    if tab:
        tab1 = Panel(child=p1, title=title)
        return tab1
    else:
        return p1


def get_nodes_specs(positions, metadata, barcodes, communities,
                    community_col='community',
                    other_cluster_col=None, palette=None):
    """Create a dataframe of x,y position of nodes
    
    Parameters
    ----------
    positions : dict
        Mapping of each node's id to its x, y position
    metadata : pandas.DataFrame
        Any additional data describing each node
    barcodes : list
        Identifiers of each node, must correspond to the index of metadata, and
        must be in the exact order found in `positions`
    communities : list
        Cluster assignments of each cell
    community_col : str, optional
        What to name the column in the dataframe for where the communities are 
        stored
    other_cluster_col : str, optional
        If not None, then this is a name of a colujmn that exists in `metadata`
        that you'd like to also create colors for
    palette : str, optional
        Name of a color palette to use. Default is "deep"
    """
    layout = pd.DataFrame(positions, index=['xs', 'ys']).T
    layout[community_col] = communities
    layout[community_col] = f'{community_col} #'.capitalize() \
                            + layout[community_col].astype(str)
    layout = layout.sort_values(community_col)

    layout['barcode'] = barcodes
    layout_metadata = layout.join(metadata, on='barcode')
    if other_cluster_col is not None:
        layout_metadata = layout_metadata.sort_values(other_cluster_col)
        layout_metadata['other_cluster_color'] = labels_to_colors(
            layout_metadata[other_cluster_col], palette)
    layout_metadata[f'{community_col}_color'] = labels_to_colors(
        layout_metadata[community_col], palette)
    return layout_metadata


def labels_to_colors(labels, palette=None):
    """Convert a column of labels into categorical colors"""
    categories = sorted(np.unique(labels))
    rgb = sns.color_palette(palette, n_colors=len(categories))
    colors = map(rgb2hex, rgb)
    colormap = dict(zip(categories, colors))
    return [colormap[x] for x in labels]


def get_edges_specs(graph, positions):
    """Create a dataframe of x,y positions of edges
    
    Parameters
    ----------
    graph : networkx.Graph
        A networkx object containing nodes connected to each other via edges
    positions : dict
        Mapping of each node's id to its x, y position
    
    Returns
    -------
    df : pandas.DataFrame
        A table containing the columns, 'xs', 'ys', and 'alphas' describing the
        location and transparency (alpha) of each edge
    """
    df = pd.DataFrame(columns=['xs', 'ys', 'alphas'],
                      index=range(len(graph.edges())))
    weights = [d['weight'] for u, v, d in graph.edges(data=True)]
    max_weight = max(weights)
    calc_alpha = lambda h: 0.1 + 0.6 * (h / max_weight)

    # example: { ..., ('user47', 'da_bjoerni', {'weight': 3}), ... }
    for i, (u, v, data) in enumerate(graph.edges(data=True)):
        # Assign 'xs' column
        df.iloc[i, 0] = [positions[u][0], positions[v][0]]
        # Assign 'ys' column
        df.iloc[i, 1] = [positions[u][1], positions[v][1]]
        # Assign 'alphas' column
        df.iloc[i, 2] = calc_alpha(data['weight'])
    return df


def explore_phenograph():
    """Interactively shows KNN graphs and community detection"""
    big_clusters, big_clusters_cells, big_clusters_genes = \
        macosko2015.load_big_clusters()
    amacrine, amacrine_cells, amacrine_genes = macosko2015.load_amacrine()

    # --- Read the "lowrank" or "smoothed" data from robust PCA --- #
    csv = macosko2015.BASE_URL + 'differential_clusters_lowrank.csv'
    lowrank = pd.read_csv(csv, index_col=0)

    lowrank_big_clusters = lowrank.loc[
        big_clusters.index, big_clusters.columns]

    lowrank_amacrine = lowrank.loc[amacrine.index, amacrine.columns]

    datasets = {'big clusters lowrank': lowrank_big_clusters,
                'big clusters': big_clusters,
                'amacrine lowrank': lowrank_amacrine,
                'amacrine': amacrine}

    correls = {k: correlation_to_distance(data.T.rank().corr())
               for k, data in datasets.items()}

    def plot_phenograph(dataset='big clusters',
                        primary_metric='euclidean',
                        lowrank=False,
                        k=30, min_cluster_size=10):
        key = dataset + ' lowrank' if lowrank else dataset
        corr = correls[key]

        if dataset == 'big clusters':
            metadata = big_clusters_cells
            palette = 'Set2'
            cluster_col = 'cluster_id'
        elif dataset == 'amacrine':
            metadata = amacrine_cells
            palette = 'husl'
            cluster_col = 'cluster_id'
        community_col = 'community'

        communities, graph, Q = phenograph.cluster(
            corr, k=k, primary_metric=primary_metric,
            min_cluster_size=min_cluster_size)
        network = networkx.from_scipy_sparse_matrix(graph)
        positions = networkx.spring_layout(network)

        nodes_source = ColumnDataSource(get_nodes_specs(
            positions, metadata, corr.index, communities,
            other_cluster_col=cluster_col,
            community_col=community_col, palette=palette))
        edges_source = ColumnDataSource(get_edges_specs(network, positions))

        # --- First tab: KNN clustering --- #
        tab1 = plot_graph(nodes_source, edges_source, legend_col=community_col,
                          color_col=f'{community_col}_color', tab=True,
                          title='KNN Clustering')

        # --- Second tab: Clusters from paper --- #
        tab2 = plot_graph(nodes_source, edges_source,
                          legend_col='cluster_n_celltype', tab=True,
                          color_col='other_cluster_color',
                          title="Clusters from paper")

        tabs = Tabs(tabs=[tab1, tab2])
        show(tabs)

    interact(plot_phenograph, dataset=['big clusters', 'amacrine'],
             primary_metric=['euclidean', 'manhattan'],
             k=IntSlider(start=3, stop=100, value=30, step=5),
             min_cluster_size=IntSlider(start=3, stop=100, value=10, step=5))


def plot_kmeans_interactive(min_clusters=1, max_clusters=6):
    """From https://github.com/jakevdp/sklearn_tutorial/"""
    from IPython.html.widgets import interact
    from sklearn.metrics.pairwise import euclidean_distances
    from sklearn.datasets.samples_generator import make_blobs

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')

        X, y = make_blobs(n_samples=300, centers=4,
                          random_state=0, cluster_std=0.60)

        def _kmeans_step(frame=0, n_clusters=4):
            rng = np.random.RandomState(2)
            labels = np.zeros(X.shape[0])
            centers = rng.randn(n_clusters, 2)

            nsteps = frame // 3

            for i in range(nsteps + 1):
                old_centers = centers
                if i < nsteps or frame % 3 > 0:
                    dist = euclidean_distances(X, centers)
                    labels = dist.argmin(1)

                if i < nsteps or frame % 3 > 1:
                    centers = np.array([X[labels == j].mean(0)
                                        for j in range(n_clusters)])
                    nans = np.isnan(centers)
                    centers[nans] = old_centers[nans]

            # plot the data and cluster centers
            plt.scatter(X[:, 0], X[:, 1], c=labels, s=50, cmap='rainbow',
                        vmin=0, vmax=n_clusters - 1);
            plt.scatter(old_centers[:, 0], old_centers[:, 1], marker='o',
                        c=np.arange(n_clusters),
                        s=200, cmap='rainbow')
            plt.scatter(old_centers[:, 0], old_centers[:, 1], marker='o',
                        c='black', s=50)

            # plot new centers if third frame
            if frame % 3 == 2:
                for i in range(n_clusters):
                    plt.annotate('', centers[i], old_centers[i],
                                 arrowprops=dict(arrowstyle='->', linewidth=1))
                plt.scatter(centers[:, 0], centers[:, 1], marker='o',
                            c=np.arange(n_clusters),
                            s=200, cmap='rainbow')
                plt.scatter(centers[:, 0], centers[:, 1], marker='o',
                            c='black', s=50)

            plt.xlim(-4, 4)
            plt.ylim(-2, 10)

            if frame % 3 == 1:
                plt.text(3.8, 9.5, "1. Reassign points to nearest centroid",
                         ha='right', va='top', size=14)
            elif frame % 3 == 2:
                plt.text(3.8, 9.5, "2. Update centroids to cluster means",
                         ha='right', va='top', size=14)

    return interact(_kmeans_step, frame=[0, 50],
                    n_clusters=[min_clusters, max_clusters])