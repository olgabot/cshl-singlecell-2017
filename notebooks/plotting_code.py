
import seaborn as sns

def violinplot_grid(tidy, col='gene_symbol', size=4, aspect=0.25, gridspec_kws=dict(wspace=0),
                    sharex=False, scale='width', linewidth=1, palette='husl', inner=None, 
                    cut=True, order=None):
    facetgrid = sns.FacetGrid(tidy, col=col, size=size, aspect=aspect,
                  gridspec_kws=gridspec_kws, sharex=sharex)
    facetgrid.map(sns.violinplot, 'expression_log', 'cluster_id', scale=scale, 
          linewidth=linewidth, palette=palette, inner=inner, cut=cut, order=order)
    facetgrid.set(xlabel='', xticks=[])
    facetgrid.set_titles('{col_name}')