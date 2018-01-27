################################################
# ---------- NBS Plotting Functions ---------- #
################################################
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test as multiv_lr_test

# Helper function for turning cluster assignments into color mappings (used for consensus clustering map figures)
def cluster_color_assign(cluster_assignments, name=None):
    k = max(cluster_assignments.value_counts().index)
    colors = sns.color_palette('hls', k)
    cluster_cmap = {i:colors[i-1] for i in range(1, k+1)}
    pat_colors = {}
    for pat in cluster_assignments.index:
        pat_colors[pat] = cluster_cmap[cluster_assignments.ix[pat]]
    cluster_cmap = pd.Series(pat_colors, name=name)
    return cluster_cmap

# Function for plotting consensus clustering map
# Needs both the consensus clustering similarity table and linkage map from assignment
# Actual cluster assignments on col_color_map
# Cluster assignments to be compared passed to row_color_map
# If there are multiple mappings for row_color_map, it can be passed as a dataframe with the index space of the cc_table
def plot_cc_map(cc_table, linkage, row_color_map=None, col_color_map=None, verbose=True, **save_args):
    title = 'Co-Clustering Map'
    if 'job_name' in save_args:
        title = save_args['job_name']+' Co-Clustering Map'
    plt.figure(figsize=(20,20))
    cg = sns.clustermap(cc_table, row_linkage=linkage, col_linkage=linkage, 
                        cmap='Blues', cbar_kws={'label': 'Co-Cluster Frequency'},
                        row_colors=row_color_map, col_colors=col_color_map, 
                        **{'xticklabels':'False', 'yticklabels':'False'})
    cg.cax.set_position([0.92, .11, .03, .584])
    cg.ax_heatmap.set_xlabel('')
    cg.ax_heatmap.set_xticks([])
    cg.ax_heatmap.set_ylabel('')
    cg.ax_heatmap.set_yticks([])
    cg.ax_row_dendrogram.set_visible(False)
    plt.suptitle(title, fontsize=20, x=0.6, y=0.95)
    if 'outdir' in save_args:
        if 'job_name' in save_args:
            save_cc_map_path = save_args['outdir']+str(save_args['job_name'])+'_cc_map.png'
        else:
            save_cc_map_path = save_args['outdir']+'cc_map.png'
        plt.savefig(save_cc_map_path, bbox_inches='tight')
        plt.show()
    if verbose:
        print 'Co-Clustering Map plotted'
    return

# Function for plotting Kaplan Meier plot of cluster survivals
# Requires lifelines package
# clin_data_fn is the the clinical data of TCGA cohort from Broad Firehose
# cluster_assign ias a pandas Series of the patient cluster assignments from NBS with patient ID's as the index
# tmax is the maximum plot duration for the KMplot, but the logrank test always calculates to longest survival point
def cluster_KMplot(cluster_assign, clin_data_fn, delimiter='\t', lr_test=True, tmax=-1, verbose=True, **save_args):
    title = 'KM Survival Plot'
    if 'job_name' in save_args:
        title = save_args['job_name']+' KM Survival Plot'

    # Initialize KM plotter
    kmf = KaplanMeierFitter()
    # Load and format clinical data   
    surv = pd.read_csv(clin_data_fn, sep=delimiter, index_col=0)
    # Number of clusters
    clusters = sorted(list(cluster_assign.value_counts().index))
    k = len(clusters)
    # Initialize KM Plot Settings
    fig = plt.figure(figsize=(10, 7)) 
    ax = plt.subplot(1,1,1)
    colors = sns.color_palette('hls', k)
    cluster_cmap = {clusters[i]:colors[i] for i in range(k)}
    # Plot each cluster onto KM Plot
    for clust in clusters:
        clust_pats = list(cluster_assign[cluster_assign==clust].index)
        clust_surv_data = surv.ix[clust_pats].dropna()
        kmf.fit(clust_surv_data.overall_survival, clust_surv_data.vital_status, label='Group '+str(clust)+' (n=' +  str(len(clust_surv_data)) + ')')
        kmf.plot(ax=ax, color=cluster_cmap[clust], ci_show=False)
    # Set KM plot limits to 5 years and labels
    # if tmax!=-1:
    plt.xlim((0,1825))
    plt.xlabel('Time (Days)', fontsize=16)
    plt.ylabel('Survival Probability', fontsize=16)
    # Multivariate logrank test
    if lr_test:
        cluster_survivals = pd.concat([surv, cluster_assign], axis=1).dropna().astype(int)
        p = multiv_lr_test(np.array(cluster_survivals.overall_survival), 
                           np.array(cluster_survivals[cluster_assign.name]), t_0=tmax,
                           event_observed=np.array(cluster_survivals.vital_status)).p_value
        if verbose:
            print 'Multi-Class Log-Rank P:', p
        plt.title(title+'\np='+repr(round(p, 4)), fontsize=24, y=1.02)
    else:
        plt.title(title, fontsize=24, y=1.02)
    # Save KM plot
    if 'outdir' in save_args:
        if 'job_name' in save_args:
            save_KMplot_path = save_args['outdir']+str(save_args['job_name'])+'_KM_plot.png'
        else:
            save_KMplot_path = save_args['outdir']+'KM_plot.png'
        plt.savefig(save_KMplot_path, bbox_inches='tight')
        plt.show()
    if verbose:
        print 'Kaplan Meier Plot constructed'
    if lr_test:
        return p
    else:
        return