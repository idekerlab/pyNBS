############################################################
# ---------- NBS Consensus Clustering Functions ---------- #
############################################################

import os
import pandas as pd
import numpy as np
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hclust

# Takes a list of 'H' (patient-by-k) dataframes and performs 'hard' consensus clustering
# Using hierarchical clustering and average linkage
# Returns similarity table (distance is 1-similarity) and linkage map of patients
# Also returns cluster assignment map of patients if wanted
def consensus_hclust_hard(Hlist, k=3, hclust_linkage_method='average',
    hclust_linkage_metric='euclidean', verbose=True, **save_args):
    # Make sure all H matrices are pandas DataFrames
    if not all([type(H)==pd.DataFrame for H in Hlist]):
        raise ValueError('Not all H matrices given are pandas DataFrames.')
    # Default make sure that the number of clusters is <= the number of columns in H
    if not all([Hlist[0].shape[1]==k for H in Hlist]):
        raise ValueError('All H matrices must have the same number of columns as number of clusters (k)')
    # Generate patient list
    pat_list = set()
    for H in Hlist:
        pat_list = pat_list.union(set(H.index))
    pat_list = sorted(list(pat_list))
    if verbose:
        print 'Constructing Hlist:', len(Hlist), 'cluster matrices, ', len(pat_list), 'samples'

    # Initialzie co-clustering tables
    co_clust_table = pd.DataFrame(0, index=pat_list, columns=pat_list)
    cluster_count = pd.DataFrame(0, index=pat_list, columns=pat_list)

    # Calculate patient similarities and linkage
    for H in Hlist:
        H.columns = range(1,len(H.columns)+1)
        # Update patient cluster count
        cluster_count.ix[H.index, H.index]+=1
        # Get cluster assignment for each patient
        cluster_assign = {i:[] for i in H.columns}
        for pat in H.index:
            cluster_assign[np.argmax(H.ix[pat])].append(pat)
        # Update co-clustering matrix with each cluster assignment
        for cluster in cluster_assign:
            cluster_pats = cluster_assign[cluster]
            co_clust_table.ix[cluster_pats, cluster_pats]+=1
    cc_hard_sim_table = co_clust_table.astype(float).divide(cluster_count.astype(float)).fillna(0)
    cc_hard_dist_table = 1-cc_hard_sim_table
    Z = hclust.linkage(dist.squareform(np.array(cc_hard_dist_table)), method=hclust_linkage_method, metric=hclust_linkage_metric)
    cluster_map = hclust.fcluster(Z, k, criterion='maxclust')
    cluster_assign = pd.Series({cc_hard_dist_table.index[i]:cluster_map[i] for i in range(len(cc_hard_dist_table.index))}, name='CC Hard, k='+repr(k))
    # Save consensus clustering results
    if 'outdir' in save_args:
        if 'job_name' in save_args:
            save_cc_matrix_path = save_args['outdir']+str(save_args['job_name'])+'_cc_matrix.csv'
            save_clusters_path = save_args['outdir']+str(save_args['job_name'])+'_cluster_assignments.csv'
        else:
            save_cc_matrix_path = save_args['outdir']+'cc_matrix.csv'
            save_clusters_path = save_args['outdir']+'cluster_assignments.csv'
        cc_hard_sim_table.to_csv(save_cc_matrix_path)
        cluster_assign.to_csv(save_clusters_path)
    if verbose:
        print 'Hlist consensus constructed and sample clusters assigned'
    return cc_hard_sim_table, Z, cluster_assign

# Constructs Hlist object for consensus clustering functions if NBS iterations were run in parallel and outputs saved to a folder
def Hlist_constructor_from_folder(folder, ext='.csv'):
    co_clustering_results = [folder+fn for fn in os.listdir(folder) if fn.endswith(ext)]
    # Generate list of patient clusterings from netNMF
    Hlist = [pd.read_csv(fn, index_col=0) for fn in co_clustering_results]
    return Hlist    