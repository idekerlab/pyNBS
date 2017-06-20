############################################################
# ---------- NBS Consensus Clustering Functions ---------- #
############################################################

import os
import pandas as pd
import numpy as np
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hclust

# Constructs Hlist object for consensus clustering functions if NBS iterations were run in parallel and outputs saved to a folder
def Hlist_constructor_from_folder(folder, ext='.csv', normalize_H=False, verbose=False):
    co_clustering_results = [folder+fn for fn in os.listdir(folder) if fn.endswith(ext)]
    # Generate list of patient clusterings from netNMF
    Hlist = [pd.read_csv(fn, index_col=0) for fn in co_clustering_results]
    # Normalize H matrices if needed (to make columns comparable if not already done in decomposition)
    if normalize_H:
        Hlist_norm = []
        for H in Hlist:
            H_norm = np.dot(H,np.diag(1/H.sum()))
            Hlist_norm.append(pd.DataFrame(H_norm, index=H.index))
        if verbose:
            print 'Hlist constructed and normalized'
        return Hlist_norm
    else:
        if verbose:
            print 'Hlist constructed'
        return Hlist

# Takes a list of 'H' (patient-by-k) dataframes and performs 'hard' consensus clustering
# Using hierarchical clustering and average linkage
# Returns similarity table (distance is 1-similarity) and linkage map of patients
# Also returns cluster assignment map of patients if wanted
def consensus_hclust_hard(Hlist, k, assign_cluster=False, verbose=False):
    # Generate patient list
    pat_list = set()
    for H in Hlist:
        pat_list = pat_list.union(set(H.index))
    pat_list = sorted(list(pat_list))
    if verbose:
        print 'Constructing Hlist:', len(Hlist), 'cluster matrices', len(pat_list), 'samples'

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
    Z = hclust.linkage(dist.squareform(np.array(cc_hard_dist_table)), method='average')
    if assign_cluster:
        cluster_map = hclust.fcluster(Z, k, criterion='maxclust')
        cluster_assign = pd.Series({cc_hard_dist_table.index[i]:cluster_map[i] for i in range(len(cc_hard_dist_table.index))}, name='CC Hard, k='+repr(k))
        if verbose:
            print 'Hlist consensus constructed and sample clusters assigned'
        return cc_hard_sim_table, Z, cluster_assign
    else:
        if verbose:
            print 'Hlist consensus constructed and sample clusters assigned'
        return cc_hard_sim_table, Z