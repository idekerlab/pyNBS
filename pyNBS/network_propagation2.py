#######################################################
# ---------- Network Propagation Functions ---------- #
#######################################################

import networkx as nx
import time
import numpy as np
import scipy
import pandas as pd

# Normalize network (or network subgraph) for random walk propagation
# If symmetric norm is used then the adjacency matrix is normalized as D^-0.5 * A * D^-0.5
# Otherwise the network is normalized as A*D-1
# Where D is the diagonalized degree sum of the adjacency matrix A
def normalize_network(network, symmetric_norm=False):
    adj_mat = nx.adjacency_matrix(network)
    adj_array = np.array(adj_mat.todense())
    if symmetric_norm:
        D = np.diag(1/np.sqrt(sum(adj_array)))
        adj_array_norm = np.dot(np.dot(D, adj_array), D)
    else:
        degree = np.sum(adj_mat, axis=1)
        adj_array_norm = adj_mat*1.0/degree
    return adj_array_norm
# Note about normalizing by degree, if multiply by degree_norm_array first (D^-1 * A), then do not need to return
# transposed adjacency array, it is already in the correct orientation

# Closed form random-walk propagation (as seen in HotNet2) for each subgraph: Ft = (1-alpha)*Fo * (I-alpha*norm_adj_mat)^-1
# Concatenate to previous set of subgraphs
def fast_random_walk(alpha, binary_mat, subgraph_norm, prop_data):
	term1=(1-alpha)*binary_mat
	term2=np.identity(binary_mat.shape[1])-alpha*subgraph_norm
	term2_inv = np.linalg.inv(term2)
	subgraph_prop = np.dot(term1, term2_inv)
	return np.concatenate((prop_data, subgraph_prop), axis=1)

# Wrapper for random walk propagation of full network by subgraphs
# Implementation is based on the closed form of the random walk model over networks presented by the HotNet2 paper
def network_propagation(network, binary_matrix, alpha, symmetric_norm=False, verbose=False, save_path=None):
	starttime=time.time()
	if verbose:
		print 'Performing network propagation with alpha:', alpha
	# Separate network into connected components and calculate propagation values of each sub-sample on each connected component
	subgraphs = list(nx.connected_component_subgraphs(network))
	# Initialize propagation results by propagating first subgraph
	subgraph = subgraphs[0]
	subgraph_nodes = subgraph.nodes()
	prop_data_node_order = list(subgraph_nodes)
	binary_matrix_filt = np.array(binary_matrix.T.ix[subgraph_nodes].fillna(0).T)
	subgraph_norm = normalize_network(subgraph, symmetric_norm=symmetric_norm)
	prop_data_empty = np.zeros((binary_matrix_filt.shape[0], 1))
	prop_data = fast_random_walk(alpha, binary_matrix_filt, subgraph_norm, prop_data_empty)
	# Get propagated results for remaining subgraphs
	for subgraph in subgraphs[1:]:
		subgraph_nodes = subgraph.nodes()
		prop_data_node_order = prop_data_node_order + subgraph_nodes
		binary_matrix_filt = np.array(binary_matrix.T.ix[subgraph_nodes].fillna(0).T)
		subgraph_norm = normalize_network(subgraph, symmetric_norm=symmetric_norm)
		prop_data = fast_random_walk(alpha, binary_matrix_filt, subgraph_norm, prop_data)
	# Return propagated result as dataframe
	prop_data_df = pd.DataFrame(data=prop_data[:,1:], index = binary_matrix.index, columns=prop_data_node_order)
	if save_path is not None:
		if save_path.endswith('.hdf'):
			hdf_key = save_path.split('/')[-1].split('.')[0]
			prop_data_df.to_hdf(save_path, key=hdf_key, mode='w')
		else:
			prop_data_df.to_csv(save_path)
	if verbose:
		print 'Network Propagation Complete:', time.time()-starttime, 'seconds'				
	return prop_data_df

# Wrapper for propagating binary mutation matrix over network by subgraph given network propagation kernel
# The network propagation kernel can be pre-computed using the network_propagation function and a identity matrix data frame of the network
# Pre-calculating the kernel for many runs of NBS saves a significant amount of time
def network_kernel_propagation(network, network_kernel, binary_matrix, verbose=False, save_path=None):
    starttime=time.time()
    if verbose:
        print 'Performing network propagation with network kernel'
    # Separate network into connected components and calculate propagation values of each sub-sample on each connected component
    subgraph_nodelists = list(nx.connected_components(network))
    # Initialize propagation results by propagating first subgraph
    prop_nodelist = list(subgraph_nodelists[0])
    prop_data = np.dot(binary_matrix[prop_nodelist], network_kernel.ix[prop_nodelist][prop_nodelist])
    # Get propagated results for remaining subgraphs
    for nodelist in subgraph_nodelists[1:]:
        subgraph_nodes = list(nodelist)
        prop_nodelist = prop_nodelist + subgraph_nodes
        subgraph_prop_data = np.dot(binary_matrix[subgraph_nodes], network_kernel.ix[subgraph_nodes][subgraph_nodes])
        prop_data = np.concatenate((prop_data, subgraph_prop_data), axis=1)
    # Return propagated result as dataframe
    prop_data_df = pd.DataFrame(data=prop_data, index = binary_matrix.index, columns=prop_nodelist)
    if save_path is not None:
        if save_path.endswith('.hdf'):
            hdf_key = save_path.split('/')[-1].split('.')[0]
            prop_data_df.to_hdf(save_path, key=hdf_key, mode='w')
        else:
            prop_data_df.to_csv(save_path)
    if verbose:
        print 'Network Propagation Complete:', time.time()-starttime, 'seconds'				
    return prop_data_df    		