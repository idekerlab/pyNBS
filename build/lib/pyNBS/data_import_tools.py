###############################################
# ---------- Data Import Functions ---------- #
###############################################

import pandas as pd
import networkx as nx
import time
import os
import random

# Load network from file as unweighted network
# Can set delimiter, but default delimiter is tab
# Only will read edges as first two columns, all other columns will be ignored
# There are also options to shuffle the network to be loaded if desired (testing randomized network controls)
def load_network_file(network_file_path, delimiter='\t', degree_shuffle=False, label_shuffle=False, verbose=False):
    network = nx.read_edgelist(network_file_path, delimiter=delimiter, data=False)
    if verbose:
        print 'Network File Loaded:', network_file_path
    if degree_shuffle:
        network = degree_shuffNet(network, verbose=verbose)
    if label_shuffle:
        network = label_shuffNet(network, verbose=verbose)
    return network

# Shuffle network by preserving node-degree
def degree_shuffNet(network, verbose=False):
	shuff_time = time.time()
	edge_len=len(network.edges())
	shuff_net=network.copy()
	try:
		nx.double_edge_swap(shuff_net, nswap=edge_len, max_tries=edge_len*10)
	except:
		if verbose:
			print 'Note: Maximum number of swap attempts ('+repr(edge_len*10)+') exceeded before desired swaps achieved ('+repr(edge_len)+').'
	if verbose:
		# Evaluate Network Similarity
		shared_edges = len(set(network.edges()).intersection(set(shuff_net.edges())))
		print 'Network shuffled:', time.time()-shuff_time, 'seconds. Edge similarity:', shared_edges/float(edge_len)
	return shuff_net

# Shuffle network by permuting network node labels
def label_shuffNet(network, verbose=False):
    shuff_time = time.time()
    edge_len=len(network.edges())
    # Permute node labels
    network_nodes = network.nodes()
    shuff_nodes = list(network_nodes)
    for i in range(10):
        random.shuffle(shuff_nodes)
    network_relabel_map = {network_nodes[i]:shuff_nodes[i] for i in range(len(network_nodes))}    
    shuff_net = nx.relabel_nodes(network, network_relabel_map, copy=True)
    if verbose:
        # Evaluate Network Similarity
        shared_edges = len(set(network.edges()).intersection(set(shuff_net.edges())))
        print 'Network shuffled:', time.time()-shuff_time, 'seconds. Edge similarity:', shared_edges/float(edge_len)
    return shuff_net    	

# Filter extended sif file where all edges are weighted by a specific quantile
# Return the filtered network edge list and save it to a file if desired (for import by load_network_file)
# The input weighted network file may be any table format of edge list, but the columns for Node A, Node B, and weight must be specified
def filter_weighted_network(network_file_path, nodeA_col=0, nodeB_col=1, score_col=2, q=0.9, delimiter='\t', verbose=False, save_path=None):
    data = pd.read_csv(network_file_path, sep=delimiter, header=-1, low_memory=False)
    # Filter edges by score quantile
    q_score = data[score_col].quantile(q)
    if verbose:
        print str(round(q*100,2))+'%', 'score:', q_score
    data_filt = data[data[score_col]>q_score][data.columns[[nodeA_col, nodeB_col, score_col]]]
    data_filt.columns = ['nodeA', 'nodeB', 'edgeScore']
    if verbose:
        print data_filt.shape[0], '/', data.shape[0], 'edges retained'
    if save_path is not None:
        data_filt.to_csv(save_path, sep='\t', header=False, index=False)
    return data_filt

# Convert and save MAF from Broad Firehose
# Can produce 2 types of filetypes: 'matrix' or 'list', matrix is a full samples-by-genes binary csv, 'list' is a sparse representaiton of 'matrix'
# This is a conversion tool, so the result must be saved (most tools will require a path to a processed MAF file and load it separately)
# Gene naming can be 'Symbol' or 'Entrez'
def process_TCGA_MAF(maf_file, save_path, filetype='matrix', gene_naming='Symbol', verbose=False):
	loadtime = time.time()
	# Load MAF File
	TCGA_MAF = pd.read_csv(maf_file,sep='\t',low_memory=False)
	# Get all patient somatic mutation (sm) pairs from MAF file
	if gene_naming=='Entrez':
		TCGA_sm = TCGA_MAF.groupby(['Tumor_Sample_Barcode', 'Entrez_Gene_Id']).size()
	else:
		TCGA_sm = TCGA_MAF.groupby(['Tumor_Sample_Barcode', 'Hugo_Symbol']).size()
	# Turn somatic mutation data into binary matrix
	TCGA_sm_mat = TCGA_sm.unstack().fillna(0)
	TCGA_sm_mat = (TCGA_sm_mat>0).astype(int)
	# Trim TCGA barcodes
	TCGA_sm_mat.index = [pat[:12] for pat in TCGA_sm_mat.index]
	# Filter samples with duplicate IDs
	non_dup_IDs = list(TCGA_sm_mat.index.value_counts().index[TCGA_sm_mat.index.value_counts()==1])
	dup_IDs = list(TCGA_sm_mat.index.value_counts().index[TCGA_sm_mat.index.value_counts()>1])
	# Save file as binary matrix or sparse list
	if filetype=='list':
		# Now try to construct two-column/sparse representation of binary sm data
		# Get list of all patient somatic mutations
		index_list = list(TCGA_sm.index)
		# Filter list of patient somatic mutations of duplicate patient barcodes
		index_list_filt = [i for i in index_list if not any([True if barcode in i[0] else False for barcode in dup_IDs])]
		# Save patient somatic mutations list to file
		f = open(save_path, 'w')
		for sm in index_list_filt:
			f.write(sm[0][:12]+'\t'+sm[1]+'\n')
		f.close()
		if verbose:
			print 'Binary somatic mutations list saved'
	else:
		# Save non-duplicate patients' binary TCGA somatic mutation matrix to csv
		TCGA_sm_mat_filt = TCGA_sm_mat.ix[non_dup_IDs]
		# Remove all genes that have no more mutations after patient filtering
		nonempty_cols = [col for col in TCGA_sm_mat_filt.columns if not all(TCGA_sm_mat_filt[col]==0)]
		TCGA_sm_mat_filt2 = TCGA_sm_mat_filt[nonempty_cols]
		# Remove columns with bad names like '0'
		named_cols = [col for col in TCGA_sm_mat_filt.columns if col!='0']
		TCGA_sm_mat_filt3 = TCGA_sm_mat_filt2[nonempty_cols]
		TCGA_sm_mat_filt3.to_csv(save_path)
		if verbose:
			print 'Binary somatic mutation matrix saved'
	if verbose:
		print 'MAF file processed:', maf_file, round(time.time()-loadtime, 2), 'seconds.'
	return

# Load binary mutation data with 2 file types (filetype= 'matrix' or 'list')
# filetype=='matrix' is a csv or tsv style matrix with row and column headers, rows are samples/patients, columns are genes
# filetype=='list' is a 2 columns text file separated by the delimiter where 1st column is sample/patient, 2nd column is one gene mutated in that patient
# Line example in 'list' file: 'Patient ID','Gene Mutated'
def load_binary_mutation_data(filename, filetype='matrix', delimiter=',', verbose=False):
	if filetype=='list':
		f = open(filename)
		binary_mat_lines = f.read().splitlines()
		binary_mat_data = [(line.split('\t')[0], line.split('\t')[1]) for line in binary_mat_lines]
		binary_mat_index = pd.MultiIndex.from_tuples(binary_mat_data, names=['Tumor_Sample_Barcode', 'Hugo_Symbol'])
		binary_mat_2col = pd.DataFrame(1, index=binary_mat_index, columns=[0])[0]
		binary_mat = binary_mat_2col.unstack().fillna(0)
	else:
		binary_mat = pd.read_csv(filename, delimiter=delimiter, index_col=0).astype(int)
	if verbose:
	   print 'Binary Mutation Matrix Loaded:', filename
	return binary_mat