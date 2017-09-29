##################################################################
# Command line script to run Python Network Based Stratification #
##################################################################
# Inputs:
#  - Binary somatic mutation data for each patient
#  - Unweighted/Undirected molecular network

# Outputs:
#  - Patient co-clustering matrix
#  - Patient cluster assignments
#  - Patient co-clustering map (optional)

from pyNBS import data_import_tools as dit
from pyNBS import network_propagation as prop
from pyNBS import pyNBS_core as core
from pyNBS import pyNBS_single
from pyNBS import consensus_clustering as cc
from pyNBS import pyNBS_plotting as plot
import argparse
import os
import time
import pandas as pd
import numpy as np

# Valid file path check (Does not check file formatting, but checks if given path exists and is readable)
def valid_infile(in_file):
    if not os.path.isfile(in_file):
        raise argparse.ArgumentTypeError("{0} is not a valid input file path".format(in_file))  
    if os.access(in_file, os.R_OK):
        return in_file
    else:
        raise argparse.ArgumentTypeError("{0} is not a readable input file".format(in_file))

# Valid output directory path check (Checks if the output directory path can be found and written to by removing given filename from full path)
# Note: This uses '/' character for splitting pathnames on Linux and Mac OSX. The character may need to be changed to '\' for Windows executions
def valid_outfile(out_file):
    outdir = '/'.join(out_file.split('/')[:-1])
    if not os.path.isdir(outdir):
        raise argparse.ArgumentTypeError("{0} is not a valid output directory".format(outdir))
    if os.access(outdir, os.W_OK):
        return out_file
    else:
        raise argparse.ArgumentTypeError("{0} is not a writable output directory".format(outdir))

# Checking valid integer values (for all values that must be >0)
def positive_int(x):
    x = int(x)
    if x <= 0:
         raise argparse.ArgumentTypeError("%s must be a positive integer" % x)
    return x

# Checking valid alpha and p values (Range is 0.0-1.0 exclusive)
def restricted_float(x):
    x = float(x)
    if x <= 0.0 or x >= 1.0:
        raise arg    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform network based stratification of tumors by binary somatic mutation profile on given molecular network.')

    # Required files from user
    parser.add_argument("sm_data_file", type=valid_infile, 
        help='Path to binary mutation matrix file. May be a csv or 2-column list where each line is a sample and the gene mutated separated by a common delimiter.')    
    parser.add_argument("network_path", type=valid_infile, 
        help='Path to molecular network file. File must be table where each line is a gene interaction separated by a common delimiter and the first 2 columns represent interacting proteins.')

    # Other optional run parameters
    parser.add_argument('-v', '--verbose', default=False, action="store_true", required=False,
        help='Verbosity flag for reporting on patient similarity network construction steps.')

    # Parameters for loading somatic mutation data
    parser.add_argument('-mf', '--mut_filetype', type=str, default='list', choices=['matrix', 'list'], required=False,
        help='File structure of binary mutation data. 2 options: "matrix" (e.g. csv or tsv) or "list" (2-column list). Typically reading a "list" is faster.')
    parser.add_argument('-md', '--mut_filedelim', type=str, default='\t', required=False,
        help='Delimiter used in binary mutation file. Default is tab white space.')    

    # Parameters for loading molecular network
    parser.add_argument('-nd', '--net_filedelim', type=str, default='\t', required=False,
        help='Delimiter used in network file between columns. Default is tab white space.') 
    parser.add_argument('-sdp', '--degree_preserved_shuffle', default=False, action="store_true", required=False,
        help='Determination of whether or not to shuffle the network edges (while preserving node degree) when loading network.')    
    parser.add_argument('-snl', '--node_label_shuffle', default=False, action="store_true", required=False,
        help='Determination of whether or not to shuffle the network node labels (while preserving network topology) when loading network.')

    # Parameters for calculating regularization network
    parser.add_argument('-reg', '--regularize_network', action="store_true", required=False,
        help='Determination of whether or not to calculate influence matrix regularization network for regularized NMF step.')
    parser.add_argument('-g', '--gamma', type=float, default=0.01, required=False,
        help='Value of adjustment on propagation network graph laplacian to calculate influence matrix for (via Vandin 2011). Used if -reg is True')
    parser.add_argument('-kn', '--k_nearest_neighbors', type=positive_int, default=11, required=False,
        help='Number of nearest neighbors to add to the regularization network during construction. Used if -reg is True.')
    parser.add_argument('-s_knngl', '--save_knn_glap', type=valid_outfile, default=None, required=False,
        help='File path of where to save graph laplacian for k-nearest-neighbor network constructed from propagation network influence matrix. No path given as default, automatically saves pandas hdf file if file path given.')
    parser.add_argument('-rngl', '--regularization_network_graph_laplacian_file', type=valid_infile, default=None, required=False,
        help='Path to regularization network graph laplacian matrix if previously calculated. Required if -reg is False.')

    # Parameters for sub-sampling data
    parser.add_argument('-n', '--niter', type=positive_int, default=1000, required=False,
        help='Number of iterations to perform sub-sampling and network-regularized NMF before consensus clustering.')
    parser.add_argument('-p', '--pats_subsample_p', type=float, default=0.8, required=False,
        help='Proportion of samples to sub-sample')
    parser.add_argument('-q', '--gene_subsample_q', type=float, default=0.8, required=False,
        help='Proportion of mutated genes to sub-sample')    
    parser.add_argument('-mm', '--min_muts', type=positive_int, default=10, required=False,
        help='Minimum number of mutations for a sample to contain after sub-sampling to be considered for further analysis.')

    # Parameters for network propagation
    parser.add_argument('-prop', '--propagate_data', type=bool, default=True, required=False,
        help='Determination of whether or not to propagate sub-sampled binary mutation data over given molecular network.')
    parser.add_argument('-cpk', '--calculate_propagation_kernel', default=False, action="store_true", required=False,
        help='Determination of whether or not to pre-calculate network kernel for network propagation. Highly recommended if no network kernel file is given already and niter > 10.') 
    parser.add_argument('-kernel', '--propagation_kernel_file', type=valid_infile, default=None, required=False,
        help='Path to pre-calculated propagation kernel of network. This will save time in the propagation step.')    
    parser.add_argument('-a', '--alpha', type=restricted_float, default=0.7, required=False,
        help='Propagation constant to use in the propagation of mutations over molecular network. Range is 0.0-1.0 exclusive.')
    parser.add_argument('-norm', '--symmetric_network_normalization', type=bool, default=False, required=False,
        help='Network degree normalization method for random walk-propagation. See network_propagation() comments for more details.')

    # Parameters for data quantile normalization
    parser.add_argument("-qnorm", "--quantile_normalize_data", type=bool, default=True, required=False,
        help='Determination of whether or not to qunatile normalize mutation profiles.')

    # Parameters for network-regularized NMF (netNMF) decomposition
    parser.add_argument('-k', '--K', type=positive_int, default=4, required=True,
        help='Number of components to decompose patient mutation data into. Same as the number of clusters of patients to separate data into.')
    parser.add_argument('-netNMF_g', '--netNMF_gamma', type=positive_int, default=200, required=False,
        help='Regularization constant to scale network regularization term in netNMF.')
    parser.add_argument('-netNMF_ug', '--netNMF_update_gamma', type=bool, default=False, required=False,
        help='Determination of whether or not to constantly update regularization constant based on balance between reconstruction error and regularization term.')
    parser.add_argument('-netNMF_gf', '--netNMF_gamma_factor', type=positive_int, default=1, required=False,
        help='Scaling factor for regularization constant updates if -ug is True.')
    parser.add_argument('-e', '--netNMF_eps', type=float, default=1e-15, required=False,
        help='Epsilon error value to adjust 0 values during multiplicative matrix updates in netNMF')
    parser.add_argument('-netNMF_n', '--netNMF_niter', type=positive_int, default=250, required=False,
        help='Maximum umber of multiplicative updates to perform within network-regularized NMF if result does not converge.')    
    parser.add_argument('-tol', '--netNMF_err_tol', type=float, default=1e-4, required=False,
        help='Minimum error tolerance for matrix reconstruction of original data for convergence.')
    parser.add_argument('-tol_d', '--netNMF_err_delta_tol', type=float, default=1e-4, required=False,
        help='Minimum error tolerance for l2 norm of difference in matrix reconstructions between iterations of netNMF for convergence.')
    parser.add_argument('-s_H', '--save_H', type=valid_outfile, default=None, required=False,
        help='File path of where to save decomposed patient profiles. No path given as default, automatically saves csv file if file path given.')    

    # Parameters for consensus clustering
    parser.add_argument('-cc', '--consensus_cluster', type=bool, default=True, required=False,
        help='Determination of whether or not to perform consensus clustering on decompositions of patient profiles.')
    parser.add_argument('-ac', '--assign_clusters', type=bool, default=True, required=False,
        help='Determination of whether or not to assign numerical clusters to patients based on consensus clustering of patient profiles.')    
    parser.add_argument('-s_cc', '--save_co_cluster_matrix', type=valid_outfile, default=None, required=False,
        help='File path of where to save patient co-clustering matrix. No path given as default, automatically saves csv file if file path given.')    
    parser.add_argument('-s_ca', '--save_cluster_assignments', type=valid_outfile, default=None, required=False,
        help='File path of where to save patient cluster assignments. No path given as default, automatically saves csv file if file path given.')    

    # Parameters for consensus clustering map plotting
    parser.add_argument('-map', '--plot_co_cluster_map', type=bool, default=True, required=False,
        help='Determination of whether or not to plot the co-clustering matrix. Requires consensus clustering and cluster assignments.')
    parser.add_argument('-t', '--plot_title', type=str, default=None, required=False,
        help='Title of co-clustering matrix map if desired.')    
    parser.add_argument('-s_map', '--save_co_cluster_map', type=valid_outfile, default=None, required=False,
        help='File path of where to save co-clustering matrix plot. No path given as default, automatically saves pdf file if file path given.')

    # KM plot options (this block is edited and differ from what is on Github)
    parser.add_argument('-km_plot', '--Kaplan_Meier_plot', type=bool, default=None, required=False,
        help='Determination of whether or not to save Kaplan-Meier plot.')
    parser.add_argument('-clin', '--clinical_data', type=str, default=None, required=False,
        help='File path of patient clinical data. This is required if Kaplan_Meier_plot is True.')
    parser.add_argument('-km_title', '--Kaplan_Meier_plot_title', type=str, default=None, required=False,
        help='Title of Kaplan-Meier plot if desired.')
    parser.add_argument('-s_km', '--save_Kaplan_Meier_plot', type=valid_outfile, default=None, required=False,
        help='File path of where to save Kaplan-Meier plot. No path given as default, automatically saves pdf file if file path given.') 



    # Error checking parameter conditions
    args = parser.parse_args()
    # 1. Regularization network path must be provided if not regularizing given network
    if (not args.regularize_network) & (args.regularization_network_graph_laplacian_file is None):
        parser.error('Regularization network path is required if not regularizing given molecular network.')
    # 2. Consensus clustering required to assign clusters
    if (not args.consensus_cluster) & (args.assign_clusters):
        parser.error('Consensus clustering required to assign patient clusters.')
    # 3. Must be performing consensus clustering and cluster assignments if we want to plot co-clustering matrix
    if args.plot_co_cluster_map:
        if (not args.consensus_cluster) | (not args.assign_clusters):
            parser.error('Consensus clustering and cluster assignments must be done to plot co-clustering matrix.')

    # Begin NBS
    if args.verbose:
        print 'Loading Data'
    # Load somatic mutation data
    sm_mat = dit.load_binary_mutation_data(args.sm_data_file, filetype=args.mut_filetype, delimiter=args.mut_filedelim, verbose=args.verbose)
    # Load network
    network = dit.load_network_file(args.network_path, delimiter=args.net_filedelim, degree_shuffle=args.degree_preserved_shuffle, 
                                    label_shuffle=args.node_label_shuffle, verbose=args.verbose)
    
    # Get knnGlap
    if args.regularize_network:
        knnGlap = core.network_inf_KNN_glap(network, gamma=args.gamma, kn=args.k_nearest_neighbors, verbose=args.verbose, save_path=args.save_knn_glap)
    else:
        # Load propatagion kernel
        if args.regularization_network_graph_laplacian_file.endswith('.hdf'):
            knnGlap = pd.read_hdf(args.regularization_network_graph_laplacian_file)
        else:
            knnGlap = pd.read_csv(args.regularization_network_graph_laplacian_file)
        if args.verbose:
            print 'Pre-calculated regularization network graph laplacian loaded'
    
    # Get network propagation kernel
    if args.propagation_kernel_file is not None:
        # Load propagation kernel
        if args.propagation_kernel_file.endswith('.hdf'):
            kernel = pd.read_hdf(args.propagation_kernel_file)
        else:
            kernel = pd.read_csv(args.propagation_kernel_file)
        if args.verbose:
            print 'Pre-calculated network kernel loaded'
    else:
        if args.calculate_propagation_kernel:
            # Calculate propagation kernel by propagating identity matrix of network
            network_nodes = network.nodes()
            network_I = pd.DataFrame(np.identity(len(network_nodes)), index=network_nodes, columns=network_nodes)
            kernel = prop.network_propagation(network, network_I, args.alpha, verbose=True)  
            if args.verbose:
                print 'Network kernel calculated'
        else:
            kernel = None
            if args.verbose:
                print 'No network kernel established'

    # Construct options dictionary for decomposition
    NBS_options = {'pats_subsample_p' : args.pats_subsample_p, 
                   'gene_subsample_p' : args.gene_subsample_q, 
                   'min_muts' : args.min_muts,
                   'prop_data' : args.propagate_data, 
                   'prop_alpha' : args.alpha, 
                   'prop_symmetric_norm' : args.symmetric_network_normalization, 
                   'qnorm_data' : args.quantile_normalize_data,
                   'netNMF_k' : args.K, 
                   'netNMF_gamma' : args.netNMF_gamma, 
                   'netNMF_update_gamma' : args.netNMF_update_gamma, 
                   'netNMF_gamma_factor' : args.netNMF_gamma_factor,
                   'netNMF_niter' : args.netNMF_niter, 
                   'netNMF_eps' : args.netNMF_eps, 
                   'netNMF_err_tol' : args.netNMF_err_tol, 
                   'netNMF_err_delta_tol' : args.netNMF_err_delta_tol}

    # Sub-sampling and netNMF decomposition
    Hlist = []
    for i in range(args.niter):
        netNMF_time = time.time()
        Hlist.append(pyNBS_single.NBS_single(sm_mat, NBS_options, propNet=network, propNet_kernel=kernel, 
                                             regNet_glap=knnGlap, verbose=False, save_path=args.save_H))
        if args.verbose:
            print 'NBS iteration:', i+1, 'complete:', time.time()-netNMF_time, 'seconds'

    # Consensus Clustering
    if args.consensus_cluster:
        NBS_cc_table, NBS_cc_linkage, NBS_cluster_assign = cc.consensus_hclust_hard(Hlist, args.K, assign_cluster=args.assign_clusters)
        if args.verbose:
            print 'Consensus Clustering complete'        
        if args.save_co_cluster_matrix is not None:
            NBS_cc_table.to_csv(args.save_co_cluster_matrix)
            if args.verbose:
                print 'Co-clustering matrix saved'
        if args.save_cluster_assignments is not None:
            NBS_cluster_assign.to_csv(args.save_cluster_assignments)
            if args.verbose:
                print 'Cluster assignments saved'

    # Plot Consensus Cluster Map
    if args.plot_co_cluster_map:
        NBS_cluster_assign_cmap = plot.cluster_color_assign(NBS_cluster_assign, name='Cluster Assignment')
        plot.plot_cc_map(NBS_cc_table, NBS_cc_linkage, title=args.plot_title, row_color_map=None, 
                         col_color_map=NBS_cluster_assign_cmap, save_path=args.save_co_cluster_map)
        if args.verbose:
            print 'Consensus Clustering map saved'


    # Kaplan-Meier Plot (This block of code is edited and differs from what is on Github)
    if args.Kaplan_Meier_plot:
        plot.cluster_KMplot(NBS_cluster_assign,clin_data_fn=args.clinical_data,title=args.Kaplan_Meier_plot_title, save_path = args.save_Kaplan_Meier_plot)
