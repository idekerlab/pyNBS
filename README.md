# pyNBS

This is a python package replicating the network-based stratification algorithm used in the Nature Methods Hofree et al 2013 paper.

To install:  
1. Clone the repository
2. cd to new respository
3. Execute following command:  
python setup.py install  

# Running the code 
There are two ways to execute the NBS code:
1. From a Jupyter Notebook. See the 'OV_run_pyNBS.ipynb' or 'UT_run_pyNBS' jupyter notebook for details
2. From the command line. See below instruction and 'run_NBS.py' for details
## Running from command line
The table below shows the required input to run pyNBS from terminal

Reuired Input|Help
-|-
sm_data_file|Path to binary mutation matrix file. May be a csv or 2-column list where each line is a sample and the gene mutated separated by a common delimiter.
network_path|Path to molecular network file. File must be table where each line is a gene interaction separated by a common delimiter and the first 2 columns represent interacting proteins.

The table below shows the options for running pyNBS

Option|Flag|Default|Help
-|-|-|-
verbose|-v|False|Verbosity flag for reporting on patient similarity network construction steps.
mut_filetype|-mf|matrix|File structure of binary mutation data. 2 options: "matrix" (e.g. csv or tsv) or "list" (2-column list). Typically reading a "list" is faster.
mut_filedelim|-md|'\t'|Delimiter used in binary mutation file. 
net_filedelim|-nd|'\t'|Delimiter used in network file between columns. 
degree_preserved_shuffle|sdp|False|Determination of whether or not to shuffle the network edges (while preserving node degree) when loading network.
node_label_shuffle|-snl|False|Determination of whether or not to shuffle the network node labels (while preserving network topology) when loading network.
regularize_network|-reg|False|Determination of whether or not to calculate influence matrix regularization network for regularized NMF step.
reg_net_gamma|-g|0.01|Value of adjustment on propagation network graph laplacian to calculate influence matrix for (via Vandin 2011).
k_nearest_neighbors|kn|11|Number of nearest neighbors to add to the regularization network during construction.
save_knn_glap|-s_knngl|None|File path of where to save graph laplacian for k-nearest-neighbor network constructed from propagation network influence matrix. No path given as default, automatically saves pandas hdf file if file path given.
regularization_network_graph_laplacian_file|-rngl|None|Path to regularization network graph laplacian matrix if previously calculated. Required if 'regularize_network' is False.
niter|-n|1000|Number of iterations to perform sub-sampling and network-regularized NMF before consensus clustering.
pats_subsample_p|-p|0.8|Proportion of samples to sub-sample
gene_subsample_q|-q|0.8|Proportion of mutated genes to sub-sample
min_muts|-nm|10|Minimum number of mutations for a sample to contain after sub-sampling to be considered for further analysis.
propagate_data|-prop|True|Determination of whether or not to propagate sub-sampled binary mutation data over given molecular network.
calculate_propagation_kernel|-cpk|False|Determination of whether or not to pre-calculate network kernel for network propagation. Highly recommended if no network kernel file is given already and niter > 10.
propagation_kernel_file|-kernel|None|Path to pre-calculated propagation kernel of network. This will save time in the propagation step.
alpha|-a|0.7|Propagation constant to use in the propagation of mutations over molecular network. Range is 0.0-1.0 exclusive.
symmetric_network_normalization|-norm|False|Network degree normalization method for random walk-propagation.
quantile_normalize_data|-qnorm|True|Determination of whether or not to qunatile normalize mutation profiles.
K|-k|4 |Number of components to decompose patient mutation data into. Same as the number of clusters of patients to separate data into. 
netNMF_gamma |-netNMF_g| 200 |Regularization constant to scale network regularization term in netNMF. 
netNMF_update_gamma |-netNMF_ug|False |Determination of whether or not to constantly update regularization constant based on balance between reconstruction error and regularization term.
netNMF_gamma_factor |-netNMF_gf |1 |Scaling factor for regularization constant updates if 'netNMF_update_gamma' is True.
netNMF_eps|-et |1e-15  |Epsilon error value to adjust 0 values during multiplicative matrix updates in netNMF 
netNMF_niter  |-netNMF_n  |250 |Maximum umber of multiplicative updates to perform within network-regularized NMF if result does not converge. 
netNMF_err_tol  |-tol  |1e-4  |Minimum error tolerance for matrix reconstruction of original data for convergence. 
netNMF_err_delta_tol  |-tol_d |1e-4  |Minimum error tolerance for l2 norm of difference in matrix reconstructions between iterations of netNMF for convergence. 
save_H|-s_H|None|File path of where to save decomposed patient profiles. No path given as default, automatically saves csv file if file path given.
consensus_cluster|-cc|True|Determination of whether or not to perform consensus clustering on decompositions of patient profiles.
assign_clusters|-ac|True|Determination of whether or not to assign numerical clusters to patients based on consensus clustering of patient profiles.
save_co_cluster_matrix|-s_cc|None|File path of where to save patient co-clustering matrix. No path given as default, automatically saves csv file if file path given.
save_cluster_assignments|-s_ca|None|File path of where to save patient cluster assignments. No path given as default, automatically saves csv file if file path given.
plot_co_cluster_map|-map|True|Determination of whether or not to plot the co-clustering matrix. Requires consensus clustering and cluster assignments.
plot_title|-t|None|Title of co-clustering matrix map if desired.
save_co_cluster_map|-s_map|None|File path of where to save co-clustering matrix plot. No path given as default, automatically saves pdf file if file path given.

