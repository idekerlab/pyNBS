##############################################
# ---------- NBS Wrapper Function ---------- #
##############################################

from pyNBS import pyNBS_core as core
from pyNBS import network_propagation as prop
import networkx as nx
import numpy as np
import pandas as pd
import random

# Wrapper function to run a single instance of network-regularized NMF on given somatic mutation data and network
# sm_mat = binary mutation matrix of data to perform network NMF on
# options = dictionary of options to set for various parts of netNMF execuction
# propNet = NetworkX graph object to propagate binary mutations over
# regNet_glap = Pandas DataFrame graph laplacian of network from influence matrix of propNet
def NBS_single(sm_mat, options, propNet=None, propNet_kernel=None, regNet_glap=None, verbose=True, save_path=None):
    # Set default NBS netNMF options
    NBS_options = {'pats_subsample_p':0.8, 
                   'gene_subsample_p':0.8, 
                   'min_muts':10,
                   'prop_data':True, 
                   'prop_alpha':0.7, 
                   'prop_symmetric_norm':False, 
                   'qnorm_data':True,
                   'netNMF_k':4, 
                   'netNMF_gamma':200, 
                   'netNMF_niter':250, 
                   'netNMF_eps':1e-15, 
                   'netNMF_err_tol':1e-4, 
                   'netNMF_err_delta_tol':1e-4}

    # Update NBS netNMF options
    for option in options:
        NBS_options[option] = options[option]
    if verbose:
        print 'NBS options set:'
        for option in NBS_options:
            print '\t', option+':', NBS_options[option]
    
    # Check for correct input data
    if NBS_options['prop_data']:
        if type(propNet)!=nx.Graph:
            raise TypeError('Networkx graph object required for propNet')
    if (NBS_options['netNMF_gamma']!=0):
        if type(regNet_glap)!=pd.DataFrame:
            raise TypeError('netNMF regularization network laplacian (regNet_glap) must be given as Pandas DataFrame')

    # Subsample Data
    sm_mat_subsample = core.subsample_sm_mat(sm_mat, propNet=propNet, 
                                             pats_subsample_p=NBS_options['pats_subsample_p'], 
                                             gene_subsample_p=NBS_options['gene_subsample_p'], 
                                             min_muts=NBS_options['min_muts'])
    if verbose:
        print 'Somatic mutation data sub-sampling complete'

    # Propagate Data
    if NBS_options['prop_data']:
        if propNet_kernel is None:
            prop_sm_data = prop.network_propagation(propNet, sm_mat_subsample, 
                                                    symmetric_norm=NBS_options['prop_symmetric_norm'], 
                                                    alpha=NBS_options['prop_alpha'], verbose=verbose)
        else:
            prop_sm_data = prop.network_kernel_propagation(propNet, propNet_kernel, sm_mat_subsample, 
                                                           verbose=verbose, save_path=None)
        if verbose:
            print 'Somatic mutation data propagated'
    else:
        prop_sm_data = sm_mat_subsample
        print 'Somatic mutation data not propagated'

    # Quantile Normalize Data
    if NBS_options['qnorm_data']:
        prop_data_qnorm = core.qnorm(prop_sm_data)
        if verbose:
            print 'Somatic mutation data quantile normalized'
    else:
        prop_data_qnorm = prop_sm_data
        print 'Somatic mutation data not quantile normalized'

    # Prepare data for mixed netNMF function
    propNet_nodes = propNet.nodes()
    data_arr = np.array(prop_data_qnorm.T.ix[propNet_nodes])
    regNet_glap_arr = np.array(regNet_glap.ix[propNet_nodes][propNet_nodes])
    # Mixed netNMF Result
    W, H, numIter, finalResid = core.mixed_netNMF(data_arr, regNet_glap_arr, NBS_options['netNMF_k'], W_init=None, H_init=None, 
                                                  gamma=NBS_options['netNMF_gamma'], niter=NBS_options['netNMF_niter'], 
                                                  eps=NBS_options['netNMF_eps'], err_tol=NBS_options['netNMF_err_tol'], 
                                                  err_delta_tol=NBS_options['netNMF_err_delta_tol'], verbose=verbose, debug_mode=False)
    
    # Save netNMF Result
    H_df = pd.DataFrame(H.T, index=prop_data_qnorm.index)
    if save_path is not None:
        H_df.to_csv(save_path)
        if verbose:
            print 'netNMF result saved:', save_path
    return H_df