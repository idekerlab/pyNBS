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
# propNet = NetworkX graph object to propagate binary mutations over
# regNet_glap = Pandas DataFrame graph laplacian of network from influence matrix of propNet
# params = dictionary of parameters to set for various parts of netNMF setup and execuction
def NBS_single(sm_mat, propNet=None, propNet_kernel=None, regNet_glap=None, params=None):
    # Set default NBS internal step parameters
    verbose = False
    prop_data = True
    qnorm_data = True
    netNMF_gamma = 200
    save_H = False
    iteration = '1'
    if (params is not None) and (type(params)==dict):
        if 'verbose' in params:
            verbose = bool(params['verbose'])
        if 'prop_data' in params:
            prop_data = bool(params['prop_data'])
        if 'qnorm_data' in params:
            qnorm_data = bool(params['qnorm_data'])
        if 'netNMF_gamma' in params:
            netNMF_gamma = float(params['netNMF_gamma'])            
        if 'save_H' in params:
            save_H = bool(params['save_H'])
        if 'pyNBS_iteration' in params:
            iteration = str(params['pyNBS_iteration']) 

    # Check for correct input data
    if type(sm_mat)!=pd.DataFrame:
        raise TypeError('Somatic mutation data must be given as Pandas DataFrame')
    if prop_data:
        if type(propNet)!=nx.Graph:
            raise TypeError('Networkx graph object required for propNet')
    if (netNMF_gamma!=0):
        if type(regNet_glap)!=pd.DataFrame:
            raise TypeError('netNMF regularization network laplacian (regNet_glap) must be given as Pandas DataFrame')

    # Subsample Data
    sm_mat_subsample = core.subsample_sm_mat(sm_mat, propNet=propNet, params=params)
    if verbose:
        print 'Somatic mutation data sub-sampling complete'

    # Propagate Data
    # Throw exception if subsampling returned empty dataframe
    if sm_mat_subsample.shape[0]==0:
        raise ValueError('Subsampled somatic mutation matrix contains no patients.')
    if prop_data:
        if propNet_kernel is None:
            prop_sm_data = prop.network_propagation(propNet, sm_mat_subsample, params=params)
        else:
            prop_sm_data = prop.network_kernel_propagation(propNet, propNet_kernel, sm_mat_subsample, 
                                                           params=params)
        if verbose:
            print 'Somatic mutation data propagated'
    else:
        prop_sm_data = sm_mat_subsample
        if verbose:
          print 'Somatic mutation data not propagated'

    # Quantile Normalize Data
    if qnorm_data:
        prop_data_qnorm = core.qnorm(prop_sm_data)
        if verbose:
            print 'Somatic mutation data quantile normalized'
    else:
        prop_data_qnorm = prop_sm_data
        if verbose:
          print 'Somatic mutation data not quantile normalized'

    # Prepare data for mixed netNMF function (align propagated profile columns with regularization network laplacian rows)
    propNet_nodes = propNet.nodes()
    data_arr = np.array(prop_data_qnorm.T.ix[propNet_nodes])
    regNet_glap_arr = np.array(regNet_glap.ix[propNet_nodes][propNet_nodes])
    # Mixed netNMF Result
    W, H, numIter, finalResid = core.mixed_netNMF(data_arr, regNet_glap_arr, params=params)
    
    # Save netNMF Result
    H_df = pd.DataFrame(H.T, index=prop_data_qnorm.index)
    if save_H:
        save_H_path = params['outdir']+params['job_name']+'_H_'+iteration+'.csv'
        H_df.to_csv(save_H_path)
        if verbose:
            print 'netNMF result saved:', save_path
    return H_df