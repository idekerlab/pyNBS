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
# regNet_glap = Pandas DataFrame graph laplacian of network from influence matrix of propNet
# ^ This is a required parameter! We cannot perform network-regularized NMF without the regularization network
# sm_mat = binary mutation matrix of data to perform network NMF on
# propNet = NetworkX graph object to propagate binary mutations over
# kwargs = dictionary of parameters to set for various parts of netNMF setup and execuction
def NBS_single(sm_mat, regNet_glap, propNet=None, propNet_kernel=None, 
    prop_data=True, qnorm_data=True, k=3, verbose=False, **kwargs):
    # Check for correct input data
    if type(sm_mat)!=pd.DataFrame:
        raise TypeError('Somatic mutation data must be given as Pandas DataFrame')
    if propNet is not None:
        if type(propNet)!=nx.Graph:
            raise TypeError('Networkx graph object required for propNet')
    if regNet_glap is not None:
        if type(regNet_glap)!=pd.DataFrame:
            raise TypeError('netNMF regularization network laplacian (regNet_glap) must be given as Pandas DataFrame')

    # Load or set subsampling parameters
    pats_subsample_p, gene_subsample_p, min_muts = 0.8, 0.8, 10
    if 'pats_subsample_p' in kwargs:
        pats_subsample_p = float(kwargs['pats_subsample_p'])
    if 'gene_subsample_p' in kwargs:
        gene_subsample_p = float(kwargs['gene_subsample_p'])
    if 'min_muts' in kwargs:
        min_muts = int(kwargs['min_muts'])  

    # Subsample Data
    sm_mat_subsample = core.subsample_sm_mat(sm_mat, propNet=propNet, 
        pats_subsample_p=pats_subsample_p, gene_subsample_p=gene_subsample_p, min_muts=min_muts)
    if verbose:
        print 'Somatic mutation data sub-sampling complete'

    # Throw exception if subsampling returned empty dataframe
    if sm_mat_subsample.shape[0]==0:
        raise ValueError('Subsampled somatic mutation matrix contains no patients.')

    # Propagate data if network object is provided
    if prop_data is not None:
        # Determine if propagation is can be based on pre-computed propagation kernel
        if propNet_kernel is None:
            # If kernel is not given and some propagation parameters are given in kwargs, set propagation parameters
            # Otherwise set default values
            alpha, symmetric_norm, save_prop = 0.7, False, False
            if 'alpha' in kwargs:
                alpha = float(kwargs['alpha'])
            if 'symmetric_norm' in kwargs:
                symmetric_norm = bool(kwargs['symmetric_norm'])
            if 'save_prop' in kwargs:
                save_prop = bool(kwargs['save_prop'])
            # Save propagation step data if desired (indicated in kwargs)
            if save_prop:
                prop_sm_data = prop.network_propagation(propNet, sm_mat_subsample, alpha=alpha, symmetric_norm=symmetric_norm, **kwargs)
            else:
                prop_sm_data = prop.network_propagation(propNet, sm_mat_subsample, alpha=alpha, symmetric_norm=symmetric_norm)
        else:
            # Save propagation step data if desired (indicated in kwargs)
            save_prop = False
            if 'save_prop' in kwargs:
                save_prop = bool(kwargs['save_prop'])
            if save_prop:
                prop_sm_data = prop.network_kernel_propagation(propNet, propNet_kernel, sm_mat_subsample, **kwargs)
            else:
                prop_sm_data = prop.network_kernel_propagation(propNet, propNet_kernel, sm_mat_subsample)
        if verbose:
            print 'Somatic mutation data propagated'
    else:
        prop_sm_data = sm_mat_subsample
        if verbose:
          print 'Somatic mutation data not propagated'

    # Quantile Normalize Data
    qnorm_data = True
    if 'qnorm_data' in kwargs:
        qnorm_data = bool(kwargs['qnorm_data'])
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

    # Set netNMF parameters from kwargs if given, otherwise use defaults
    netNMF_gamma, netNMF_maxiter, netNMF_verbose = 200, 250, False
    netNMF_eps, netNMF_err_tol, netNMF_err_delta_tol = 1e-15, 1e-4, 1e-4
    if 'netNMF_gamma' in kwargs:
        netNMF_gamma = float(kwargs['netNMF_gamma'])
    if 'netNMF_maxiter' in kwargs:
        netNMF_maxiter = int(kwargs['netNMF_maxiter'])
    if 'netNMF_eps' in kwargs:
        netNMF_eps = float(kwargs['netNMF_eps'])
    if 'netNMF_err_tol' in kwargs:
        netNMF_err_tol = float(kwargs['netNMF_err_tol'])
    if 'netNMF_err_delta_tol' in kwargs:
        netNMF_err_delta_tol = float(kwargs['netNMF_err_delta_tol']) 
    if 'netNMF_verbose' in kwargs:
        netNMF_verbose = bool(kwargs['netNMF_verbose'])

    # Mixed netNMF Result
    W, H, numIter, finalResid = core.mixed_netNMF(data_arr, regNet_glap_arr, k=k, 
        gamma=netNMF_gamma, maxiter=netNMF_maxiter, eps=netNMF_eps, 
        err_tol=netNMF_err_tol, err_delta_tol=netNMF_err_delta_tol, verbose=netNMF_verbose)
    
    # Return netNMF result (dimension-reduced propagated patient profiles)
    H_df = pd.DataFrame(H.T, index=prop_data_qnorm.index)

    # Save netNMF result
    # Saving the propagation result
    if 'outdir' in kwargs:
        if 'job_name' in kwargs:
            if 'iteration_label' in kwargs:
                save_path = kwargs['outdir']+str(kwargs['job_name'])+'_H_'+str(kwargs['iteration_label'])+'.csv'
            else:
                save_path = kwargs['outdir']+str(kwargs['job_name'])+'_H.csv'
        else:
            if 'iteration_label' in kwargs:
                save_path = kwargs['outdir']+'H_'+str(kwargs['iteration_label'])+'.csv'
            else:
                save_path = kwargs['outdir']+'H.csv'
        H_df.to_csv(save_path)
        if verbose:
            print 'H matrix saved:', save_path
    else:
        pass
    if verbose:
        print 'pyNBS iteration complete'
    return H_df