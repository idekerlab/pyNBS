#! /bin/bash

# Valid values for 'cancer_type' for these examples are: 'BLCA', 'COAD', 'HNSC', 'OV', 'UCEC'
# Use one of the above codes when calling this script as the first parameter
# For more information on pyNBS options, please visit the GitHub Documentation Wiki for more detailed documentation
cancer_type=$1
echo
echo 'Running pyNBS command line example for: '$cancer_type
# The second parameter is the number of clusters to separate all samples into (typically a small integer)
k=$2
echo 'Number of clusters: '$k
# The third parameter is the number of iterations to run pyNBS
niter=$3
echo  'Number of pyNBS iterations: '$niter

# This example script assumes you are in the Command_Line_Tools directory of the pyNBS GitHub repository
# Please edit the following paths to point to the desired files for your own usage
pyNBS_script=$PWD'/run_pyNBS.py'
echo 'Python Script: '$pyNBS_script
network_file1='/cellar/users/jkhuang/Data/Projects/pyNBS/Data/CancerSubnetwork_v1_old.txt'
echo 'Network File: '$network_file
mutation_file=$PWD'/Example_Data/Mutation_Files/'$cancer_type'_sm_data.txt'
echo 'Somatic Mutation Data File: '$mutation_file
survival_data=$PWD'/Example_Data/Clinical_Files/'$cancer_type'.clin.merged.surv.txt'
echo 'Survival Data File: '$survival_data
params_file=$PWD'/'$cancer_type'_run_pyNBS_params.csv'
echo 'pyNBS Parameters File: '$params_file

# Setting optional parameters for pyNBS run
outdir=$PWD'/Results/via_script/'$cancer_type'/'
echo 'Output directory: '$outdir
job_name1=$cancer_type'_CancerSubnetwork_old_k'$k
echo 'Job Name: '$job_name1

python $pyNBS_script $mutation_file $network_file1 -params $params_file -o $outdir -j $job_name1 -k $k -n $niter -surv $survival_data

# # Run pyNBS again with "new" Cancer Subnetwork
# network_file2=$PWD'/Example_Data/Network_Files/CancerSubnetwork.txt'
# echo 'Network File: '$network_file
# job_name2=$cancer_type'_CancerSubnetwork_new_k'$k
# echo 'Job Name: '$job_name2
# python $pyNBS_script $mutation_file $network_file2 -params $params_file -o $outdir -j $job_name2 -k $k -n $niter -surv $survival_data