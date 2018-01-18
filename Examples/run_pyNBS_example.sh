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
network_file=$PWD'/Example_Data/Network_Files/CancerSubnetwork.txt'
mutation_file=$PWD'/Example_Data/Mutation_Files/'$cancer_type'_sm_data.txt'
survival_data=$PWD'/Example_Data/Clinical_Files/'$cancer_type'.clin.merged.surv.txt'

# Setting optional parameters for pyNBS run
outdir='./Results/via_command_line_tool/'$cancer_type'/'
job_name=$cancer_type'_pyNBS'

python $pyNBS_script $mutation_file $network_file -o $outdir -j $job_name -k $k -n $niter -surv $survival_data

