#! /bin/bash

# Which shell to use
#$ -S /bin/bash
# Which queue to use
#$ -l long
# Transfer all variables to job script (e.g. PATH, LD_LIBRARY_PATH, etc.)
#$ -V

date
# Get task parameters
param_file=$1
task_params=($(sed "${SGE_TASK_ID}q;d" $param_file))
echo '-------------'
mutation_file=${task_params[0]}
echo 'Somatic Mutation Data File: '$mutation_file
network_file=${task_params[1]}
echo 'Network File: '$network_file
params_file=${task_params[2]}
echo 'pyNBS Parameters File: '$params_file
outdir=${task_params[3]}
echo 'Output Directory: '$outdir
job_name=${task_params[4]}
echo 'Job Name: '$job_name
a=${task_params[5]}
echo 'Alpha: '$a
k=${task_params[6]}
echo 'Nubmer of Clusters: '$k
n=${task_params[7]}
echo 'pyNBS Iterations: '$n
survival_file=${task_params[8]}
echo 'Survival Data File: '$survival_file
threads=${task_params[9]}
echo 'Number of threads: '$threads

#Run python script
pyNBS_script='/cellar/users/jkhuang/Data/Projects/pyNBS/pyNBS/Examples/run_pyNBS.py'
echo 'Python Script: '$pyNBS_script
echo 'pyNBS Call:'
echo python $pyNBS_script $mutation_file $network_file -params $params_file -o $outdir -j $job_name -a $a -k $k -n $n -surv $survival_file -t $threads
echo '-------------'
python $pyNBS_script $mutation_file $network_file -params $params_file -o $outdir -j $job_name -a $a -k $k -n $n -surv $survival_file -t $threads
date
