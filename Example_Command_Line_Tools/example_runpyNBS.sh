#! /bin/bash

# cancer_type in the examples are 'BLCA' 'COAD' 'HNSC' 'OV' 'UCEC'
# To use this bash script, please modify all file path below to indicate your local directory
# For more information on pyNBS options, please refer to run_pyNBS.py for detailed documentation

cancer_type='BLCA'
echo $cancer_type
pyNBS_script='~/pyNBS/Example_Command_Line_Tools/run_pyNBS.py'
network_file='~/pyNBS/data/Network_Files/CancerSubnetwork.txt'
mutation_file='~/Data/Cancer_Mutations/'$cancer_type'_sm_data.txt'

network_file_basename=${network_file##*/}
network_name=${network_file_basename%.*}


clincal_file='~/pyNBS/data/Clinical_Files/'$cancer_type'.clin.merged.surv.txt'
result_dir='~/pyNBS/Example_Command_Line_Tools/results/'$cancer_type'/'
mkdir -p $result_dir

cc_mat_file=$result_dir$cancer_type'_'$network_name'_k_'$K'_cc_mat.csv'
cluster_assign_file=$result_dir$cancer_type'_'$network_name'_k_'$K'_cluster_assign.csv'
map_file=$result_dir$cancer_type'_'$network_name'_k_'$K'_map.png'
km_file=$result_dir$cancer_type'_'$network_name'_k_'$K'_km.png'

python $pyNBS_script $mutation_file $network_file -v -mf list -md '\t' -reg -n 100 -k $K -netNMF_n 250 -s_cc $cc_mat_file -s_ca $cluster_assign_file -s_map $map_file -km_plot True -clin $clincal_file -km_title $cancer_type' Survival Plot' -s_km $km_file

