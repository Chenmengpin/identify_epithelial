#!/bin/bash

#  commands to monitor jobs:
#while :; do sleep 5; tail seurat_CID4290A/Output/Logs/CID4290A.infercnv/2.infercnv.Rout; done
#while :; do sleep 5; qstat -j 2194509 | grep usage; done

#sample_ids=( CID4463 CID3941 CID4067 CID4461 CID4465 CID4398 CID3921 CID3963 
#	CID45171 CID44041 CID4523 CID4290A CID4471 CID4535 CID3586 CID4066 CID4495 CID44971 
#	CID44991 CID4513 CID4530 CID4515 )
sample_ids=( CID3948 )

project_name="identify_epithelial"
subproject_name="brca_mini_atlas_130819"
ncores=5
subset_data="FALSE"
HMM="FALSE"
run_mode="subpop"
combined_sample="FALSE"
split_dataset="FALSE"
throw_unknowns="TRUE"
pool_cells="TRUE"
draw_metabric_annotations="FALSE"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/single_cell/$project_name"
seurat_path="$project_dir/results/seurat/$subproject_name"
script_dir="$project_dir/scripts/pipeline_v3"

for sample in ${sample_ids[@]}; do
  id=$(echo $sample | sed "s/seurat_//")
  short_id=$(echo $id | sed "s/CID//")
  echo "Calling InferCNV script for $id..."
  log_dir="$seurat_path/seurat_$id/Output/Logs/$id.infercnv/"
  mkdir -p $log_dir
  echo "Logs are in $log_dir"
  echo qsub -wd $log_dir -pe smp $ncores -N c10.$short_id -b y -j y -V -P TumourProgression \
    "${R} CMD BATCH  --no-save '--args $ncores $id $subset_data $HMM $run_mode $combined_sample $subproject_name $split_dataset $throw_unknowns $pool_cells $draw_metabric_annotations' $script_dir/2.infercnv.R"
  #/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/2.infercnv.R $ncores $id $subset_data $HMM $run_mode $combined_sample $subproject_name $split_dataset $throw_unknowns $pool_cells $draw_metabric_annotations
done;

