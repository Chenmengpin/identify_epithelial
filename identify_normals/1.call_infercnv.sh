#!/bin/bash

#  commands to monitor jobs:
#while :; do sleep 5; tail seurat_CID4290A/Output/Logs/CID4290A.infercnv/2.infercnv.Rout; done
#while :; do sleep 5; qstat -j 2194509 | grep usage; done

#sample_ids=(
#  CID3586 CID3921 CID3941 CID3948 CID3963 CID4066 CID4067 CID4290A CID4398 
#  CID44041 CID4461 CID4463 CID4465 CID4471 CID4495 CID44971 CID44991 CID4513 
#  CID4515 CID45171 CID4523 CID4530 CID4535
#)

R_dir="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin"

sample_ids=( CID4463 )

project_name="identify_epithelial"
subproject_name="identify_normals"
ncores=6
subset_data="FALSE"
draw_bulk_annotations="FALSE"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/single_cell/$project_name"
seurat_path="$project_dir/results/$subproject_name"
script_dir="$project_dir/scripts/$subproject_name"

for sample in ${sample_ids[@]}; do
  id=$(echo $sample | sed "s/seurat_//")
  short_id=$(echo $id | sed "s/CID//")
  echo "Calling InferCNV script for $id..."
  log_dir="$seurat_path/seurat_$id/Output/Logs/$id.infercnv/"
  mkdir -p $log_dir
  echo "Logs are in $log_dir"
  qsub -wd $log_dir -pe smp $ncores -N c$ncores.$short_id.infercnv -b y -j y -V -P TumourProgression \
  "$R_dir/R CMD BATCH  --no-save '--args $subproject_name $id $subset_data $ncores $draw_bulk_annotations' $script_dir/2.infercnv.R"
  #/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/2.infercnv.R $subproject_name $id $subset_data $ncores $draw_bulk_annotations
done;

