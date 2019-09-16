#!/bin/bash

project_name="identify_epithelial"
subproject_name="brca_mini_atlas_030719"
sample_id="CID3963"
ncores=3
h_vmem="5M"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/single_cell/$project_name"
seurat_path="$project_dir/results/seurat/$subproject_name"
script_dir="$project_dir/scripts/"

echo "Calling InferCNV script for $sample_id with $ncores cores and $h_vmem memory..."
log_dir="$seurat_path/seurat_$sample_id/Output/Logs/$sample_id.subset.test.$ncores.cores.$h_vmem.memory/"
mkdir -p $log_dir
echo "Logs are in $log_dir"

qsub -wd $log_dir -pe smp $ncores -l h_vmem=$h_vmem -N $sample_id.subset.test.$ncores.cores.$h_vmem_memory -b y -j y -V -P TumourProgression \
    "${R} CMD BATCH  --no-save '--args $sample_id $ncores $h_vmem' $script_dir/subset_for_test_infercnv_memory.R"

