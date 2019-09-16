#!/bin/bash

anchor_no="25000"
dims=( "75" "100" "125" "150" )
numcores=$(($anchor_no/800))
hvmem=$(($anchor_no/25))


home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/single_cell/identify_epithelial"
script_dir="$project_dir/scripts"
R_dir="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin"

for dim in ${dims[@]}; do
	qsub -N cca.$anchor_no.$dim -b y \
	-wd $script_dir -j y -R y -pe smp $numcores -l h_vmem="$hvmem"G -V \
	"$R_dir/Rscript \
	--vanilla $script_dir/combine_seurat_objects.R $anchor_no $dim"
done;