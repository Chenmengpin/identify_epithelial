# run from /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/results/seurat/brca_mini_atlas
sample_names=( seurat_IND5 )
script_name="seurat_v3_individual_processing/seurat_V3_HPC_with_garnett_classifier_and_infercnv.R"
for s in ${sample_names[@]}; do 
	echo $s
	cd $s
	id=$(echo $s | sed "s/seurat_//")
	echo $id
	qsub -cwd -pe smp 32 -l h_vmem=300G -b y -j y -V -P TumourProgression "${R} CMD BATCH --no-save '--args $s human $(pwd) /share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4404/count_4404_primary_GRCh38/outs/raw_gene_bc_matrices/GRCh38/ /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/seurat_v3_individual_processing/seurat_gene_input_file.csv /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/seurat_v3_individual_processing/seurat_V3_HPC_params_file.csv /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/refs/' /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/$script_name"
	cd ..
done