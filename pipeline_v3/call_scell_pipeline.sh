# run from /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/results/seurat/brca_mini_atlas
R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R"
sample_names=( $(echo $1 | sed "s/\_/\ /g")  )
echo ${sample_names[@]}
#script_name="seurat_v3_individual_processing/seurat_V3_HPC_with_garnett_classifier_and_infercnv_working_copy.R"
script_name="seurat_V3_HPC_with_garnett_classifier_and_infercnv.R"
for s in ${sample_names[@]}; do 
	sname="seurat_$s"
	echo $sname
	cd $sname
	"${R} CMD BATCH --no-save '--args $sname human $(pwd) /share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4404/count_4404_primary_GRCh38/outs/raw_gene_bc_matrices/GRCh38/ /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/seurat_v3_individual_processing/seurat_gene_input_file.csv /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/seurat_v3_individual_processing/seurat_V3_HPC_params_file.csv /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/refs/' /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/normal_threshold/$script_name"
	cd ..
done