# testing parameters:
infercnv_obj <- initial_infercnv_object
cutoff <- 0.1
min_cells_per_gene = 3
window_length = 101
smooth_method = "pyramidinal"
num_ref_groups = NULL
ref_subtract_use_mean_bounds = TRUE
cluster_by_groups = FALSE
k_obs_groups = 1
hclust_method = "ward.D2"
max_centered_threshold = 3
scale_data = FALSE
HMM = FALSE
HMM_transition_prob = 1e-06
HMM_report_by = "subcluster"
HMM_type = "i6"
HMM_i3_pval = 0.05
HMM_i3_use_KS = TRUE
BayesMaxPNormal = 0.5
sim_method = "meanvar"
sim_foreground = FALSE
analysis_mode = "subclusters"
tumor_subcluster_partition_method = "random_trees"
tumor_subcluster_pval = 0.1
denoise = FALSE
noise_filter = NA
sd_amplifier = 1.5
noise_logistic = FALSE
outlier_method_bound = "average_bound"
outlier_lower_bound = NA
outlier_upper_bound = NA
final_scale_limits = NULL
final_center_val = NULL
debug = FALSE
num_threads = 10
plot_steps = FALSE
resume_mode = TRUE
png_res = 300
plot_probabilities = TRUE
diagnostics = FALSE
remove_genes_at_chr_ends = FALSE
prune_outliers = FALSE
mask_nonDE_genes = FALSE
mask_nonDE_pval = 0.05
test.use = "wilcoxon"
require_DE_all_normals = "any"
hspike_aggregate_normals = FALSE
no_plot = FALSE
no_prelim_plot = FALSE


#infercnv:::run
test_run <- function (infercnv_obj, cutoff = 1, min_cells_per_gene = 3, out_dir = ".",
    window_length = 101, smooth_method = c("pyramidinal", "runmeans"),
    num_ref_groups = NULL, ref_subtract_use_mean_bounds = TRUE,
    cluster_by_groups = FALSE, k_obs_groups = 1, hclust_method = "ward.D2",
    max_centered_threshold = 3, scale_data = FALSE, HMM = FALSE,
    HMM_transition_prob = 1e-06, HMM_report_by = c("subcluster",
        "consensus", "cell"), HMM_type = c("i6", "i3"), HMM_i3_pval = 0.05,
    HMM_i3_use_KS = TRUE, BayesMaxPNormal = 0.5, sim_method = "meanvar",
    sim_foreground = FALSE, analysis_mode = c("samples", "subclusters",
        "cells"), tumor_subcluster_partition_method = c("random_trees",
        "qnorm", "pheight", "qgamma", "shc"), tumor_subcluster_pval = 0.1,
    denoise = FALSE, noise_filter = NA, sd_amplifier = 1.5, noise_logistic = FALSE,
    outlier_method_bound = "average_bound", outlier_lower_bound = NA,
    outlier_upper_bound = NA, final_scale_limits = NULL, final_center_val = NULL,
    debug = FALSE, num_threads = 4, plot_steps = FALSE, resume_mode = TRUE,
    png_res = 300, plot_probabilities = TRUE, diagnostics = FALSE,
    remove_genes_at_chr_ends = FALSE, prune_outliers = FALSE,
    mask_nonDE_genes = FALSE, mask_nonDE_pval = 0.05, test.use = "wilcoxon",
    require_DE_all_normals = "any", hspike_aggregate_normals = FALSE,
    no_plot = FALSE, no_prelim_plot = FALSE)
{
    smooth_method <<- match.arg(smooth_method)
    HMM_report_by = match.arg(HMM_report_by)
    analysis_mode = match.arg(analysis_mode)
    tumor_subcluster_partition_method = match.arg(tumor_subcluster_partition_method)
    HMM_type = match.arg(HMM_type)
    if (debug) {
        flog.threshold(DEBUG)
    }
    else {
        flog.threshold(INFO)
    }
    flog.info(paste("::process_data:Start", sep = ""))
    infercnv.env$GLOBAL_NUM_THREADS <- num_threads
    if (out_dir != "." & !file.exists(out_dir)) {
        dir.create(out_dir)
    }
    step_count = 0
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %d: incoming data\n", step_count))
    infercnv_obj_file = file.path(out_dir, sprintf("%02d_incoming_data.infercnv_obj",
        step_count))
    if (!(resume_mode & file.exists(infercnv_obj_file))) {
        saveRDS(infercnv_obj, infercnv_obj_file)
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Removing lowly expressed genes\n",
        step_count))
    infercnv_obj_file = file.path(out_dir, sprintf("%02d_reduced_by_cutoff.infercnv_obj",
        step_count))
    if (resume_mode & file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s",
            infercnv_obj_file))
        infercnv_obj = readRDS(infercnv_obj_file)
    }
    else {
        infercnv_obj <- require_above_min_mean_expr_cutoff(infercnv_obj,
            cutoff)
        infercnv_obj <- require_above_min_cells_ref(infercnv_obj,
            min_cells_per_gene = min_cells_per_gene)
        saveRDS(infercnv_obj, file = infercnv_obj_file)
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: normalization by sequencing depth\n",
        step_count))
    resume_file_token = ifelse((HMM), paste0("HMM", HMM_type),
        "")
    infercnv_obj_file = file = file.path(out_dir, sprintf("%02d_normalized_by_depth%s.infercnv_obj",
        step_count, resume_file_token))
    if (resume_mode & file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s",
            infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    }
    else {
        infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj)
        if (HMM && HMM_type == "i6") {
            infercnv_obj <- .build_and_add_hspike(infercnv_obj,
                sim_method = sim_method, aggregate_normals = hspike_aggregate_normals)
            if (sim_foreground) {
                infercnv_obj <- .sim_foreground(infercnv_obj,
                  sim_method = sim_method)
            }
        }
        saveRDS(infercnv_obj, infercnv_obj_file)
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: log transformation of data\n",
        step_count))
    infercnv_obj_file = file.path(out_dir, sprintf("%02d_logtransformed%s.infercnv_obj",
        step_count, resume_file_token))
    if (resume_mode & file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s",
            infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    }
    else {
        infercnv_obj <- log2xplus1(infercnv_obj)
        saveRDS(infercnv_obj, file = infercnv_obj_file)
        if (plot_steps) {
            plot_cnv(infercnv_obj = infercnv_obj, k_obs_groups = k_obs_groups,
                cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                title = sprintf("%02d_log_transformed_data",
                  step_count), output_filename = sprintf("infercnv.%02d_log_transformed",
                  step_count), write_expr_matrix = TRUE, png_res = png_res)
        }
    }
    if (scale_data) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: scaling all expression data\n",
            step_count))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_scaled%s.infercnv_obj",
            step_count, resume_file_token))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            infercnv_obj <- scale_infercnv_expr(infercnv_obj)
            saveRDS(infercnv_obj, file = infercnv_obj_file)
            if (plot_steps) {
                plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_scaled", step_count),
                  output_filename = sprintf("infercnv.%02d_scaled",
                    step_count), write_expr_matrix = TRUE, png_res = png_res)
            }
        }
    }
    if (!is.null(num_ref_groups)) {
        if (!has_reference_cells(infercnv_obj)) {
            stop("Error, no reference cells defined. Cannot split them into groups as requested")
        }
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: splitting reference data into %d clusters\n",
            step_count, num_ref_groups))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_split_%02d_refs%s.infercnv_obj",
            step_count, resume_file_token, num_ref_groups))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            infercnv_obj <- split_references(infercnv_obj, num_groups = num_ref_groups,
                hclust_method = hclust_method)
            saveRDS(infercnv_obj, file = infercnv_obj_file)
        }
    }
    if (analysis_mode == "subclusters" & tumor_subcluster_partition_method ==
        "random_trees") {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n",
            step_count, tumor_subcluster_partition_method))
        resume_file_token = paste0(resume_file_token, ".rand_trees")
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_tumor_subclusters%s.%s.infercnv_obj",
            step_count, resume_file_token, tumor_subcluster_partition_method))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            infercnv_obj <- define_signif_tumor_subclusters_via_random_smooothed_trees(infercnv_obj,
                p_val = tumor_subcluster_pval, hclust_method = hclust_method)
            saveRDS(infercnv_obj, file = infercnv_obj_file)
            if (plot_steps) {
                plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_tumor_subclusters.%s",
                    step_count, tumor_subcluster_partition_method),
                  output_filename = sprintf("infercnv.%02d_tumor_subclusters.%s",
                    step_count, tumor_subcluster_partition_method),
                  write_expr_matrix = TRUE, png_res = png_res)
            }
        }
    }
    else if (analysis_mode != "subclusters") {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Clustering samples (not defining tumor subclusters)\n",
            step_count))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_no_subclustering%s.infercnv_obj",
            step_count, resume_file_token))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            infercnv_obj <- define_signif_tumor_subclusters(infercnv_obj,
                p_val = tumor_subcluster_pval, hclust_method = hclust_method,
                partition_method = "none")
            saveRDS(infercnv_obj, file = infercnv_obj_file)
        }
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (before smoothing)\n",
        step_count))
    infercnv_obj_file = file.path(out_dir, sprintf("%02d_remove_ref_avg_from_obs_logFC%s.infercnv_obj",
        step_count, resume_file_token))
    if (resume_mode & file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s",
            infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    }
    else {
        infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj,
            inv_log = FALSE, use_bounds = ref_subtract_use_mean_bounds)
        saveRDS(infercnv_obj, file = infercnv_obj_file)
        if (plot_steps) {
            plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                title = sprintf("%02d_remove_average", step_count),
                output_filename = sprintf("infercnv.%02d_remove_average",
                  step_count), write_expr_matrix = TRUE, png_res = png_res)
        }
    }
    if (!is.na(max_centered_threshold)) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: apply max centered expression threshold: %s\n",
            step_count, max_centered_threshold))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_apply_max_centered_expr_threshold%s.infercnv_obj",
            step_count, resume_file_token))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            threshold = max_centered_threshold
            if (is.character(max_centered_threshold) && max_centered_threshold ==
                "auto") {
                threshold = mean(abs(get_average_bounds(infercnv_obj)))
                flog.info(sprintf("Setting max centered threshoolds via auto to: +- %g",
                  threshold))
            }
            infercnv_obj <- apply_max_threshold_bounds(infercnv_obj,
                threshold = threshold)
            saveRDS(infercnv_obj, file = infercnv_obj_file)
            if (plot_steps) {
                plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_apply_max_centered_expr_threshold",
                    step_count), output_filename = sprintf("infercnv.%02d_apply_max_centred_expr_threshold",
                    step_count), write_expr_matrix = TRUE, png_res = png_res)
            }
        }
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Smoothing data per cell by chromosome\n",
        step_count))
    infercnv_obj_file = file.path(out_dir, sprintf("%02d_smoothed_by_chr%s.infercnv_obj",
        step_count, resume_file_token))
    if (resume_mode & file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s",
            infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    }
    else {
        if (smooth_method == "runmeans") {
            infercnv_obj <- smooth_by_chromosome_runmeans(infercnv_obj,
                window_length)
        }
        else if (smooth_method == "pyramidinal") {
            infercnv_obj <- smooth_by_chromosome(infercnv_obj,
                window_length = window_length, smooth_ends = TRUE)
        }
        else {
            stop(sprintf("Error, don't recognize smoothing method: %s",
                smooth_method))
        }
        saveRDS(infercnv_obj, file = infercnv_obj_file)
        if (plot_steps) {
            plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                title = sprintf("%02d_smoothed_by_chr", step_count),
                output_filename = sprintf("infercnv.%02d_smoothed_by_chr",
                  step_count), write_expr_matrix = TRUE, png_res = png_res)
        }
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: re-centering data across chromosome after smoothing\n",
        step_count))
    infercnv_obj_file = file.path(out_dir, sprintf("%02d_recentered_cells_by_chr%s.infercnv_obj",
        step_count, resume_file_token))
    if (resume_mode & file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s",
            infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    }
    else {
        infercnv_obj <- center_cell_expr_across_chromosome(infercnv_obj,
            method = "median")
        saveRDS(infercnv_obj, file = infercnv_obj_file)
        if (plot_steps) {
            plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                title = sprintf("%02d_centering_of_smoothed",
                  step_count), output_filename = sprintf("infercnv.%02d_centering_of_smoothed",
                  step_count), write_expr_matrix = TRUE, png_res = png_res)
        }
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (after smoothing)\n",
        step_count))
    infercnv_obj_file = file.path(out_dir, sprintf("%02d_remove_ref_avg_from_obs_adjust%s.infercnv_obj",
        step_count, resume_file_token))
    if (resume_mode & file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s",
            infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    }
    else {
        infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj,
            inv_log = FALSE, use_bounds = ref_subtract_use_mean_bounds)
        saveRDS(infercnv_obj, file = infercnv_obj_file)
        if (plot_steps) {
            plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                title = sprintf("%02d_remove_average", step_count),
                output_filename = sprintf("infercnv.%02d_remove_average",
                  step_count), write_expr_matrix = TRUE, png_res = png_res)
        }
    }
    if (remove_genes_at_chr_ends == TRUE) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: removing genes at chr ends\n",
            step_count))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_remove_gene_at_chr_ends%s.infercnv_obj",
            step_count, resume_file_token))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            infercnv_obj <- remove_genes_at_ends_of_chromosomes(infercnv_obj,
                window_length)
            saveRDS(infercnv_obj, file = infercnv_obj_file)
            if (plot_steps) {
                plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_remove_genes_at_chr_ends",
                    step_count), output_filename = sprintf("infercnv.%02d_remove_genes_at_chr_ends",
                    step_count), write_expr_matrix = TRUE, png_res = png_res)
            }
        }
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: invert log2(FC) to FC\n",
        step_count))
    infercnv_obj_file = file.path(out_dir, sprintf("%02d_invert_log_transform%s.infercnv_obj",
        step_count, resume_file_token))
    if (resume_mode & file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s",
            infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    }
    else {
        infercnv_obj <- invert_log2(infercnv_obj)
        saveRDS(infercnv_obj, file = infercnv_obj_file)
        if (plot_steps) {
            plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                title = sprintf("%02d_invert_log_transform log(FC)->FC",
                  step_count), output_filename = sprintf("infercnv.%02d_invert_log_FC",
                  step_count), write_expr_matrix = TRUE, png_res = png_res)
        }
    }
    if (analysis_mode == "subclusters" & tumor_subcluster_partition_method !=
        "random_trees") {
        resume_file_token = paste0(resume_file_token, ".", tumor_subcluster_partition_method)
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n",
            step_count, tumor_subcluster_partition_method))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_tumor_subclusters%s.infercnv_obj",
            step_count, resume_file_token))
        if (resume_mode && file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            infercnv_obj <- define_signif_tumor_subclusters(infercnv_obj,
                p_val = tumor_subcluster_pval, hclust_method = hclust_method,
                partition_method = tumor_subcluster_partition_method)
            saveRDS(infercnv_obj, file = infercnv_obj_file)
            if (plot_steps) {
                plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_tumor_subclusters", step_count),
                  output_filename = sprintf("infercnv.%02d_tumor_subclusters",
                    step_count), write_expr_matrix = TRUE, png_res = png_res)
            }
        }
    }
    infercnv_obj_prelim <- infercnv_obj
    infercnv_obj_file = file.path(out_dir, "preliminary.infercnv_obj")
    saveRDS(infercnv_obj_prelim, file = infercnv_obj_file)
    if (!(no_prelim_plot | no_plot)) {
        prelim_heatmap_png = "infercnv.preliminary.png"
        if (!file.exists(file.path(out_dir, prelim_heatmap_png))) {
            plot_cnv(infercnv_obj_prelim, k_obs_groups = k_obs_groups,
                cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                title = "Preliminary infercnv (pre-noise filtering)",
                output_filename = "infercnv.preliminary", write_expr_matrix = TRUE,
                png_res = png_res)
        }
    }
    if (prune_outliers) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Removing outliers\n",
            step_count))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_remove_outlier%s.infercnv_obj",
            step_count, resume_file_token))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            infercnv_obj = remove_outliers_norm(infercnv_obj,
                out_method = outlier_method_bound, lower_bound = outlier_lower_bound,
                upper_bound = outlier_upper_bound)
            saveRDS(infercnv_obj, file = infercnv_obj_file)
            if (plot_steps) {
                plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_removed_outliers", step_count),
                  output_filename = sprintf("infercnv.%02d_removed_outliers",
                    step_count), write_expr_matrix = TRUE, png_res = png_res)
            }
        }
    }
    if (HMM) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: HMM-based CNV prediction\n",
            step_count))
        hmm_resume_file_token = paste0(resume_file_token, ".hmm_mode-",
            analysis_mode)
        hmm.infercnv_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred%s.infercnv_obj",
            step_count, hmm_resume_file_token))
        if (resume_mode & file.exists(hmm.infercnv_obj_file)) {
            flog.info(sprintf("-restoring hmm.infercnv_obj from %s",
                hmm.infercnv_obj_file))
            hmm.infercnv_obj <- readRDS(hmm.infercnv_obj_file)
        }
        else {
            if (HMM_type == "i6") {
                hmm_center = 3
                hmm_state_range = c(0, 6)
            }
            else {
                hmm_center = 2
                hmm_state_range = c(1, 3)
            }
            if (analysis_mode == "subclusters") {
                if (HMM_type == "i6") {
                  hmm.infercnv_obj <- predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                    t = HMM_transition_prob)
                }
                else if (HMM_type == "i3") {
                  hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                    i3_p_val = HMM_i3_pval, t = HMM_transition_prob,
                    use_KS = HMM_i3_use_KS)
                }
                else {
                  stop("Error, not recognizing HMM_type")
                }
            }
            else if (analysis_mode == "cells") {
                if (HMM_type == "i6") {
                  hmm.infercnv_obj <- predict_CNV_via_HMM_on_indiv_cells(infercnv_obj,
                    t = HMM_transition_prob)
                }
                else if (HMM_type == "i3") {
                  hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_indiv_cells(infercnv_obj,
                    i3_p_val = HMM_i3_pval, t = HMM_transition_prob,
                    use_KS = HMM_i3_use_KS)
                }
                else {
                  stop("Error, not recognizing HMM_type")
                }
            }
            else {
                if (HMM_type == "i6") {
                  hmm.infercnv_obj <- predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj,
                    t = HMM_transition_prob)
                }
                else if (HMM_type == "i3") {
                  hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                    i3_p_val = HMM_i3_pval, t = HMM_transition_prob,
                    use_KS = HMM_i3_use_KS)
                }
                else {
                  stop("Error, not recognizing HMM_type")
                }
            }
            saveRDS(hmm.infercnv_obj, file = hmm.infercnv_obj_file)
            generate_cnv_region_reports(hmm.infercnv_obj, output_filename_prefix = sprintf("%02d_HMM_preds",
                step_count), out_dir = out_dir, ignore_neutral_state = hmm_center,
                by = HMM_report_by)
            if (!no_plot) {
                plot_cnv(infercnv_obj = hmm.infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_HMM_preds", step_count),
                  output_filename = sprintf("infercnv.%02d_HMM_pred%s",
                    step_count, hmm_resume_file_token), write_expr_matrix = TRUE,
                  x.center = hmm_center, x.range = hmm_state_range,
                  png_res = png_res)
            }
        }
        if (HMM_type == "i6" & BayesMaxPNormal > 0) {
            step_count = step_count + 1
            flog.info(sprintf("\n\n\tSTEP %02d: Run Bayesian Network Model on HMM predicted CNV's\n",
                step_count))
            mcmc_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj",
                step_count, hmm_resume_file_token))
            if (resume_mode & file.exists(mcmc_obj_file)) {
                flog.info(sprintf("-restoring mcmc_obj from %s",
                  mcmc_obj_file))
                mcmc_obj <- readRDS(mcmc_obj_file)
            }
            else {
                mcmc_obj <- infercnv::inferCNVBayesNet(infercnv_obj = infercnv_obj_prelim,
                  HMM_states = hmm.infercnv_obj@expr.data, file_dir = out_dir,
                  postMcmcMethod = "removeCNV", out_dir = file.path(out_dir,
                    sprintf("BayesNetOutput.%s", hmm_resume_file_token)),
                  quietly = TRUE, CORES = num_threads, plotingProbs = plot_probabilities,
                  diagnostics = diagnostics)
                saveRDS(mcmc_obj, file = mcmc_obj_file)
            }
            mcmc_obj <- infercnv::filterHighPNormals(MCMC_inferCNV_obj = mcmc_obj,
                HMM_states = hmm.infercnv_obj@expr.data, BayesMaxPNormal = BayesMaxPNormal)
            hmm.infercnv_obj <- infercnv::returningInferCNV(mcmc_obj,
                hmm.infercnv_obj@expr.data)
            mcmc.infercnv_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj",
                step_count, hmm_resume_file_token, BayesMaxPNormal))
            saveRDS(hmm.infercnv_obj, file = mcmc.infercnv_obj_file)
            if (!no_plot) {
                plot_cnv(infercnv_obj = hmm.infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_HMM_preds_Bayes_Net",
                    step_count), output_filename = sprintf("infercnv.%02d_HMM_pred.Bayes_Net.Pnorm_%g",
                    step_count, BayesMaxPNormal), write_expr_matrix = TRUE,
                  x.center = 3, x.range = c(0, 6), png_res = png_res)
            }
        }
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Converting HMM-based CNV states to repr expr vals\n",
            step_count))
        hmm.infercnv_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.repr_intensities%s.Pnorm_%g.infercnv_obj",
            step_count, hmm_resume_file_token, BayesMaxPNormal))
        if (resume_mode & file.exists(hmm.infercnv_obj_file)) {
            flog.info(sprintf("-restoring hmm.infercnv_obj from %s",
                hmm.infercnv_obj_file))
            hmm.infercnv_obj <- readRDS(hmm.infercnv_obj_file)
        }
        else {
            if (HMM_type == "i6") {
                hmm.infercnv_obj <- assign_HMM_states_to_proxy_expr_vals(hmm.infercnv_obj)
            }
            else if (HMM_type == "i3") {
                hmm.infercnv_obj <- i3HMM_assign_HMM_states_to_proxy_expr_vals(hmm.infercnv_obj)
            }
            saveRDS(hmm.infercnv_obj, file = hmm.infercnv_obj_file)
            if (!no_plot) {
                plot_cnv(infercnv_obj = hmm.infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_HMM_preds.repr_intensities",
                    step_count), output_filename = sprintf("infercnv.%02d_HMM_pred%s.repr_intensities",
                    step_count, hmm_resume_file_token), write_expr_matrix = TRUE,
                  x.center = 1, x.range = c(-1, 3), png_res = png_res)
            }
            generate_cnv_region_reports(hmm.infercnv_obj, output_filename_prefix = sprintf("HMM_CNV_predictions.%s.Pnorm_%g",
                hmm_resume_file_token, BayesMaxPNormal), out_dir = out_dir,
                ignore_neutral_state = 1, by = HMM_report_by)
        }
    }
    if (mask_nonDE_genes) {
        if (!has_reference_cells(infercnv_obj)) {
            stop("Error, cannot mask non-DE genes when there are no normal references set")
        }
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Identify and mask non-DE genes\n",
            step_count))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_mask_nonDE%s.infercnv_obj",
            step_count, resume_file_token))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            infercnv_obj <- mask_non_DE_genes_basic(infercnv_obj,
                p_val_thresh = mask_nonDE_pval, test.use = test.use,
                center_val = mean(infercnv_obj@expr.data), require_DE_all_normals = require_DE_all_normals)
            saveRDS(infercnv_obj, file = infercnv_obj_file)
            if (plot_steps) {
                plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  title = sprintf("%02d_mask_nonDE", step_count),
                  output_filename = sprintf("infercnv.%02d_mask_nonDE",
                    step_count), write_expr_matrix = TRUE, png_res = png_res)
            }
        }
    }
    if (denoise) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Denoising\n", step_count))
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_denoise%s.NF_%s.SD_%g.NL_%s.infercnv_obj",
            step_count, resume_file_token, noise_filter, sd_amplifier,
            noise_logistic))
        if (resume_mode & file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s",
                infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        }
        else {
            if (!is.na(noise_filter)) {
                if (noise_filter > 0) {
                  flog.info(paste("::process_data:Remove noise, noise threshold at: ",
                    noise_filter))
                  infercnv_obj <- clear_noise(infercnv_obj, threshold = noise_filter,
                    noise_logistic = noise_logistic)
                }
                else {
                }
            }
            else {
                flog.info(paste("::process_data:Remove noise, noise threshold defined via ref mean sd_amplifier: ",
                  sd_amplifier))
                infercnv_obj <- clear_noise_via_ref_mean_sd(infercnv_obj,
                  sd_amplifier = sd_amplifier, noise_logistic = noise_logistic)
            }
            saveRDS(infercnv_obj, file = infercnv_obj_file)
            if (!no_plot) {
                plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups,
                  cluster_by_groups = cluster_by_groups, out_dir = out_dir,
                  color_safe_pal = FALSE, title = sprintf("%02d_denoised",
                    step_count), output_filename = sprintf("infercnv.%02d_denoised",
                    step_count), write_expr_matrix = TRUE, png_res = png_res)
            }
        }
    }
    saveRDS(infercnv_obj, file = file.path(out_dir, "run.final.infercnv_obj"))
    if (!no_plot) {
        if (is.null(final_scale_limits)) {
            final_scale_limits = "auto"
        }
        if (is.null(final_center_val)) {
            final_center_val = 1
        }
        flog.info("\n\n## Making the final infercnv heatmap ##")
        plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups,
            out_dir = out_dir, x.center = final_center_val, x.range = final_scale_limits,
            title = "inferCNV", output_filename = "infercnv",
            write_expr_matrix = TRUE, png_res = png_res)
    }
    return(infercnv_obj)
}