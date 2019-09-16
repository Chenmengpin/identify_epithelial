run_infercnv <- function(infercnv_object, ncores, outdir, cutoff, window, 
  threshold, denoise_value, typeHMM, run_mode="sample") {
  if (run_mode == "subpop") {
    infercnv::run(infercnv_object,
      num_threads=ncores,
      out_dir=outdir,
      cutoff=cutoff,
      window_length=window,
      max_centered_threshold=threshold,
      cluster_by_groups=T,
      plot_steps=F,
      denoise=T,
      sd_amplifier=denoise_value,
      HMM=T,
      HMM_type=typeHMM,
      analysis_mode="subclusters",
      HMM_report_by="subcluster"
    )
  } else if (run_mode == "sample") {
    infercnv::run(infercnv_object,
      num_threads=ncores,
      out_dir=outdir,
      cutoff=cutoff,
      window_length=window,
      max_centered_threshold=threshold,
      cluster_by_groups=T,
      plot_steps=F,
      denoise=T,
      sd_amplifier=denoise_value,
      HMM=T,
      HMM_type=typeHMM,
      analysis_mode="samples",
      HMM_report_by="consensus"
    )
  } else if (run_mode == "cell") {
    infercnv::run(infercnv_object,
      num_threads=ncores,
      out_dir=outdir,
      cutoff=cutoff,
      window_length=window,
      max_centered_threshold=threshold,
      cluster_by_groups=T,
      plot_steps=F,
      denoise=T,
      sd_amplifier=denoise_value,
      HMM=T,
      HMM_type=typeHMM,
      analysis_mode="cells",
      HMM_report_by="cell"
    )
  }
}
