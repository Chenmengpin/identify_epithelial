run_infercnv <- function(infercnv_object, ncores, outdir, cutoff, window, 
  threshold, denoise_value, runHMM, typeHMM) {
  if (runHMM) {
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
            HMM=runHMM,
            HMM_type=typeHMM,
            analysis_mode="subclusters"
      )
  } else {
    infercnv::run(infercnv_object,
            num_threads=ncores,
            out_dir=outdir,
            cutoff=cutoff,
            window_length=window,
            max_centered_threshold=threshold,
            cluster_by_groups=T,
            plot_steps=F,
            denoise=T,
            sd_amplifier=denoise_value
    )
  }
}
