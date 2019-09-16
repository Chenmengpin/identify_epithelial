run_infercnv <- function(infercnv_object, ncores, outdir, cutoff, window, 
  threshold, denoise_value) {
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
