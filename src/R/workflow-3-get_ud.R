#' @title UD functions

get_ud_path <- function(sim, path,
                        spat, win, sigma = spatstat.explore::bw.diggle,
                        overwrite = TRUE) {
  out_file <- here_alg(sim, "path", "ud.tif")
  if (overwrite | !file.exists(out_file)) {
    ud_path <-
      map_dens(
        .map = spat,
        .owin = win,
        .coord = path[, .(x, y)],
        .discretise = TRUE,
        .plot = FALSE,
        .verbose = FALSE,
        sigma = sigma
      )$ud
    write_rast(ud_path, out_file)
  } else {
    ud_path <- terra::rast(out_file)
  }
  ud_path
}

get_ud_coa <- function(sim,
                       detections, delta_t,
                       spat, win, sigma = spatstat.explore::bw.diggle) {
  out_file  <- here_alg(sim, "coa", delta_t, "ud.tif")
  if (secs(max(detections$timestamp), min(detections$timestamp)) <=
      as.numeric(lubridate::duration(delta_t))) {
    return(FALSE)
  }
  out_coa <- coa(.map = spat,
                 .acoustics = detections[obs == 1L, ],
                 .delta_t = delta_t,
                 .plot_weights = FALSE)
  # Use .discretise = TRUE for consistency
  ud_coa  <-
    map_dens(.map = spat,
             .owin = win,
             .coord = out_coa[, .(x, y)],
             .discretise = TRUE,
             .plot = FALSE,
             .verbose = FALSE,
             sigma = sigma)$ud
  write_rast(ud_coa, out_file)
  TRUE
}

get_ud_rsp <- function(sim, spat, spat_ll_dbb, tm, type = c("default", "custom"),
                       overwrite = TRUE) {
  # Define outfile
  type <- match.arg(type)
  out_file <- here_alg(sim, "rsp", type, "ud.tif")
  if (!overwrite && file.exists(out_file)) {
    return("exists")
  }
  # Read actel
  act <- read_actel(sim)
  # Examine locations
  if (FALSE) {
    water <- spat_ll_dbb > 0
    names(water) <- "layer"
    RSP::plotRaster(input = act,
                    base.raster = water,
                    coord.x = "Longitude", coord.y = "Latitude")
  }
  # Define RSP args
  # * Arbitrarily set max.time (days) to any time greater than study duration
  args <- list(input = act,
               t.layer = tm,
               coord.x = "Longitude", coord.y = "Latitude",
               max.time = 100)
  if (type == "custom") {
    # dynBBMM() is sensitive to raster area
    # * Larger er.ad appears to require larger areas
    # * Hence small modification of er.ad argument only
    # args$distance   <- 100
    # args$time.step  <- 0.4
    # args$min.time   <- 2
    args$er.ad        <- 250 * 0.10
  }
  # Run RSP
  t1_rsp <- Sys.time()
  out_rsp <- tryCatch(
    do.call(RSP::runRSP, args),
    error = function(e) e)
  if (inherits(out_rsp, "error")) {
    warning("RSP::runRSP failure!")
    warning(out_rsp$message)
    return("runRSP")
  }
  t2_rsp <- Sys.time()
  # Generate UD
  t1_ud <- Sys.time()
  dbb <- tryCatch(
    RSP::dynBBMM(input = out_rsp,
                 base.raster = spat_ll_dbb,
                 UTM = "1"),
    error = function(e) e)
  if (inherits(dbb, "error")) {
    warning("RSP::dynBBMM() failure!")
    warning(dbb$message)
    return("dynBBMM")
  }
  dbb <- dbb$dbbmm[[1]]
  dbb <- terra::rast(dbb)
  if (terra::nlyr(dbb) > 1L) {
    warning(paste("Multi-layered UD: sim$id", sim$id))
  }
  # Resample RSP onto spat grid for consistency
  ud_rsp <- terra::project(dbb, terra::crs(spat))
  ud_rsp <- terra::resample(ud_rsp, spat)
  ud_rsp <- ud_rsp / terra::global(ud_rsp, "sum")[1, 1]
  # stopifnot(all.equal(1, terra::global(ud_rsp, "sum")[1, 1]))
  t2_ud <- Sys.time()

  # Save outputs
  time <- data.table(id = sim$id,
                     rsp = mins(t2_rsp, t1_rsp),
                     ud = mins(t2_ud, t1_ud))
  qs::qsave(time, here_alg(sim, "rsp", type, "time.qs"))
  write_rast(ud_rsp, out_file)
  "success"
}

get_ud_patter <- function(sim,
                          timeline, yobs,
                          algorithm = c("acpf", "acdcpf"),
                          spat, win, sigma = spatstat.explore::bw.diggle,
                          overwrite = TRUE) {

  #### (optional) Prior convergence check
  out_file_convergence <-
    here_alg(sim, "patter", algorithm, sim$alg_par, "convergence.rds")
  if (!overwrite && file.exists(out_file_convergence)) {
    return(readRDS(out_file_convergence))
  }

  #### Filter args
  # Observations
  yobs # TO DO
  # Models
  model_obs <- c("ModelObsAcousticLogisTrunc", "ModelObsDepthUniform")
  model_move <-
    move_xy(dbn_length =
              glue("truncated(Gamma({sim$shape},
                                    {sim$scale}),
                              upper = {sim$mobility})"),
            dbn_angle = "Uniform(-pi, pi)")
  # List
  args <- list(.map = spat,
               .timeline = timeline,
               .state = "StateXY",
               .yobs = yobs,
               .model_obs = model_obs,
               .model_move = model_move,
               .n_particle = sim$n_particles,
               .n_resample = sim$n_particles,
               .verbose = FALSE
               )

  #### Forward filter
  t1_pff   <- Sys.time()
  out_pff  <- do.call(pf_filter, args, quote = TRUE)
  t2_pff   <- Sys.time()
  pff_mins <- mins(t2_pff, t1_pff)
  if (!out_pff$convergence) {
    saveRDS(FALSE, out_file_convergence)
    return(FALSE)
  }

  #### Backward filter
  t1_pfb  <- Sys.time()
  args$direction <- "backward"
  out_pfb  <- do.call(pf_filter, args, quote = TRUE)
  t2_pfb   <- Sys.time()
  pfb_mins <- mins(t2_pfb, t1_pfb)
  if (!out_pfb$convergence) {
    saveRDS(FALSE, out_file_convergence)
    return(FALSE)
  }

  #### Backward smoother
  # Record successful convergence
  saveRDS(TRUE, out_file_convergence)
  # Implement sampler
  pfbs_mins   <- NA_real_
  run_sampler <- TRUE
  if (run_sampler) {
    t1_pfbs  <- Sys.time()
    # (optional) TO DO
    # * Improve speed by pre-defining box in Julia
    out_smo <- pf_smoother_two_filter(.map = spat,
                                      .mobility = sim$mobility,
                                      .n_particle = 1000L)
    t2_pfbs  <- Sys.time()
    pfbs_mins <- mins(t2_pfbs, t1_pfbs)
  }

  #### Map arguments
  # Use .discretise = TRUE for speed
  map_args <- list(.map = spat,
                   .owin = win,
                   .coord = out_pff$states,
                   .discretise = TRUE,
                   .plot = FALSE,
                   .verbose = FALSE,
                   sigma = sigma)

  #### Mapping (forward run)
  t1_udf          <- Sys.time()
  udf             <- do.call(map_dens, map_args)$ud
  t2_udf          <- Sys.time()
  udf_mins        <- mins(t2_udf, t1_udf)

  #### Mapping (backward run)
  # (For speed, this is not currently implemented)

  #### Mapping (backward sampler)
  uds_mins          <- NA_real_
  if (run_sampler) {
    map_args$.coord <- out_smo$states
    t1_uds          <- Sys.time()
    uds             <- do.call(map_dens, map_args)$ud
    t2_uds          <- Sys.time()
    uds_mins        <- mins(t2_uds, t1_uds)
  }

  #### Outputs
  # Timings
  time <- data.table(id = sim$id,
                     pff = pff_mins,
                     pfbs = pfbs_mins,
                     udf = udf_mins,
                     uds = uds_mins)

  qs::qsave(time, here_alg(sim, "patter", algorithm, sim$alg_par, "time.qs"))
  # Particle samples (for checks)
  # qs::qsave(out_pff$states, here_alg(sim, "patter", algorithm, sim$alg_par, "particles.qs"))
  # UD
  write_rast(udf, here_alg(sim, "patter", algorithm, sim$alg_par, "ud-f.tif"))
  if (run_sampler) {
    write_rast(uds, here_alg(sim, "patter", algorithm, sim$alg_par, "ud-s.tif"))
  }
  return(TRUE)
}
