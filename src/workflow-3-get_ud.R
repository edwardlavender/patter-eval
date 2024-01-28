#' @title UD functions

get_ud_path <- function(sim,
                        spat, path,
                        im, win,
                        sigma = spatstat.explore::bw.diggle,
                        overwrite = TRUE) {
  out_file <- here_alg(sim, "path", "ud.tif")
  if (overwrite | !file.exists(out_file)) {
    # Use .discretise = TRUE for consistency
    # (The path is already defined on the grid)
    ud_path <-
      map_dens(.map = spat,
               .im = im,
               .win = win,
               .coord = path[, .(x, y)],
               .discretise = TRUE,
               .plot = FALSE,
               .verbose = FALSE,
               sigma = sigma)
    write_rast(ud_path, out_file)
  } else {
    ud_path <- terra::rast(out_file)
  }
  ud_path
}

get_ud_coa <- function(sim,
                       spat, dlist, delta_t,
                       im, win,
                       sigma = spatstat.explore::bw.diggle) {
  out_file  <- here_alg(sim, "coa", delta_t, "ud.tif")
  acoustics <- dlist$data$acoustics
  if (secs(max(acoustics$timestamp), min(acoustics$timestamp)) <=
      as.numeric(lubridate::duration(delta_t))) {
    return(FALSE)
  }
  out_coa <- coa(.dlist = dlist,
                 .delta_t = delta_t,
                 .plot_weights = FALSE)
  # Use .discretise = TRUE for consistency
  ud_coa  <-
    map_dens(.map = spat,
             .im = im,
             .win = win,
             .coord = out_coa[, .(x = coa_x, y = coa_y)],
             .discretise = TRUE,
             .plot = FALSE,
             .verbose = FALSE,
             sigma = sigma)
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
  args <- list(input = act,
               t.layer = tm,
               coord.x = "Longitude", coord.y = "Latitude")
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
  # Resample RSP onto spat grid for consistency
  ud_rsp <- terra::project(terra::rast(dbb), terra::crs(spat))
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
                          obs, dlist, algorithm = c("acpf", "acdcpf"),
                          spat, im, win,
                          sigma = spatstat.explore::bw.diggle,
                          overwrite = TRUE) {

  #### (optional) Prior convergence check
  out_file_convergence <-
    here_alg(sim, "patter", algorithm, sim$alg_par, "convergence.rds")
  if (!overwrite && file.exists(out_file_convergence)) {
    return(readRDS(out_file_convergence))
  }

  #### Define simulation arguments
  # Proposal function
  # * Stochastic kicks are up to sim$mobility in length
  # * In directed sampling we account for the discretisation error
  rargs <- list(.shape = sim$shape, .scale = sim$scale,
                .mobility = sim$mobility)
  dargs <- list(.shape = sim$shape, .scale = sim$scale,
                .mobility = sim$mobility + sr)
  # Likelihood functions
  algorithm <- match.arg(algorithm)
  if (algorithm == "acpf") {
    lik <-  list(acs_filter_container = acs_filter_container,
                 pf_lik_ac = pf_lik_ac)
  } else if (algorithm == "acdcpf") {
    lik <-  list(pf_lik_dc = pf_lik_dc,
                 acs_filter_container = acs_filter_container,
                 pf_lik_ac = pf_lik_ac)
  }
  # Record opts
  record <- pf_opt_record(.save = TRUE)

  #### Forward simulation
  t1_pff  <- Sys.time()
  ssf()
  out_pff  <- pf_forward(.obs = obs,
                         .dlist = dlist,
                         .rargs = rargs,
                         .dargs = dargs,
                         .likelihood = lik,
                         .n = sim$n_particles,
                         .record = record,
                         .verbose = FALSE)
  t2_pff   <- Sys.time()
  pff_mins <- mins(t2_pff, t1_pff)
  saveRDS(out_pff$convergence, out_file_convergence)
  if (!out_pff$convergence) {
    return(FALSE)
  }

  #### Backward killer
  pfbk_mins  <- NA_real_
  run_killer <- TRUE
  if (run_killer) {
    t1_pfbk   <- Sys.time()
    out_pfbk  <- pf_backward_killer(.history = out_pff$history,
                                    .record = record,
                                    .verbose = FALSE)
    t2_pfbk   <- Sys.time()
    pfbk_mins <- mins(t2_pfbk, out_pfbk)
  }

  #### Backward sampler
  pfbs_mins  <- NA_real_
  run_sample <- TRUE
  if (run_sampler) {
    t1_pfbs  <- Sys.time()
    ssf()
    dlist$algorithm$sim <- sim
    out_pfbs <- pf_backward_sampler_p(.history = out_pff$history,
                                      .dpropose = pf_dpropose_read,
                                      .obs = obs,
                                      .dlist = dlist,
                                      .dargs = list(),
                                      .record = record,
                                      .verbose = FALSE)
    t2_pfbs  <- Sys.time()
    pfbs_mins <- mins(t2_pfbs, t1_pfbs)
  }

  #### Map arguments
  # Use .discretise = TRUE for speed
  map_args <- list(.map = spat,
                   .im = im,
                   .win = win,
                   .discretise = TRUE,
                   .plot = FALSE,
                   .verbose = FALSE,
                   sigma = sigma)

  #### Mapping (forward run)
  map_args$.coord <- pf_coord(.history = out_pff$history, .bathy = spat)
  t1_udf          <- Sys.time()
  do.call(map_dens, map_args)
  t2_udf          <- Sys.time()
  udf_mins        <- mins(t2_udf, t1_udf)

  #### Mapping (backward killer)
  udk_mins          <- NA_real_
  if (run_killer) {
    map_args$.coord <- NULL
    map_args$.coord <- pf_coord(.history = out_pfbk$history, .bathy = spat)
    t1_udk          <- Sys.time()
    do.call(map_dens, map_args)
    t2_udk          <- Sys.time()
    udk_mins        <- mins(t2_udk, t1_udk)
  }

  #### Mapping (backward sampler)
  uds_mins          <- NA_real_
  if (run_sampler) {
    map_args$.coord <- NULL
    map_args$.coord <- pf_coord(.history = out_pfbs$history, .bathy = spat)
    t1_uds          <- Sys.time()
    do.call(map_dens, map_args)
    t2_uds          <- Sys.time()
    uds_mins        <- mins(t2_uds, t1_uds)
  }

  #### Outputs
  # Timings
  time <- data.table(id = sim$id,
                     pff = pff_mins,
                     pfbk = pfbk_mins,
                     pfbs = pfbs_mins,
                     udf = udf_mins
                     udk = udk_mins,
                     uds = uds_mins)

  qs::qsave(time, here_alg(sim, "patter", algorithm, sim$alg_par, "time.qs"))
  # Particle samples (for checks)
  # qs::qsave(out_pfbk, here_alg(sim, "patter", algorithm, sim$alg_par, "particles-k.qs"))
  # UD
  write_rast(udf, here_alg(sim, "patter", algorithm, sim$alg_par, "ud-f.tif"))
  if (run_killer) {
    write_rast(udk, here_alg(sim, "patter", algorithm, sim$alg_par, "ud-k.tif"))
  }
  if (run_sampler) {
    write_rast(uds, here_alg(sim, "patter", algorithm, sim$alg_par, "ud-s.tif"))
  }
  return(TRUE)
}
