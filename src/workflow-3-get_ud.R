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
  stopifnot(delta_t %in% c("30 mins", "120 mins"))
  delta_t <- ifelse(delta_t == "120 mins", "2 hours", delta_t)
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
                          timeline, acoustics, archival = NULL, model_move,
                          algorithm = c("acpf", "acdcpf"),
                          spat, win, sigma = spatstat.explore::bw.diggle,
                          overwrite = TRUE,
                          performance = TRUE,
                          test = FALSE) {

  #### (optional) TO DO
  # Split get_ud_patter() into two functions
  # * get_ud_patter_0() runs the algorithms (and saves particles);
  # * get_ud_patter_1() estimates the UDs only, reading particles from file
  # * Then we can easily re-run get_ud_patter_1()
  # * ... with the same particles
  # * ... but different estimation (e.g., sigma) parameters, if needed
  # * This option may also be faster:
  # - multi-thread algorithm runs in Julia
  # - multi-thread UD estimation in R
  # - through writing to disk also has a (high) speed penalty
  # * But requires much more storage space (see run-patter.R)
  # * Currently, we do not save particle samples to minimise storage costs

  #### Check inputs
  if (test) {
    warning("Test is true!", immediate. = TRUE)
  }

  #### (optional) Prior run checks
  # Return FALSE, if previous convergence failure
  out_file_convergence <-
    here_alg(sim, "patter", algorithm, sim$alg_par, "convergence.rds")
  if (!overwrite && file.exists(out_file_convergence)) {
    return(readRDS(out_file_convergence))
  }

  #### Baseline filter args
  args <- list(.map = spat,
               .timeline = timeline,
               .state = "StateXY",
               .xinit = NULL,
               .xinit_pars = list(mobility = sim$mobility),
               .yobs = NULL,
               .model_obs = NULL,
               .model_move = model_move,
               .n_particle = NULL,
               .n_resample = NULL,
               .verbose = FALSE
               )

  #### Forward filter arguments
  algorithm <- match.arg(algorithm)
  args      <- assemble_args(sim = sim, args = args,
                             acoustics = acoustics, archival = archival, direction = "forward",
                             algorithm = algorithm)

  #### Forward filter implementation
  # out_file_pff <- here_alg(sim, "patter", algorithm, sim$alg_par, "out_pff.qs")
  t1_pff       <- Sys.time()
  out_pff      <- do.call(pf_filter, args, quote = TRUE)
  t2_pff       <- Sys.time()
  pff_mins     <- mins(t2_pff, t1_pff)
  # qs::qsave(out_pff, out_file_pff)
  if (!out_pff$convergence) {
    saveRDS(FALSE, out_file_convergence)
    return(FALSE)
  }

  #### Backward filter arguments
  args <- assemble_args(sim = sim, args = args,
                        acoustics = acoustics, archival = archival, direction = "backward",
                        algorithm = algorithm)

  #### Backward filter implementation
  # out_file_pfb <- here_alg(sim, "patter", algorithm, sim$alg_par, "out_pfb.qs")
  t1_pfb          <- Sys.time()
  out_pfb         <- do.call(pf_filter, args, quote = TRUE)
  t2_pfb          <- Sys.time()
  pfb_mins        <- mins(t2_pfb, t1_pfb)
  if (!out_pfb$convergence) {
    saveRDS(FALSE, out_file_convergence)
    return(FALSE)
  }
  # qs::qsave(out_pfb, out_file_pfb)
  saveRDS(TRUE, out_file_convergence)

  #### Mapping
  # This code does not need to be run if we just want to check convergence
  map <- TRUE
  if (map) {

    #### Backward smoother
    # (optional) TO DO
    # * Improve speed by pre-defining box in Julia
    # out_file_smo <- here_alg(sim, "patter", algorithm, sim$alg_par, "out_smo.qs")
    t1_pfbs      <- Sys.time()
    out_smo      <- pf_smoother_two_filter(.map = spat,
                                           .mobility = sim$mobility,
                                           .n_particle = 1000L,
                                           .verbose = FALSE)
    t2_pfbs      <- Sys.time()
    pfbs_mins    <- mins(t2_pfbs, t1_pfbs)
    # qs::qsave(out_smo, out_file_smo)

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
    # We only implement this for the performance simulations (for speed)
    if (performance) {
      t1_udf   <- Sys.time()
      udf      <- do.call(map_dens, map_args)$ud
      t2_udf   <- Sys.time()
      udf_ok   <- !is.null(udf)
      if (udf_ok) {
        udf_mins <- mins(t2_udf, t1_udf)
      } else {
        udf_mins <- NA_real_
      }
    } else {
      udf      <- NULL
      udf_mins <- NA_real_
    }

    #### Mapping (backward run)
    # (For speed, this is not currently implemented)

    #### Mapping (smoother)
    map_args$.coord <- out_smo$states
    t1_uds          <- Sys.time()
    uds             <- do.call(map_dens, map_args)$ud
    t2_uds          <- Sys.time()
    uds_ok          <- !is.null(uds)
    if (uds_ok) {
      uds_mins      <- mins(t2_uds, t1_uds)
    } else {
      uds_mins <- NA_real_
    }

    #### Outputs
    time <- data.table(id = sim$id,
                       pff = pff_mins,
                       pfbs = pfbs_mins,
                       udf = udf_mins,
                       uds = uds_mins)
    qs::qsave(time, here_alg(sim, "patter", algorithm, sim$alg_par, "time.qs"))
    if (performance) {
      write_rast(udf, here_alg(sim, "patter", algorithm, sim$alg_par, "ud-f.tif"))
    }
    write_rast(uds, here_alg(sim, "patter", algorithm, sim$alg_par, "ud-s.tif"))

  } else {
    stop("`map` is FALSE! This is only for testing convergence!")
  }

  #### Return outputs
  if (test) {

    return(
      list(
        particles = list(pff = out_pff,
                         pfb = out_pfb,
                         smo = out_smo),
        uds = list(udf = udf,
                   uds = uds)
      )
    )

  } else {

    return(TRUE)

  }

}
