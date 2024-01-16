#' @title UD functions

get_ud_path <- function(sim,
                        spat, path,
                        im, win,
                        sigma = spatstat.explore::bw.diggle,
                        overwrite = TRUE) {
  out_file <- here_alg(sim, "path", "ud.tif")
  if (overwrite | !file.exists(out_file)) {
    ud_path <-
      map_dens(.map = spat,
               .im = im,
               .win = win,
               .coord = path[, .(x, y)],
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
  out_file <- here_alg(sim, "coa", delta_t, "ud.tif")
  acoustics <- dlist$data$acoustics
  if (difftime(max(acoustics$timestamp), min(acoustics$timestamp), units = "mins") <= delta_t) {
    return(NULL)
  }
  out_coa <- coa(.dlist = dlist,
                 .delta_t = delta_t,
                 .plot_weights = FALSE)
  ud_coa  <-
    map_dens(.map = spat,
             .im = im,
             .win = win,
             .coord = out_coa[, .(x = coa_x, y = coa_y)],
             .plot = FALSE,
             .verbose = FALSE,
             sigma = sigma)
  write_rast(ud_coa, out_file)
  ud_coa
}

get_ud_rsp <- function() {

}

get_ud_patter <- function(sim,
                          obs, dlist, algorithm = c("acpf", "acdcpf"),
                          spat, im, win,
                          sigma = spatstat.explore::bw.diggle) {

  #### Define simulation arguments
  # Proposal function
  margs <- list(.shape = sim$shape, .scale = sim$scale, .mobility = sim$mobility)
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
  t1_pff <- Sys.time()
  ssf()
  out_pff <- pf_forward(.obs = obs,
                        .dlist = dlist,
                        .rargs = margs,
                        .dargs = margs,
                        .likelihood = lik,
                        .n = sim$n_particles,
                        .record = record,
                        .verbose = TRUE)
  t2_pff <- Sys.time()
  if (!out_pff$convergence) {
    return(0)
  }

  #### Backward killer
  t1_pfbk <- Sys.time()
  out_pfbk <- pf_backward_killer(out_pff$history,
                                 .record = record,
                                 .verbose = FALSE)
  t2_pfbk <- Sys.time()

  #### Backward sampler
  # TO DO

  #### Smoothing (backward killer)
  t1_udk   <- Sys.time()
  udk <- map_dens(.map = spat,
                  .im = im,
                  .win = win,
                  .coord = pf_coord(.history = out_pfbk$history, .bathy = spat),
                  .plot = FALSE,
                  .verbose = FALSE,
                  sigma = sigma)
  t2_udk <- Sys.time()

  #### Smoothing (backward sampler)
  # TO DO

  #### Outputs
  # Timings
  time <- data.table(id = sim$id,
                     pff = mins(t2_pff, t1_pff),
                     pfbk = mins(t2_pfbk, t1_pfbk),
                     pfbs = NA,
                     ud = mins(t2_udk, t1_udk))
  qs::qsave(time, here_alg(sim, "patter", algorithm, sim$alg_par, "time.qs"))
  # Particle samples (for checks)
  qs::qsave(out_pfbk, here_alg(sim, "patter", algorithm, sim$alg_par, "particles-k.qs"))
  # UD
  write_rast(udk, here_alg(sim, "patter", algorithm, sim$alg_par, "ud-k.tif"))
  return(1)
}
