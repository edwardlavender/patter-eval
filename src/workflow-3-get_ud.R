#' @title UD functions

get_ud_path <- function(sim, grid, path, sigma = spatstat.explore::bw.diggle) {
  out_file <- here_alg(sim, "path", "ud.tif")
  ud_path <-
    pf_dens(.xpf = grid,
            .coord = path[, .(x, y)],
            .plot = FALSE,
            .verbose = FALSE,
            sigma = sigma)
  write_rast(ud_path, out_file)
  ud_path
}

get_ud_coa <- function(sim, grid, acoustics, delta_t, sigma = spatstat.explore::bw.diggle) {
  out_file <- here_alg(sim, "coa", delta_t, "ud.tif")
  if (difftime(max(acoustics$timestamp), min(acoustics$timestamp), units = "mins") <= delta_t) {
    return(NULL)
  }
  out_coa <- coa(acoustics,
                 .delta_t = delta_t,
                 .plot_weights = FALSE)
  ud_coa  <- pf_dens(.xpf = grid,
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
                           obs, grid,
                           array, overlaps, kernels, update_ac = NULL,
                           sigma = spatstat.explore::bw.diggle) {

  #### Forward simulation
  t1_pff <- Sys.time()
  set.seed(1)
  out_pff <- pf_forward_2(obs,
                          .bathy = grid,
                          .moorings = array,
                          .detection_overlaps = overlaps,
                          .detection_kernels = kernels,
                          .update_ac = update_ac,
                          .kick = pf_kick,
                          .shape = sim$shape, .scale = sim$scale, .mobility = sim$mobility,
                          .n = sim$n_particles,
                          .save_history = TRUE,
                          .verbose = FALSE)
  t2_pff <- Sys.time()
  if (!inherits(out_pff, "pf")) {
    return(0)
  }

  #### Backward pass
  t1_pfb <- Sys.time()
  out_pfb <- pf_backward(out_pff$history,
                         .save_history = TRUE,
                         .verbose = FALSE)
  t2_pfb <- Sys.time()

  #### Smoothing
  t1_ud   <- Sys.time()
  out_pfp <- pf_coords(out_pfb$history, grid)
  ud_acpf <- pf_dens(.xpf = grid,
                     .coord = out_pfp,
                     .plot = FALSE,
                     .verbose = FALSE,
                     sigma = sigma)
  t2_ud <- Sys.time()

  #### Outputs
  # Timings
  folder <- ifelse(is.null(update_ac), "acpf", "acdcpf")
  time <- data.table(id = sim$id,
                     pff = mins(t2_pff, t1_pff),
                     pfb = mins(t2_pfb, t1_pfb),
                     ud = mins(t2_ud, t1_ud))
  qs::qsave(time, here_alg(sim, "patter", folder, sim$alg_par, "time.qs"))
  # Particle samples (for checks)
  qs::qsave(out_pfb, here_alg(sim, "patter", folder, sim$alg_par, "particles.qs"))
  # UD
  write_rast(ud_acpf, here_alg(sim, "patter", folder, sim$alg_par, "ud.tif"))
  return(1)
}

update_ac <- function(.particles, .bathy, .obs, .t, ...) {
  .particles$bathy <- terra::extract(.bathy, as.matrix(.particles[, c("x_now", "y_now")]))
  (.particles$bathy >= .obs$depth_shallow[.t] & .particles$bathy <= .obs$depth_deep[.t]) + 0
}
