#' @title Workflows

workflow_path <- function(sim, spat, im, win) {
  path    <- read_path(sim)
  stopifnot(nrow(path) != 0L)
  ud_path <- get_ud_path(sim = sim, spat = spat,
                         path = path, im = im, win = win)
  if (is.null(ud_path)) {
    stop(paste0("`map_dens()` failed for simulation ID ", sim$id, "."))
  }
  NULL
}

workflow_coa <- function(sim, spat, im, win) {
  dlist <- read_dlist(sim)
  get_ud_coa(sim = sim, spat = spat,
             dlist = dlist, delta_t = "30 mins",
             im = im, win = win)
  get_ud_coa(sim = sim, spat = spat,
             dlist = dlist, delta_t = "120 mins",
             im = im, win = win)
  NULL
}

workflow_rsp <- function(sim, spat, spat_ll_dbb, tm) {
  s1 <- get_ud_rsp(sim = sim, spat = spat, spat_ll_dbb = spat_ll_dbb, tm = tm,
                   type = "default")
  s2 <- get_ud_rsp(sim = sim, spat = spat, spat_ll_dbb = spat_ll_dbb, tm = tm,
                   er.ad = "custom")
  data.frame(row = sim$row, rsp_1 = s1, rsp_2 = s2)
}

workflow_patter <- function(sim, spat, im, win) {
  # Define data list
  # (optional) TO DO
  # * Pre-prepare dlist as required with detection_overlaps & detection_kernels elements
  dlist <- read_dlist(sim)
  dlist$spatial$bathy <- spat
  dlist$algorithm$detection_overlaps <- read_overlaps(sim)
  dlist$algorithm$detection_kernels  <- read_kernels(sim)
  # Define observation timeline
  # * Note that sim$mobility has been adjusted for spat resolution (see sim-data.R)
  obs <- pf_setup_obs(.dlist = dlist,
                      .trim = TRUE,
                      .step = paste(sim$step, "mins"),
                      .mobility = sim$mobility,
                      .receiver_range = dlist$data$moorings$receiver_range[1])
  obs[, depth_shallow := obs$depth - 5]
  obs[, depth_deep := obs$depth + 5]
  # ACPF algorithm
  success_acpf <- get_ud_patter(sim = sim,
                                obs = obs, dlist = dlist, algorithm = "acpf",
                                spat = spat, im = im, win = win)
  # ACDCPF algorithm
  success_acdcpf <- get_ud_patter(sim = sim,
                                  obs = obs, dlist = dlist, algorithm = "acdcpf",
                                  spat = spat, im = im, win = win)
  data.table(row = sim$row, acpf = success_acpf, acdcpf = success_acdcpf)
}
