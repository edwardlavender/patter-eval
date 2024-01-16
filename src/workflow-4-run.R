#' @title Workflows

workflow_path <- function(sim, grid, im, win) {
  path    <- read_path(sim)
  stopifnot(nrow(path) != 0L)
  ud_path <- get_ud_path(sim = sim, grid = grid,
                         path = path, im = im, win = win)
  if (is.null(ud_path)) {
    stop(paste0("`map_dens()` failed for simulation ID ", sim$id, "."))
  }
  NULL
}

workflow_coa <- function(sim, grid, im, win) {
  dlist <- read_dlist(sim)
  get_ud_coa(sim = sim, grid = grid,
             dlist = dlist, delta_t = "30 mins",
             im = im, win = win)
  get_ud_coa(sim = sim, grid = grid,
             dlist = dlist, delta_t = "120 mins",
             im = im, win = win)
  NULL
}

workflow_patter <- function(sim, grid, im, win) {
  # Define data list
  # (optional) TO DO
  # * Pre-prepare dlist as required with detection_overlaps & detection_kernels elements
  dlist <- read_dlist(sim)
  dlist$spatial$bathy <- grid
  dlist$algorithm$detection_overlaps <- read_overlaps(sim)
  dlist$algorithm$detection_kernels  <- read_kernels(sim)
  # Define observation timeline
  # * Note that sim$mobility has been adjusted for grid resolution (see sim-data.R)
  obs <- pf_setup_obs(.dlist = dlist,
                      .trim = TRUE,
                      .step = paste(sim$step, "mins"),
                      .mobility = sim$mobility,
                      .receiver_range = dlist$data$moorings$receiver_range[1])
  obs[, depth_shallow := obs$depth - 5]
  obs[, depth_deep := obs$depth + 5]
  # ACPF algorithm
  get_ud_patter(sim = sim,
                obs = obs, dlist = dlist, algorithm = "acpf",
                grid = grid, im = im, win = win)
  # ACDCPF algorithm
  get_ud_patter(sim = sim,
                obs = obs, dlist = dlist, algorithm = "acdcpf",
                grid = grid, im = im, win = win)
  NULL
}
