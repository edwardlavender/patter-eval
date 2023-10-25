#' @title Workflows

workflow_path <- function(sim, grid) {
  path    <- read_path(sim)
  ud_path <- get_ud_path(sim = sim, grid = grid, path = path)
  if (is.null(ud_path)) {
    stop(paste0("`pf_dens()` failed for simulation ID ", sim$id, "."))
  }
  NULL
}

workflow_coa <- function(sim, grid) {
  acoustics <- read_acoustics(sim)
  get_ud_coa(sim = sim, grid = grid, acoustics = acoustics, delta_t = "30 mins")
  get_ud_coa(sim = sim, grid = grid, acoustics = acoustics, delta_t = "120 mins")
  NULL
}

workflow_patter <- function(sim, grid) {
  # Read data
  acoustics <- read_acoustics(sim)
  path      <- read_path(sim)
  array     <- read_array(sim)
  overlaps  <- read_overlaps(sim)
  kernels   <- read_kernels(sim)
  # Algorithm preparation
  # * Note that sim$mobility has been adjusted for grid resolution (see sim-data.R)
  obs <- acs_setup_obs(acoustics,
                       .archival = path,
                       .step = paste(sim$step, "mins"),
                       .mobility = sim$mobility)
  obs[, depth_shallow := obs$depth - 5]
  obs[, depth_deep := obs$depth + 5]
  # ACPF algorithm
  get_ud_patter(sim = sim,
                obs = obs, grid = grid,
                array = array, overlaps = overlaps, kernels = kernels, update_ac = NULL)
  # ACDCPF algorithm
  get_ud_patter(sim = sim,
                obs = obs, grid = grid,
                array = array, overlaps = overlaps, kernels = kernels, update_ac = update_ac)
  NULL
}
