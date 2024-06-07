#' @title Workflows

# TO DO
# To revise in line with _ud() function updates!

workflow_path <- function(sim, spat, win) {
  path    <- read_path(sim)
  stopifnot(nrow(path) != 0L)
  ud_path <- get_ud_path(sim = sim, spat = spat,
                         path = path, im = im, win = win)
  if (is.null(ud_path)) {
    stop(paste0("`map_dens()` failed for simulation ID ", sim$id, "."))
  }
  NULL
}

workflow_coa <- function(sim, spat, win) {
  dlist <- read_dlist(sim)
  s1 <- get_ud_coa(sim = sim, spat = spat,
                   dlist = dlist, delta_t = "30 mins",
                   im = im, win = win)
  s2 <- get_ud_coa(sim = sim, spat = spat,
                   dlist = dlist, delta_t = "120 mins",
                   im = im, win = win)
  data.frame(row = sim$row, coa_1 = s1, coa_2 = s2)
}

workflow_rsp <- function(sim, spat, spat_ll_dbb, tm) {
  s1 <- get_ud_rsp(sim = sim, spat = spat, spat_ll_dbb = spat_ll_dbb, tm = tm,
                   type = "default")
  s2 <- get_ud_rsp(sim = sim, spat = spat, spat_ll_dbb = spat_ll_dbb, tm = tm,
                   type = "custom")
  data.frame(row = sim$row, rsp_1 = s1, rsp_2 = s2)
}

workflow_patter <- function(sim, spat, win) {

  # Define observations
  # TO DO
  # dlist <- read_dlist(sim)

  # ACPF algorithm
  success_acpf <- get_ud_patter(sim = sim,
                                obs = obs, dlist = dlist, algorithm = "acpf",
                                spat = spat, im = im, win = win)

  # ACDCPF algorithm
  success_acdcpf <- get_ud_patter(sim = sim,
                                  obs = obs, dlist = dlist, algorithm = "acdcpf",
                                  spat = spat, im = im, win = win)
  # Outputs
  data.table(row = sim$row, acpf = success_acpf, acdcpf = success_acdcpf)
}
