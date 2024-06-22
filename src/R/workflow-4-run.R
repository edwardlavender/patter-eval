#' @title Workflows

# spat <- terra::unwrap(spatw)

workflow_path <- function(sim, spat, win) {
  path    <- read_path(sim)
  stopifnot(nrow(path) != 0L)
  ud_path <- get_ud_path(sim = sim, path = path,
                         spat = spat, win = win)
  # terra::plot(ud_path)
  # points(path$x, path$y, pch = ".")
  if (is.null(ud_path)) {
    stop(paste0("`map_dens()` failed for simulation ID ", sim$id, "."))
  }
  NULL
}

workflow_coa <- function(sim, spat, win) {
  detections <- read_detections(sim)
  s1 <- get_ud_coa(sim = sim,
                   detections = detections, delta_t = "30 mins",
                   spat = spat, win = win)
  s2 <- get_ud_coa(sim = sim,
                   detections = detections, delta_t = "120 mins",
                   spat = spat, win = win)
  data.frame(row = sim$row, coa_1 = s1, coa_2 = s2)
}

workflow_rsp <- function(sim, spat, spat_ll_dbb, tm) {
  # spat = terra::unwrap(spatw)
  # spat_ll_dbb = terra::unwrap(spat_ll_dbbw)
  s1 <- get_ud_rsp(sim = sim, spat = spat, spat_ll_dbb = spat_ll_dbb, tm = tm,
                   type = "default")
  s2 <- get_ud_rsp(sim = sim, spat = spat, spat_ll_dbb = spat_ll_dbb, tm = tm,
                   type = "custom")
  data.frame(row = sim$row, rsp_1 = s1, rsp_2 = s2)
}

workflow_patter <- function(sim, spat, win, test = FALSE) {

  # Read data
  acoustics  <- read_acoustics(sim)
  archival   <- read_archival(sim)

  # Algorithm components
  t1         <- min(acoustics$timestamp)
  tT         <- max(acoustics$timestamp)
  timeline   <- seq(t1, tT, by = "2 mins")
  model_move <-
    move_xy(dbn_length =
              glue::glue("truncated(Gamma({sim$shape},
                                    {sim$scale}),
                              upper = {sim$mobility})"),
            dbn_angle = "Uniform(-pi, pi)")

  # ACPF algorithm
  success_acpf <- get_ud_patter(sim = sim,
                                timeline = timeline,
                                acoustics = acoustics,
                                archival = NULL,
                                model_move = model_move,
                                algorithm = "acpf",
                                spat = spat, win = win,
                                test = test)

  # ACDCPF algorithm
  success_acdcpf <- get_ud_patter(sim = sim,
                                  timeline = timeline,
                                  acoustics = acoustics,
                                  archival = archival,
                                  model_move = model_move,
                                  algorithm = "acdcpf",
                                  spat = spat, win = win,
                                  test = test)

  # Outputs
  if (test) {
    list(
      input = list(timeline = timeline, path = read_path(sim), acoustics = acoustics, archival = archival),
      output = list(acpf = success_acpf, acdcpf = success_acdcpf)
    )
  } else {
    data.table(row = sim$row, acpf = success_acpf, acdcpf = success_acdcpf)
  }

}
