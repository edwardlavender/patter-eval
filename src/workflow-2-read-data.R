#' @title Read data functions

read_array <- function(sim) {
  here_input("arrays", sim$array_type, sim$array_realisation, "array.qs") |>
    qs::qread()
}

read_acoustics <- function(sim) {
  here_input("acoustics",
             sim$combination,
             sim$array_type, sim$array_realisation,
             sim$path_realisation,
             "acoustics.qs") |>
    qs::qread()
}

read_path <- function(sim) {
  here_input("paths",
             sim$combination,
             sim$array_type, sim$array_realisation,
             sim$path_realisation,
             "path.qs") |>
    qs::qread()
}

read_overlaps <- function(sim) {
  here_input("ac", sim$array_type, sim$array_realisation, sim$gamma, "overlaps.rds") |>
    readRDS()
}

read_kernels <- function(sim) {
  here_input("ac", sim$array_type, sim$array_realisation, sim$gamma,
             sim$alpha, sim$beta, "kernels.rds") |>
    readRDS() |>
    unwrapr()
}

read_ud_path <- function() {

}
