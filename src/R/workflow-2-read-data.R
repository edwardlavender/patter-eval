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

read_detections <- function(sim) {
  here_input("acoustics",
             sim$combination,
             sim$array_type, sim$array_realisation,
             sim$path_realisation,
             "detections.qs") |>
    qs::qread()
}

read_archival <- function(sim) {
  here_input("acoustics",
             sim$combination,
             sim$array_type, sim$array_realisation,
             sim$path_realisation,
             "archival.qs") |>
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

read_actel <- function(sim) {
  here_input("actel",
             sim$combination,
             sim$array_type, sim$array_realisation,
             sim$path_realisation,
             "actel.qs") |>
    qs::qread()
}
