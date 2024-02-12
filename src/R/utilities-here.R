#' @title `here::here()` helpers

here_input <- function(...) {
  dv::here_data("sims", "input", ...)
}

here_output <- function(...) {
  dv::here_data("sims", "output", ...)
}

here_alg <- function(sim, ...) {
  here_output("run",
              sim$combination, sim$array_type, sim$array_realisation,
              sim$path_realisation, ...)
}
