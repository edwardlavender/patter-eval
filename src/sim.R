#' @title Data generation/simulation helpers

gen_depth <- function(.xy) {
  250 - 100 * cos(sqrt(.xy$x^2 + .xy$y^2) / (500 * 2 * pi))
}

sim_depth <- function(.seabed) {
  z <- runif(length(.seabed), .seabed - 5, .seabed + 5)
  z[z < 0] <- 0
  z
}
