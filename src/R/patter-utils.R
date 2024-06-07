ddetlogistic <- function(.x, .alpha, .beta, .gamma) {
  ifelse(.x <= .gamma, plogis(.alpha + .beta * .x), 0)
}

rtruncgamma <- function(.n, .shape, .scale, .mobility) {
  truncdist::rtrunc(.n, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}

dtruncgamma <- function(.x, .shape, .scale, .mobility) {
  truncdist::dtrunc(.x, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}
