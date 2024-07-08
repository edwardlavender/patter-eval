ddetlogistic <- function(.x, .alpha, .beta, .gamma) {
  ifelse(.x <= .gamma, plogis(.alpha + .beta * .x), 0)
}

rtruncgamma <- function(.n, .shape, .scale, .mobility) {
  truncdist::rtrunc(.n, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}

dtruncgamma <- function(.x, .shape, .scale, .mobility, .scale_to_one = FALSE) {
  y <- truncdist::dtrunc(.x, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
  if (.scale_to_one) {
    y <- y / max(y)
  }
  y
}

check_multithreading <- function(.multithread = c("R", "Julia")) {
  .multithread <- match.arg(.multithread)
  nthreads <- JuliaCall::julia_eval("Threads.nthreads()")
  if (.multithread == "R" & nthreads != 1L) {
    stop("R is multi-threaded but Julia is not using one thread!")
  }
  if (.multithread == "Julia" & nthreads == 1L) {
    stop("Julia is multi-threaded but Julia is only using one thread!")
  }
  nthreads
}
