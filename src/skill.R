#' @title skill_ wrappers

skill_by_alg <- function(.uds, .ud_path, .f = patter::skill_mb) {
  sapply(.uds, \(.ud) {
    if (is.null(.ud)) {
      return(NA_real_)
    }
    .f(.ud, .ud_path)
  }) |>
    as.numeric()
}

skill_mew <- function(.obs, .mod) {
  if (is.null(.obs) | is.null(.mod)) {
    return(NA_real_)
  }
  skill_me(.obs = .obs, .mod = .mod)
}
