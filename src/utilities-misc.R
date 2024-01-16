#' @title Utilities

secs <- function(time1, time2) {
  as.numeric(difftime(time1, time2, units = "secs"))
}

mins <- function(time1, time2) {
  as.numeric(difftime(time1, time2, units = "mins"))
}
