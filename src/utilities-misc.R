#' @title Utilities

secs <- function(time) {
  as.numeric(lubridate::duration(time))
}

mins <- function(time1, time2) {
  as.numeric(difftime(time1, time2, units = "mins"))
}
