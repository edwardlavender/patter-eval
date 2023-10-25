#' @title Utilities

mins <- function(time1, time2) {
  as.numeric(difftime(time1, time2, units = "mins"))
}
