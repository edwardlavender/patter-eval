#' @title Redirect console outputs to file
#' @description This is useful for redirecting progress bars to file for code run via R CMD BATCH.

log_open <- function(log.txt) {
  if (!interactive()) {
    log.con <- file(log.txt, open = "wt")
    sink(log.con, type = "output")
    sink(log.con, type = "message")
  }
  invisible(log.con)
}

log_close <- function(log.con) {
  if (!interactive()) {
    sink(type = "output")
    sink(type = "message")
    close(log.con)
  }
  invisible(NULL)
}
