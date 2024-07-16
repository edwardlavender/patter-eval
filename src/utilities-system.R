#' Get the user name
user <- function() {
  unname(Sys.info()["login"])
}

#' Is the username lavended?
lavended <- function() {
  user() == "lavended"
}
