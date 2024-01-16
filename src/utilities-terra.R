
#' @title Read/write SpatRasters

read_rast <- function(x) {
  if (file.exists(x)) {
    terra::rast(x)
  } else {
    NULL
  }
}

write_rast <- function(x, filename, overwrite = TRUE) {
  if (!is.null(x)) {
    terra::writeRaster(x, filename, overwrite = overwrite)
  }
}

#' @title Quick plots from file
terra_qplot <- function(file) {
  file |>
    terra::rast() |>
    terra::plot()
}

#' @title Wrap/unwrap SpatRasters in a list recusively

wrapr <- function(.x) {
  f <- function(.e) {
    if (inherits(.e, "SpatRaster")) {
      terra::wrap(.e)
    } else {
      .e
    }
  }
  rapply(.x, f, how = "replace")
}

unwrapr <- function(.x) {
  f <- function(.e) {
    if (inherits(.e, "PackedSpatRaster")) {
      terra::unwrap(.e)
    } else {
      .e
    }
  }
  rapply(.x, f, how = "replace")
}

#' @examples
if (FALSE) {
  l <- list(spat = spat, elm = list(g = spat, b = spat, c = NULL), d = 1)
  wrapr(l)
  unwrapr(wrapr(l))
}
