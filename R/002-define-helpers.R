###########################
###########################
#### define_helpers.R

#### Aims
# 1) Define helper functions

#### Prerequisites
# 1) NA


###########################
###########################
#### Define functions

#' @title Utilities

mins <- function(time1, time2) {
  as.numeric(difftime(time1, time2, units = "mins"))
}

#' @title `here::here()` helpers

here_input <- function(...) {
  dv::here_data("sims", "input", ...)
}

here_output <- function(...) {
  dv::here_data("sims", "output", ...)
}

here_alg <- function(...) {
  dv::here_data("sims", "output", "algorithm", ...)
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
  l <- list(grid = grid, elm = list(g = grid, b = grid, c = NULL), d = 1)
  wrapr(l)
  unwrapr(wrapr(l))
}

#' @title Data generation/simulation helpers

gen_depth <- function(.xy) {
  250 - 100 * cos(sqrt(.xy$x^2 + .xy$y^2) / (500 * 2 * pi))
}

sim_depth <- function(.seabed) {
  z <- runif(length(.seabed), .seabed - 5, .seabed + 5)
  z[z < 0] <- 0
  z
}

#' @title pf helpers
#' @details
#' * `pf_coords()` extracts particle coordinates from the `history` element of a `pf` object

pf_coords <- function(.history, .bathy) {
  .history |>
    rbindlist() |>
    select(cell = cell_now) |>
    mutate(xy = terra::xyFromCell(.bathy, cell),
           x = xy[, 1],
           y = xy[, 2]) |>
    select(x, y) |>
    as.data.table()
}

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


#### End of code.
###########################
###########################
