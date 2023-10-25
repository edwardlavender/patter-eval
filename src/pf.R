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
