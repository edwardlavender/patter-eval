#' @title Pre-calculate movement densities around a point

spatDens <- function(.spat, .xy, .shape, .scale, .mobility, .file) {

  # Note that .mobility has been adjusted to account for grid resolution

  # Define reachable zone
  zone <-
    .xy |>
    terra::vect(crs = crs) |>
    terra::buffer(width = .mobility, quadsegs = 1e3L) |>
    sf::st_as_sf()

  # Define probability density of movements into reachable cells
  # Identify reachable cells
  exactextractr::exact_extract(.spat, zone,
                               include_cell = TRUE,
                               include_cols = NULL,
                               include_xy = TRUE)[[1]] |>
    # Compute densities between cells
    mutate(cell = as.integer(.data$cell),
           pt_x = .xy[, 1],
           pt_y = .xy[, 2],
           len = patter:::dist_2d(.data$pt_x, .data$pt_y, .data$x, .data$y),
           dens = dtruncgamma(.data$len,
                              .shape = .shape,
                              .scale = .scale,
                              .mobility = .mobility)) |>
    # Retain only essential columns
    select("cell", "dens") |>
    as.data.table() |>
    # Write to file
    qs::qsave(.file)

  NULL
}
