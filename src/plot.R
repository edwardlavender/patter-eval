#' @title Map SpatRasters with out axis labels

# Make a map without axes
spatMap <- function(x, ...) {
  b <- terra::ext(x)[] |> as.numeric()
  terra::plot(x, xlim = b[1:2], ylim = b[3:4], axes = FALSE, ...)
}

# Add box
# * This needs to be implemented after adding points to a plot
# * ... or points aren't added properly, presumably due to how terra
# * ... sets the plotting window
spatAxes <- function(x) {
  b <- terra::ext(x)[] |> as.numeric()
  axis(side = 1, pos = b[1], tck = 0.025, labels = FALSE)
  axis(side = 2, pos = b[3], tck = 0.025, labels = FALSE)
  axis(side = 3, pos = b[4], tck = 0.025, labels = FALSE)
  axis(side = 4, pos = b[2], tck = 0.025, labels = FALSE)
}
