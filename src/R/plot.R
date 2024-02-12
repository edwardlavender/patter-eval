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

#' @title Add a smoother to a plot

add_smoother <- function(x, y, ...) {
  if (all(is.na(y))) {
    return(NULL)
  }
  mod  <- mgcv::gam(y ~ s(x, k = 3))
  data <- data.frame(x = seq(min(x, na.rm = TRUE),
                             max(x, na.rm = TRUE),
                             length.out = 100))
  data$y <- predict(mod, newdata = data, type = "response")
  lines(data$x, data$y, ...)
}
