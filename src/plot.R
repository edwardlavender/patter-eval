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

# Add points to a map
spatPoints <- function(x, y, ...) {
  # Use terra::plot() to add points to a map for safety
  # (with points(), points may be misplaced when par() is set)
  terra::plot(terra::vect(cbind(x, y), crs = "EPSG:32629"),
              ..., add = TRUE)
}

# Add paths to a map
spatPath <- function(x, y, ...) {
  n <- length(x)
  p1 <- cbind(x, y)[1:(n - 1), ]
  p2 <- cbind(x, y)[2:n, ]
  x <- terra::vect(p1, crs = "EPSG:32629")
  y <- terra::vect(p2, crs = "EPSG:32629")
  terra::lines(x, y, arrows = TRUE, ...)
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
