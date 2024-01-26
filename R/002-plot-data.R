#########################
#########################
#### plot-data.R

#### Aims
# 1) Plot simulated data

#### Prerequisites
# 1) NA


#########################
#########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)
library(prettyGraphics)

#### Load data
detection_pars <- readRDS(here_input("detection_pars.rds"))
path_pars      <- readRDS(here_input("path_pars.rds"))


#########################
#########################
#### Process models

#### Set up figure
# png()
png(here_fig("systems.png"),
    height = 4, width = 8, units = "in", res = 300)
# Window
pp <- par(mfrow = c(1, 2),
          oma = c(2, 2, 2, 3),
          mar = c(2, 2, 2, 2))
# Parameters
lwd <- 2
cols <- c("orange", "royalblue")
cex.axis <- 1.2
cex.label <- 1

#### Plot movement (step length) densities
# Limits
xlim <- c(0, 800)
ylim <- c(0, 0.008)
# Base plot
pretty_plot(0,
            xlim = xlim, ylim = ylim,
            xlab = "", ylab = "",
            type = "n")
# Densities
for (i in seq_len(nrow(detection_pars))) {
  x <- seq(0, path_pars$mobility[i] + 1, length.out = 1e3L)
  y <- dtruncgamma(x,
                   .shape = path_pars$shape[i],
                   .scale = path_pars$scale[i],
                   .mobility = path_pars$mobility[i])
  lines(x, y, col = cols[i], lwd = lwd)
}
# Axes
mtext(side = 3, "A", cex = cex.label, font = 2, adj = 0.05)

#### Plot detection densities
# Limits
xlim <- c(0, 1000)
ylim <- c(0, 1)
# Base plot
pretty_plot(0,
            xlim = xlim, ylim = ylim,
            xlab = "", ylab = "",
            type = "n")
# Densities
for (i in seq_len(nrow(detection_pars))) {
  x <- seq(0, detection_pars$gamma[i] + 1, length.out = 1e3L)
  pr <- ddetlogistic(x,
                     .alpha = detection_pars$alpha[i],
                     .beta = detection_pars$beta[i],
                     .gamma = detection_pars$gamma[i])
  y <- dbinom(1L, size = 1L, prob = pr)
  lines(x, y, col = cols[i], lwd = lwd)
}
# Axes
mtext(side = 3, "B", cex = cex.label, font = 2, adj = 0.05)

#### Add global objects
legend("topright",
       legend = c("System 1", "System 2"),
       col = cols,
       lwd = lwd,
       box.lwd = 0.5)
mtext(side = 1, "Distance (m)", cex = cex.axis, outer = TRUE, line = 1)
mtext(side = 2, "Probability density", cex = cex.axis, outer = TRUE, line = 1)

par(pp)
dev.off()


#########################
#########################
#### Arrays

png(here_fig("bathymetry.png"),
    height = 3, width = 3, units = "in", res = 300)
# terra::plot(spat)
spatMap(spat, legend = TRUE,
        col = viridis::mako(255),
        range = c(150, 350),
        cex.axis = 5)
dev.off()

png(here_fig("arrays.png"),
    height = 6, width = 6, units = "in", res = 300)
pp <- par(mfcol = c(5, 4),
          oma = c(2, 1, 1, 1),
          mar = c(0, 0, 0, 0))

lapply(arrays, function(m) {
  # m <- arrays[[1]]
  # Base map
  spatMap(spat, legend = FALSE, col = viridis::mako(255),
          mar = NA, reset = TRUE)
  # Add receivers
  points(m$receiver_easting, m$receiver_northing,
         pch = 21,
         bg = scales::alpha("black", 0.75),
         col = scales::alpha("black", 0.75),
         cex = 0.5)
  # Add axes
  spatAxes(spat)
})
par(pp)
dev.off()


#### End of code.
#########################
#########################
