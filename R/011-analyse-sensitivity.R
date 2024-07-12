#########################
#########################
#### analyse-sensitivity.R

#### Aims:
# (1) Analyses algorithm sensitivity

#### Prerequisites
# 1) Implement simulations


#########################
#########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(collapse)
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(prettyGraphics)
library(testthat)
library(tictoc)
dv::src()

#### Load data
spat <- terra::rast(here_input("spat.tif"))
sims_performance   <- readRDS(here_input("sims-performance.rds"))
sims_senstivity    <- readRDS(here_input("sims-sensitivity-minimal.rds"))
# skills <- readRDS(here_output("synthesis", "skill.rds"))


#########################
#########################
#### Data preparation

selected_pars <- c("shape", "scale", "mobility",
                   "alpha", "beta", "gamma")

if (FALSE) {

  #### Aims
  # Here, we prepare a data.table of simulations:
  # * One row for each simulation
  # * A column distinguishing which parameter is changing (alpha, beta, etc.)
  # * A column with the degree of mis-specification
  # * Columns for the Mean Error
  # * We can use this to plot the error ~ degree of mis-specification for each parameter

  #### Define alg_pars
  # We will use alg_pars to define the simulations that are relevant for each parameter
  # This is necessary b/c we did not build the sims data.table with this information included
  # First, we read in alg_pars and then write an Excel file
  # This is modified:
  # * Add a parameter column that defines which parameter is changing
  # * Delete the extra row for alpha
  # Then we read in alg_pars & use it below to collate the sims data.table
  # alg_pars <- readRDS(here_input("alg_pars.rds"))
  # openxlsx::write.xlsx(alg_pars[[1]], here_input("alg_pars_1.xlsx"))
  alg_pars <-
    "alg_pars_1.xlsx" |>
    here_input() |>
    openxlsx::read.xlsx() |>
    as.data.table()
  default <- alg_pars[parameter == "default", ]

  #### Define performance simulations
  # We define the performance simulations for which we have implemented sensitivity analyses
  # It will be convenient, for each parameter, to include this information in the data.table below
  sims_performance <-
    sims_performance |>
    lazy_dt() |>
    filter(system_type == 1L) |>
    filter(n_receiver %in% c(20, 40, 60, 80, 100)) |>
    as.data.table()

  #### Build the sims data.table
  sims <-
    lapply(selected_pars, function(parameter_name) {

      # Define alg_par indices for the relevant parameter
      # (This will enable us to define the simulations where the value of that parameter changed)
      # parameter_name <- "alpha"
      pars <- alg_pars |> filter(parameter %in% c("default", parameter_name))

      # Define true value
      # (This is used to define the degree of mis-specification of the parameter, below)
      truth <- default[[parameter_name]]

      # Define the relevant simulations for that parameter
      sims_for_par <-
        # Join performance & sensitivity simulations
        # (i.e., We include the correct value for each parameter as well as all the mis-specified values)
        rbind(sims_performance, sims_senstivity) |>
        # Focus on the relevant parameter using the alg_par index
        filter(alg_par %in% pars$alg_par) |>
        # Select relevant columns & conveniently arrange rows (by parameter value)
        select(combination, system_type, path_type, array_type, arrangement, n_receiver,
               array_realisation, path_realisation, alg_par,
               alpha, beta, gamma, mobility, shape, scale) |>
        arrange(combination, system_type, path_type, array_type, arrangement, n_receiver,
                array_realisation, path_realisation, !!sym(parameter_name)) |>
        # (optional) We drop system_type and path_type as they are unused
        select(-system_type, -path_type) |>
        mutate(
          # Label the parameter of interest (to faciliate plotting)
          parameter = as.character(parameter_name),
          # Define the degree of mis-specification
          degree = !!sym(parameter_name) / truth,
          # Flag performance simulations (i.e., correct parameter values)
          performance = if_else(degree == 1, TRUE, FALSE)) |>
        as.data.table()

      # Check for duplicate rows
      if (nrow(sims_for_par) != nrow(distinct(sims_for_par))) {
        warning("There are duplicate rows!")
      }

      # Return outputs
      sims_for_par

    }, .combine = rbindlist)

  #### Add error statistics (~2.75 mins with 10 cl)
  # For each simulation, we compute the mean error
  # We focus on mean error as the performance skill analysis
  # ... shows this metric is most effective at distinguishing skill
  # We can compute mean error against:
  # * UD for true path
  # * UD for best algorithm implementation (i.e., the performance simulation
  # (optional) TO DO
  # * The efficiency of this code can be improved by only reading true/performance UDs once
  # * But it is easy to read this way

  sims[, index := seq_row(sims)]
  nrow(sims)
  sims <- cl_lapply(split(sims, seq_row(sims)),
                    .cl = 10L,
                    .chunk = TRUE,
                    .fun = function(sim) {

                      # Read 'true' UD
                      # sim <- sims[419, ]
                      print(sim$index)
                      path   <- terra::rast(here_alg(sim, "path", "ud.tif"))

                      # Read 'best' algorithm UDs
                      # * alg_par = "1" is the default parameters
                      acpf_best   <- read_rast(here_alg(sim, "patter", "acpf", "1", "ud-s.tif"))
                      acdcpf_best <- read_rast(here_alg(sim, "patter", "acdcpf", "1", "ud-s.tif"))

                      # Read algorithm UDs
                      acpf_alg   <- read_rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-s.tif"))
                      acdcpf_alg <- read_rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-s.tif"))

                      # Define data.table for results
                      out <- rbind(sim, sim)
                      out[, algorithm := c("ACPF", "ACDCPF")]

                      # Calculate skill against 'true' UD
                      out[, error_path := c(skill_mew(path, acpf_alg), skill_mew(path, acdcpf_alg))]

                      # Calculate skill against 'performance' algorithm implementation
                      out[, error_alg := c(skill_mew(acpf_best, acpf_alg), skill_mew(acdcpf_best, acdcpf_alg))]

                      out

                    }, .combine = rbindlist)

  #### Add convergence (~16 s, 1 cl)
  sims <-
    cl_lapply(split(sims, seq_row(sims)), function(sim) {

      # Define algorithm (acpf, acdcpf)
      algorithm <- sim$algorithm

      # Define path to convergence.rds
      convergence.rds <-
        here_alg(sim, "patter", algorithm, sim$alg_par, "convergence.rds")

      # Define convergence
      # * We use success = FALSE by default:
      # * Some simulations may fail before convergence.rds file is written
      # * (E.g., due to failure to sample initial locations)
      sim[, convergence := FALSE]

      # Record success
      if (file.exists(convergence.rds)) {
        sim[, convergence := readRDS(convergence.rds)]
      }

      sim

    }, .combine = rbindlist)


  #### Save sims
  saveRDS(sims, here_output("synthesis", "skill-sensitivity.rds"))

} else {

  sims <- readRDS(here_output("synthesis", "skill-sensitivity.rds"))

}

sims[, arrangement := factor(arrangement, levels = c("random", "regular"))]
sims[, algorithm := factor(algorithm, levels = c("ACPF", "ACDCPF"))]


#########################
#########################
#### Visualise models

#### Aims
# Visualise true & mis-specified models

#### Plot set up
# Define 'true' (default) parameters
defaults <-
  sims_performance |>
  slice(1L) |>
  select(alpha, beta, gamma, mobility, shape, scale) |>
  as.data.table()
# Define graphical parameters (e.g., cols)
gp <- data.frame(
  multiplier = c(0.1, 0.5, 1.0, 1.5, 2.0),
  label = c("0.1", "0.5", "1.0", "1.5", "2.0"),
  col = c("red", "blue", "black", "orange", "purple"),
  lty = c(2, 3, 1, 4, 5)
  )
gp_sbt     <- gp[gp$multiplier != 1.0, ]
multiplier <- gp_sbt$multiplier
cols       <- gp_sbt$col
ltys       <- gp_sbt$lty
index      <- seq_len(length(multiplier))
# Define plot limits
xlim <- c(0, max(c(sims$gamma, sims$mobility)))
ylim <- c(-0.05, 1)
# Define x values at which to calculate density
x <- seq(0, xlim[2], length.out = 1000)
# Define pretty_base() function
pretty_base <- function(x, y) {
  pretty_plot(x, y,
              xlim = xlim, ylim = ylim,
              type = "l", lwd = 4)
}

#### Visualise models
png(here_fig("sensitivity", "parameters.png"),
    height = 6, width = 8, units = "in", res = 600)
pp <- par(mfrow = c(2, 3), oma = c(2, 2, 2, 2), mar = c(1, 1, 1, 1))

# Visualise shape
y <- dtruncgamma(x, defaults$shape, defaults$scale, defaults$mobility, TRUE)
pretty_base(x, y)
shapes <- defaults$shape * multiplier
sapply(index, function(i) {
  lines(x, dtruncgamma(x, shapes[i], defaults$scale, defaults$mobility, TRUE),
        col = cols[i], lty = ltys[i])
})

# Visualise scale
pretty_base(x, y)
scales <- defaults$scale * multiplier
sapply(index, function(i) {
  lines(x, dtruncgamma(x, defaults$shape, scales[i], defaults$mobility, TRUE),
        col = cols[i], lty = ltys[i])
})

# Visualise mobility
pretty_base(x, y)
mobilities <- defaults$mobility * multiplier
sapply(index, function(i) {
  lines(x, dtruncgamma(x, defaults$shape, defaults$scale, mobilities[i], TRUE),
        col = cols[i], lty = ltys[i])
})

# Visualise alpha
y <- ddetlogistic(x, defaults$alpha, defaults$beta, defaults$gamma)
pretty_base(x, y)
alphas <- defaults$alpha * multiplier
sapply(index, function(i) {
  lines(x, ddetlogistic(x, alphas[i], defaults$beta, defaults$gamma),
        col = cols[i], lty = ltys[i])
})

# Visualise beta
pretty_base(x, y)
betas <- defaults$beta * multiplier
sapply(index, function(i) {
  lines(x, ddetlogistic(x, defaults$alpha, betas[i], defaults$gamma),
        col = cols[i], lty = ltys[i])
})

# Visualise gamma
pretty_base(x, y)
gammas <- defaults$gamma * multiplier
sapply(index, function(i) {
  lines(x, ddetlogistic(x, defaults$alpha, defaults$beta, gammas[i]),
        col = cols[i], lty = ltys[i])
})

dev.off()

#### Plot legend
png(here_fig("sensitivity", "parameters-legend.png"),
    height = 3, width = 3, units = "in", res = 600)
plot(0, 0, type = "n")
legend("topleft",
       legend = gp$label,
       col = gp$col,
       lty = gp$lty,
       lwd = c(1, 1, 4, 1, 1))
dev.off()


#########################
#########################
#### Convergence

#### Aim:
# Plot proportion of failures ~ degree of mis-specification
# > This is a coarse way of looking at 'algorithm sensitivity'

#### Define a 'convergence' data.table
# This defines the proportion of algorithm failures.
# This is evaluated across simulated path realisations, for each:
# * (System)
# * Array design
# * Parameter
# * Degree of mis-specification)

# sims |>
#   filter(arrangement == "random" & n_receiver == 20L &
#            parameter == "gamma" & degree == 0.1 & algorithm == "ACDCPF") |>
#   as.data.table() |>
#   View()

convergence <-
  sims |>
  group_by(combination, array_type, arrangement, n_receiver, array_realisation,
           parameter, degree, algorithm) |>
  summarise(n = n(),
            failure = length(which(convergence == FALSE)) / n) |>
  mutate(algorithm = factor(algorithm, levels = c("ACPF", "ACDCPF")),
         arrangement = factor(arrangement, levels = c("random", "regular")),
         parameter = factor(parameter, levels = c("shape", "scale", "mobility", "alpha", "beta", "gamma")),
         group = interaction(algorithm, arrangement, parameter, sep = ", ", lex.order = TRUE),
         group = factor(group, levels = levels(group), labels = LETTERS[1:length(levels(group))])
         ) |>
  as.data.table()

table(convergence$n)

#### Visualise the failure rate ~ degree of mis-specification
png(here_fig("sensitivity", "convergence.png"),
    height = 8, width = 12, units = "in", res = 600)
set.seed(1)
ggplot(convergence) +
  geom_jitter(aes(degree, failure, colour = factor(n_receiver), alpha = 0.5),
              width = 0.1, height = 0) +
  # facet_wrap(~algorithm + arrangement + parameter, ncol = 6) +
  facet_wrap(~group, ncol = 6) +
  scale_x_continuous(breaks = gp$multiplier, labels = gp$label) +
  # scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(face = "bold"))
dev.off()

#### Example checks
sims |>
  filter(parameter == "beta") |>
  mutate(check = beta) |>
  group_by(combination, array_type, arrangement, n_receiver, array_realisation,
           parameter, check, algorithm) |>
  summarise(failure = length(which(convergence == FALSE)) / n()) |>
  as.data.table() |>
  ggplot(aes(check, failure, colour = factor(n_receiver))) +
  geom_point() +
  facet_wrap(~algorithm + arrangement)


#########################
#########################
#### Maps

#### Aims
# For each parameter:
# * Map the 'true' & 'best' UDs compared to those from mis-specified algorithms
# * (for selected arrays).
# This provides a visual assessment of how maps are
# ... affected by parameter mis-specification.

#### Compute map area
area <- terra::cellSize(spat)

#### Make maps(~ 9 s, 6 cl; 45 s, 1 cl)
# This returns a list, with one element for each parameter
# Each parameter element is a list, with one element for each array design
# Each design element is a list with the UD for each panel
# We use UDs below to calculate home range area
# This code also produces maps
# NB: .cl must be NULL if the UDs are returned
uds <-
  cl_lapply(selected_pars,
          .cl = NULL,
          .fun = function(parameter_name) {

  # parameter_name <- "alpha"
  png(here_fig("sensitivity", paste0("map-", parameter_name, ".png")),
      height = 8, width = 12, units = "in", res = 600)
  pp <- par(mfrow = c(4, 7),
            oma = c(0, 0, 0, 0), mar = rep(1, 4))

  # For the selected arrange design...
  range(sims$n_receiver)
  designs <- CJ(n_receiver = range(sims$n_receiver),
                arrangement = c("random", "regular"))

  uds_by_design <-
    lapply(split(designs, seq_row(designs)), function(design) {

    # design <- designs[1, ]
    sims_for_maps <-
      sims |>
      filter(parameter == parameter_name) |>
      filter(arrangement == design$arrangement) |>
      filter(n_receiver %in% design$n_receiver) |>
      filter(array_realisation == 1L & path_realisation == 1L) |>
      arrange(n_receiver, arrangement) |>
      as.data.table()

    #### Read UDs
    # Baseline simulation parameters
    sim             <- sims_for_maps[1, ]
    # Selected folders
    # unique(sims$degree)
    # Use 0.5 and 2.0 degree of mis-specification as
    # ... for movement parameters simulations failed with more severe underestimation
    folder_under    <- sims_for_maps$alg_par[sims_for_maps$degree == 0.5][1]
    folder_best     <- "1"
    folder_over     <- sims_for_maps$alg_par[sims_for_maps$degree == 2.0][1]
    # UDs
    ud_path         <- read_rast(here_alg(sim, "path", "ud.tif"))
    ud_acpf_under   <- read_rast(here_alg(sim, "patter", "acpf", folder_under, "ud-s.tif"))
    ud_acpf_best    <- read_rast(here_alg(sim, "patter", "acpf", folder_best, "ud-s.tif"))
    ud_acpf_over    <- read_rast(here_alg(sim, "patter", "acpf", folder_over, "ud-s.tif"))
    ud_acdcpf_under <- read_rast(here_alg(sim, "patter", "acdcpf", folder_under, "ud-s.tif"))
    ud_acdcpf_best  <- read_rast(here_alg(sim, "patter", "acdcpf", folder_best, "ud-s.tif"))
    ud_acdcpf_over  <- read_rast(here_alg(sim, "patter", "acdcpf", folder_over, "ud-s.tif"))

    #### Define scaling parameter
    uds <- list(path = ud_path,
                acpf_under = ud_acpf_under,
                acpf_best = ud_acpf_best,
                acpf_over = ud_acpf_over,
                acdcpf_under = ud_acdcpf_under,
                acdcpf_best = ud_acdcpf_best,
                acdcpf_over = ud_acdcpf_over)
    scale <- sapply(uds, \(ud) {
      if (is.null(ud)) {
        return(NA)
      } else {
        terra::global(ud, "max")$max
      }
    }) |> max(na.rm = TRUE)
    stopifnot(!is.na(scale))

    #### Plot UDs for selected array
    m <- read_array(sim)
    lapply(seq_len(length(uds)), function(ud_index) {
      # If the UD is NULL, make a blank placeholder plot
      # ud_index
      ud <- uds[[ud_index]]
      if (is.null(ud)) {
        terra::plot(dat_gebco(), axes = FALSE, legend = FALSE, col = NA)
        return(NULL)
      }
      # Otherwise, plot the UD
      sud   <- ud / scale
      range <- c(0, 1)
      spatMap(sud, range = range, legend = FALSE, mar = NA)
      # (optional) Set speed to TRUE to check plot layout only
      speed <- FALSE
      if (speed) {
        spatAxes(ud)
        return(NULL)
      }
      # Add home range (based on unscaled UD)
      hr <- map_hr_home(ud, .add = TRUE,
                        border = "dimgrey", lwd = 0.75)
      # Plot UD (scaled)
      # terra::plot(ud)
      # Add receivers
      points(m$receiver_x, m$receiver_y,
             pch = 21,
             bg = scales::alpha("black", 0.75),
             col = scales::alpha("black", 0.75),
             cex = 0.5)
      spatAxes(ud)
      # terra::sbar(sim$gamma, lonlat = FALSE, col = "darkred")
      # Return raw ud
      ud
    })

  })

  par(pp)
  dev.off()
  uds_by_design

}) |> invisible()

#### Compute home-range area
# Build a data.table of home range areas
hr <- CJ(parameter = 1:6, design = 1:4, implementation = 1:7, area = NA_real_)
hr <- cl_lapply(split(hr, seq_row(hr)), function(d) {
  # Get UD
  ud <- uds[[d$parameter]][[d$design]][[d$implementation]]
  # Compute home-range area (km2)
  if (is.null(ud)) {
    return(d)
  }
  hr <- map_hr_home(ud, .add = FALSE)
  hr    <- terra::classify(hr, cbind(0, NA))
  hr_cz <- terra::mask(area, hr)
  # terra::plot(hr_cz)
  d[, area := terra::global(hr_cz, "sum", na.rm = TRUE)[1, 1] / 1e6]
  d
}, .combine = rbindlist)
# Define convenience variables/objects for plotting
hr[, group := .GRP, by = .(parameter, design)]
ylim <- c(0, max(hr$area, na.rm = TRUE))

#### Visualise home range area
# Columns: Parameters
# Rows: Array designs
# We will manually add the plots to the right-hand side of the relevant map

# Check column fill order
pp <- par(mfcol = c(4, 6))
for (i in 1:24) {
  plot(i, main = i)
}
par(pp)

# Make figure
png(here_fig("sensitivity", "hr-area.png"),
    height = 6, width = 10, units = "in", res = 600)
pp <- par(mfcol = c(4, 6),
          oma = c(2, 2, 2, 2),
          mar = c(1, 1, 1, 1))
lapply(split(hr, hr$group), function(d) {
  # Process data.table for convenience
  # d <- split(hr, hr$group)[[12]]
  d$x <- d$implementation
  d$y <- d$area
  # Define bar colours (distinguish algorithms + implementations)
  # path, acpf_under, acpf_best, acpf_over, acdcpf_under, acdcpf_best, acdcpf_over
  d$col <- c("dimgrey",
             alpha("red", 0.2), alpha("red", 0.5), alpha("red", 0.8),
             alpha("blue", 0.2), alpha("blue", 0.5), alpha("blue", 0.8))
  # Distinguish 'best'
  d$lwd <- c(1,
             1, 2, 1,
             1, 2, 1)
  # Create blank plot
  axis_ls <- pretty_plot(d$x,
                         d$y,
                         xlim = c(0.5, max(d$x) + 0.5),
                         ylim = ylim,
                         pretty_axis = list(axis = list(list(side = 1, at = d$x),
                                                        list(NULL))),
                         xlab = "", ylab = "",
                         type = "n")
  # Add bars
  for (j in seq_len(nrow(d))) {
    rect(xleft = d$x[j] - 0.5,
         ybottom = axis_ls[[2]]$lim[1],
         xright = d$x[j] + 0.5,
         ytop = d$y[j],
         col = d$col[j],
         lwd = d$lwd[j])

  }
}) |> invisible()
par(pp)
dev.off()


#########################
#########################
#### Boxplots

#### Choose error metric
# Use `error_path` or `error_alg`
# metric <- "algorithm" # (deprecated)
metric <- "path"
if (metric == "path") {
  sims[, error := error_path]
} else if (metric == "algorithm") {
  sims[, error := error_alg]
} else {
  stop("Metric not recognised!")
}
# Express error metrics relative to the correct implementation for a given path realisation
# * If metric = "algorithm", everything is expressed relative to zero error
# * If metric = "path", everything is relative to the best error
# sims |>
#    filter(combination == 1L & array_type == 2L &
#             arrangement == "random" & n_receiver == 20L &
#             array_realisation == 1L & path_realisation == 1L &
#             parameter == "gamma" & algorithm == "ACDCPF") |>
#    as.data.table() |>
#    View()
boxsims <- copy(sims)
if (metric == "path") {
  boxsims <-
    boxsims |>
    as.data.frame() |>
    group_by(combination, array_type,
             arrangement, n_receiver,
             array_realisation, path_realisation,
             parameter, algorithm) |>
    mutate(error = error / error[performance == TRUE]) |>
    ungroup() |>
    mutate(arrangement = factor(arrangement, c("random", "regular"), c("Random", "Regular")),
           n_receiver = factor(n_receiver, unique(sort(n_receiver)), paste("n =", unique(sort(n_receiver)))),
           group = interaction(arrangement, n_receiver,sep = ", ", lex.order = FALSE)) |>
           # group = factor(group, levels = levels(group), labels = 1:length(levels(group)))) |>
    filter(!is.na(error)) |>
    as.data.table()
}

#### Create boxplots
# ylim <- range(boxsims$error)
ylim <- c(NA, NA)
cl_lapply(selected_pars, .cl = length(selected_pars), .fun = function(parameter_name) {

  # parameter_name <- "alpha"
  print(parameter_name)
  png(here_fig("sensitivity", paste0("boxplot-", metric, "-", parameter_name, ".png")),
      height = 10, width = 6, units = "in", res = 600)

  p <-
    boxsims |>
    filter(parameter == parameter_name & degree != 1.0) |>
    as.data.frame() |>
    ggplot() +
    geom_boxplot(aes(x = factor(degree), y = error, fill = algorithm),
                 varwidth = TRUE) +
    scale_x_discrete(breaks = gp$multiplier, labels = gp$label) +
    # ylim(ylim) +
    # scale_y_continuous(labels = sci_notation) +
    # facet_wrap(~n_receiver + arrangement, ncol = 2) +
    xlab("") + ylab("") +
    facet_wrap(~group, ncol = 2) +
    theme_bw() +
    theme(strip.text = element_text(face = "bold")) +
    theme(panel.grid.minor = element_blank())

  print(p)

  dev.off()

})


#### End of code.
#########################
#########################

