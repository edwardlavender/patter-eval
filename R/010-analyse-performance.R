#########################
#########################
#### analyse-performance.R

#### Aims:
# (1) Compares algorithm performance

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
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(prettyGraphics)
library(tictoc)
dv::src()

#### Load data
# Performance simulations
sims   <- readRDS(here_input("sims-performance.rds"))
skills <- readRDS(here_data("sims", "synthesis", "skill.rds"))


#########################
#########################
#### Data preparation

# Define function that assigns png names
png_name <- function(name) {
  paste0(name, "-",  sims$combination[1], ".png")
}

#### Define selected receiver numbers
# On some plots, we only consider two different numbers of receivers
nr <- c(10, 100)
metrics <- c("mb", "me", "rmse", "R", "d")
metrics_lim <-
  list(mb = NULL,
       me = NULL,
       rmse = NULL,
       R = c(0, 1),
       d = c(0, 1))

#### Define algorithms
# Define algorithms (subset)
algs <- data.frame(alg = c("Null", "COA(30)", "COA(120)", "RSP(1)", "RSP(2)", "ACPF(F)", "ACDCPF(F)", "ACPF(S)", "ACDCPF(S)"),
                   col_name = c("dimgrey", "red", "darkred", "orange", "darkorange", "slateblue1", "slateblue4", "lightgreen", "darkgreen")
                   )
algs$col <- scales::alpha(algs$col_name, 0.5)
algs$label <- seq_len(nrow(algs))
skills <-
  skills |>
  filter(alg %in% algs$alg) |>
  left_join(algs, by = "alg", relationship = "many-to-one") |>
  as.data.table()


#########################
#########################
#### Maps of space use

# Make maps (~37 s)
if (FALSE) {

  #### Select simulations
  # Use specific combination, two receiver levels, one array/path realisation
  sims_for_maps <-
    sims |>
    filter(combination == 1L) |>
    filter(n_receiver %in% nr) |>
    filter(array_realisation == 1L & path_realisation == 1L) |>
    arrange(arrangement, n_receiver) |>
    as.data.table()

  #### Set up plot
  # The figure will be annotated outside of R
  # (optional) FLAG: adjust width & number of columns for algorithms
  dir.create(here_fig("performance"), recursive = TRUE)
  png(here_fig("performance", png_name("map")),
      height = 8, width = 12, units = "in", res = 600)
  pp <- par(mfrow = c(4, 9),
            oma = c(0, 0, 0, 0), mar = rep(1, 4))
  pbapply::pblapply(1:4, function(i) {

    #### Read UDs
    sim        <- sims_for_maps[i, ]
    # list.files(here_alg(sim), recursive = TRUE)
    ud_path    <- terra::rast(here_alg(sim, "path", "ud.tif"))
    ud_coa_1   <- terra::rast(here_alg(sim, "coa", "30 mins", "ud.tif"))
    ud_coa_2   <- terra::rast(here_alg(sim, "coa", "120 mins", "ud.tif"))
    ud_rsp_1   <- terra::rast(here_alg(sim, "rsp", "default", "ud.tif"))
    ud_rsp_2   <- terra::rast(here_alg(sim, "rsp", "custom", "ud.tif"))
    ud_acpff   <- terra::rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-f.tif"))
    ud_acdcpff <- terra::rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-f.tif"))
    ud_acpfk   <- terra::rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-k.tif"))
    ud_acdcpfk <- terra::rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-k.tif"))
    ud_acpfs   <- terra::rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-s.tif"))
    ud_acdcpfs <- terra::rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-s.tif"))
    # Select UDs
    # * (optional) FLAG: modify UDs included in this list
    uds <- list(ud_path,
                ud_coa_1, ud_coa_2,
                ud_rsp_1, ud_rsp_2,
                ud_acpff, ud_acdcpff,
                # ud_acpfk, ud_acdcpfk,
                ud_acpfs, ud_acdcpfs
    )

    #### Calculate scaling parameter
    # We will scale UDs (within rows) to a maximum value of 1 for comparison
    # Scaling is implemented below to ensure correct HR calculation
    scale <- sapply(uds, \(ud) terra::global(ud, "max")$max) |> max()
    # uds   <- lapply(uds, \(ud) ud / scale)
    # sapply(uds, \(ud) terra::global(ud, "max"))

    #### Plot UDs for selected array
    m <- read_array(sim)
    lapply(uds, function(ud) {
      # Plot UD (scaled)
      sud <- ud / scale
      spatMap(sud, range = c(0, 1), legend = FALSE, mar = NA)
      # (optional) Set speed to TRUE to check plot layout only
      speed <- FALSE
      if (speed) {
        spatAxes(ud)
        return(NULL)
      }
      # Add home range (based on unscaled UD)
      map_hr_home(ud, .add = TRUE,
                  border = "dimgrey", lwd = 0.75)
      # Plot UD (scaled)
      # terra::plot(ud)
      # Add receivers
      points(m$receiver_easting, m$receiver_northing,
             pch = 21,
             bg = scales::alpha("black", 0.75),
             col = scales::alpha("black", 0.75),
             cex = 0.5)
      # Add detection range(s)
      # * We use a scale bar instead
      if (FALSE) {
        cbind(m$receiver_easting, m$receiver_northing) |>
          terra::vect() |>
          terra::buffer(width = sim$gamma) |>
          terra::lines(lwd = 0.5, col = "dimgrey")
      }
      spatAxes(ud)
      terra::sbar(sim$gamma, lonlat = FALSE, col = "darkred")
      NULL
    }) |> invisible()

  }) |> invisible()
  dev.off()

}




#########################
#########################
#### Barplots of error statistics

#### Select simulations
# Define arrays
combs <- CJ(arrangement = c("random", "regular"),
            n_receiver = nr)
# Define simulations
sims_for_skill <-
  sims |>
  filter(combination == 1L) |>
  filter(n_receiver %in% nr) |>
  arrange(arrangement, n_receiver) |>
  as.data.table()

#### Visualise boxplots
png(here_fig("performance", png_name("barplots")),
    height = 10, width = 15, units = "in", res = 600)
pp <- par(mfrow = c(4, length(metrics)),
          oma = c(1, 3, 1, 1), mar = c(2, 2, 2, 2))
lapply(1:4, function(i) {
  # Define skill scores across path realisations
  skill <- skills[id %in%
                  sims_for_skill[arrangement == combs$arrangement[i] & n_receiver == combs$n_receiver[i], ]$id
                  ]

  # Visualise distribution of error scores (by metric)
  lapply(metrics, function(metric) {
    # metric <- metrics[1]
    x <- skill$label # skill$alg
    y <- skill[[metric]]
    # p <- min(floor(log10(abs(y))), na.rm = TRUE)
    # y <- y / 10^p
    pretty_boxplot(x, y,
                   varwidth = TRUE,
                   col = skill$col,
                   ylim = metrics_lim[[metric]],
                   pretty_axis_args = list(side = 1:2, pretty = list(list(n = 100), list(n = 5))),
                   las = TRUE,
                   xlab = "", ylab = "")
  }) |> invisible()
}) |> invisible()
dev.off()


#########################
#########################
#### Trends

# TO DO
# * use n_receivers or detection coverage here

#### (optional) Subset algorithms for improved clarity on plot
unique(skills$alg)
skills_for_trends <-
  skills |>
  filter(alg %in% c("Null", "COA(30)", "RSP(1)", "ACPF(S)", "ACDCPF(S)")) |>
  as.data.table()
combs <- unique(sims$combination)
nc    <- length(unique(sims$combination))


png(here_fig("performance", "relationships.png"),
    height = 4 * nc, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(nc * 2, length(metrics)), oma = c(1, 3, 1, 1), mar = rep(2, 4))
lapply(combs, function(combination) {

  # Identify skill datasets
  # combination <- 1
  skill <- skills_for_trends[system_type == combination, ]

  # Loop over array arrangements
  lapply(split(skill, skill$arrangement), function(sk) {
    # Visualise skill ~ n_receivers
    lapply(metrics, function(metric) {
      # Create plot, distinguishing between algorithms
      # metric <- metrics[1]
      x <- sk$n_receiver + rnorm(nrow(sk), 0, 1)
      y <- sk[[metric]]
      pretty_plot(x, y,
                  col = scales::alpha(sk$col, 0.25),
                  xlab = "", ylab = "",
                  cex = 0.5, lwd = 0.5)
      # Add trends for each algorithm
      lapply(split(sk, sk$alg, drop = TRUE), function(d) {
       # d <- split(sk, sk$alg)[[1]]
        add_smoother(d$n_receiver, d[[metric]], col = d$col, lwd = 2)
      }) |> invisible()
    }) |> invisible()

  }) |> invisible()
}) |> invisible()
dev.off()

# Check colours
algs



#### End of code.
#########################
#########################
