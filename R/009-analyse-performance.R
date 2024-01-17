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

#### Skill factors
unique(skills$alg)
skills$alg <- factor(skills$alg,
                     c("Null", "COA(30)", "COA(120)", "ACPF(K)", "ACDCPF(K)"),
                     labels = 0:4)


#########################
#########################
#### Maps of space use

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
# The figure will be annotated outside of R.
png(here_fig("performance", png_name("map")),
    height = 7, width = 10, units = "in", res = 600)
pp <- par(mfrow = c(4, 5))
pbapply::pblapply(1:4, function(i) {

  #### Read UDs
  sim      <- sims_for_maps[i, ]
  # list.files(here_alg(sim), recursive = TRUE)
  ud_path  <- terra::rast(here_alg(sim, "path", "ud.tif"))
  ud_coa_1 <- terra::rast(here_alg(sim, "coa", "30 mins", "ud.tif"))
  ud_coa_2 <- terra::rast(here_alg(sim, "coa", "120 mins", "ud.tif"))
  # ud_rsp_1
  # ud_rsp_2
  ud_acpf   <- terra::rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-k.tif"))
  ud_acdcpf <- terra::rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-k.tif"))
  uds <- list(ud_path,
              ud_coa_1, ud_coa_2,
              ud_acpf, ud_acdcpf)

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
    spatMap(sud, range = c(0, 1), legend = FALSE)
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
pp <- par(mfrow = c(4, 5), oma = c(1, 3, 1, 1), mar = c(2, 2, 2, 2))
lapply(1:4, function(i) {
  # Define skill scores across path realisations
  skill <- skills[id %in%
                  sims_for_skill[arrangement == combs$arrangement[i] & n_receiver == combs$n_receiver[i], ]$id
                  ]

  # Visualise distribution of error scores (by metric)
  lapply(metrics, function(metric) {
    # metric <- metrics[1]
    x <- skill$alg
    y <- skill[[metric]]
    # p <- min(floor(log10(abs(y))), na.rm = TRUE)
    # y <- y / 10^p
    pretty_boxplot(x, y,
                   varwidth = TRUE,
                   ylim = metrics_lim[[metric]],
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
#

png(here_fig("performance", "relationships.png"),
    height = 10, width = 10, units = "in", res = 600)
pp <- par(mfrow = c(4, 5))
lapply(split(sims, sims$combination), function(sim) {

  # Read all skill datasets
  skill <-
    here_output("skill", paste0(sim$id, ".rds")) |>
    lapply(readRDS) |>
    rbindlist()
  # TO DO
  # Add array details (number, arrangement) by merging with sims

  # Loop over array arrangements
  lapply(split(skill, skill$arrangement), function(sk) {
    # Summarise average skill across all simulations
    # * For coverage, coverage statistics may need to be binned
    sk_sry <-
      sk |>
      group_by(alg, n_receivers) |>
      mutate(across(all_of(metrics)), mean) |>
      slice(1L) |>
      as.data.table()

    # Visualise skill ~ n_receivers
    lapply(metrics, function(metric) {
      # (optional) Assign colours based on skill quantile bins at each n_receiver level
      # utils.add::find_quantile_bin()
      # Create plot, distinguishing between algorithms
      pretty_plot(sk$n_receivers, sk[[metric]],
                  colour = sk$alg)
      # Add average line across all simulations
      lapply(split(sk_sry, sk_sry$alg), function(d) {
        lines(d$n_receivers, d[[metric]], col = d$alg, lwd = 2)
      }) |> invisible()
    }) |> invisible()

  }) |> invisible()
}) |> invisible()
dev.off()


#### End of code.
#########################
#########################
