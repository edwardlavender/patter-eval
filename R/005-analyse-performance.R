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
library(tictoc)
sapply(list.files(here::here("src"), full.names = TRUE), source)

#### Load data
sims <- readRDS(here_input("sims.rds"))


#########################
#########################
#### Data preparation

#### Select results for a specific combination
# Filter sims
sims <-
  sims |>
  filter(performance) |>
  as.data.table()
# Define function that assigns png names
png_name <- function(name) {
  paste0(name, "-",  sims$combination[1], ".png")
}

#### Define selected receiver numbers
# On some plots, we only consider two different numbers of receivers
nr <- c(10, 100)
metrics <- c("mb", "me", "rmse", "R", "d")


#########################
#########################
#### Maps of space use

#### Select simulations
# Use specific combination, two receiver levels, one array/path realisation
sims_for_maps <-
  sims |>
  filter(combination == 1L) |>
  filter(n_receivers == nr) |>
  filter(array_realisation == 1L & path_realisation == 1L) |>
  arrange(arrangement, n_receivers)
  as.data.table()

#### Set up plot
# The figure will be annotated outside of R.
png(here_fig("performance", png_name("map")),
    height = 10, width = 10, units = "in", res = 600)
pp <- par(mfrow = c(4, 7))
lapply(1:4, function(i) {

  # Read UDs
  id       <- sims_for_maps$sim[i]
  ud_path  <- terra::rast(here_alg(id, "path", "ud.tif"))
  ud_coa_1 <- terra::rast(here_alg(id, "coa", "30 mins", "ud.tif"))
  ud_coa_2 <- terra::rast(here_alg(id, "coa", "30 mins", "ud.tif"))
  # ud_rsp_1
  # ud_rsp_2
  ud_acpf   <- terra::rast(here_alg(id, "patter", "acpf", "ud.tif"))
  ud_acdcpf <- terra::rast(here_alg(id, "patter", "acdcpf", "ud.tif"))
  uds <- list(ud_coa_1, ud_coa_2,
              ud_acpf, ud_acdcpf)

  # Scale UDs (within rows) to a maximum value of 1 for comparison
  scale <- sapply(uds, \(ud) terra::global(ud, "max")$max) |> max()
  uds   <- lapply(uds, \(ud) ud / scale)

  # Plot UDs for selected array
  lapply(uds, terra::plot) |> invisible()

}) |> invisible()
dev.off()


#########################
#########################
#### Barplots of error statistics

#### Select simulations
# Select specific path realisation
sims <-
  sims |>
  filter(combination == 1L) |>
  filter(n_receivers == nr) |>
  filter(array_realisation == 1L) |>
  arrange(arrangement, n_receivers) |>
  as.data.table()

#### Visualise boxplots
png(here_fig("performance", png_name("barplots")),
    height = 10, width = 10, units = "in", res = 600)
pp <- par(mfrow = c(4, 6))
lapply(1:4, function(i) {
  # Load skill scores for selected simulation run
  id <- sims$id[i]
  skill <- readRDS(here_output("skill", paste0(id, ".rds")))
  # Visualise distribution of error scores (by metric)
  lapply(metrics, function(metric) {
    pretty_boxplot(skill$alg, skill[[metric]])
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
