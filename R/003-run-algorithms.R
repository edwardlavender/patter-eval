#########################
#########################
#### run-algorithms.R

#### Aims
# 1) Implement algorithms using simulated data

#### Prerequisites
# 1) Simulate data & prepare for algorithm implementation

#### Time records
# *


#########################
#########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()
options(error = function(...) beepr::beep(7))

#### Essential packages
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tictoc)
sapply(list.files(here::here("src"), full.names = TRUE), source)

#### Load data
grid       <- terra::rast(here_input("grid.tif"))
sims       <- readRDS(here_input("sims.rds"))


#########################
#########################
#### Generic preparation

#### (optional) Set up cluster
grid <- terra::wrap(grid)


#########################
#########################
#### Performance simulations

#### Overview
# Here, the workflows required for 'performance' simulations are implemented:
# * Generate true UDs (for the path)
# * Generate COA UDs
# * Generate RSP UDs (TO DO)

#### Prepare data
# Isolate relevant simulations
sims_for_performance <-
  sims |>
  filter(performance) |>
  as.data.table()
# Define list to loop over
sims_for_performance_ls <-
  split(sims_for_performance, seq_len(nrow(sims_for_performance)))

#### Path UDs (~3 mins, one core)
pbapply::pblapply(sims_for_performance_ls, function(sim) {
  # sim <- sims_for_performance_ls[[1]]
  grid <- terra::unwrap(grid)
  workflow_path(sim, grid)
}) |> invisible()

#### COA UDs (~3 mins, one core)
pbapply::pblapply(sims_for_performance_ls, function(sim) {
  # sim <- sims_for_performance_ls[[1]]
  grid <- terra::unwrap(grid)
  workflow_coa(sim, grid)
}) |> invisible()

#### RSP UDs
# TO DO

#### Quick checks
s <- sims_for_performance_ls[[2]]
here_alg(s, "path", "ud.tif")  |> terra_qplot()
here_alg(s, "coa", "120 mins", "ud.tif") |> terra_qplot()


#########################
#########################
#### Patter simulations

#### Overview
# Here, {patter} workflows are implemented (for both performance & sensitivity simulations)

#### Time trials
# Estimated duration (days) on {cl} CPUs, assuming simulations take {guess} seconds
cl    <- 10L
guess <- 30 # 30 s
(nrow(sims) * 30)/60/60/24/cl

pbapply::pblapply(split(sims, seq_len(nrow(sims))), function(sim) {
  # sim <- sims[1, ]
  tic()
  grid <- terra::unwrap(grid)
  workflow_patter(sim, grid)
  toc()
}) |> invisible()



#### End of code.
#########################
#########################
