#########################
#########################
#### prep-sims.R

#### Aims:
# (1) Prepare algorithm inputs for simulations
# * This code is implemented prior to simulations for improved speed

#### Prerequisites
# 1) Simulate data (sim-data.R)


#########################
#########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)

#### Essential packages
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tictoc)

#### Load data
source(here_r("002-define-helpers.R"))
grid     <- terra::rast(here_input("grid.tif"))
arrays   <- readRDS(here_input("arrays.rds"))
sims     <- readRDS(here_input("sims.rds"))


#########################
#########################
#### Prepare AC* inputs

#### Isolate detection parameters
det_pars_all <-
  sims |>
  select(array_type, array_realisation, alpha, beta, gamma) |>
  distinct() |>
  as.data.table()

#### Define directories
# Clean directories
if (FALSE) {
  tic()
  unlink(here_input("ac"), recursive = TRUE)
  toc()
}
# Define folders
# * ac/{array_type}/{array_realisation}/{gamma}/overlaps.rds
# * ac/{array_type}/{array_realisation}/{gamma}/{alpha}/{beta}/kernels.rds
pbapply::pbsapply(split(det_pars_all, seq_len(nrow(det_pars_all))), function(d) {
  dir.create(here_input("ac", d$array_type, d$array_realisation,
                        d$gamma, d$alpha, d$beta),
             recursive = TRUE)
}) |> invisible()

#### Define detection overlaps (2 hr 10 using 10 cl)
# Define gamma parameters (for each array)
gc()
gammas <-
  det_pars_all |>
  select(array_type, array_realisation, gamma) |>
  distinct() |>
  as.data.table()
tic()
pbapply::pblapply(split(gammas, seq_len(nrow(gammas))), cl = 10L, function(d) {
  # Isolate array data
  # d = split(gammas, seq_len(nrow(gammas)))[[1]]
  out_file <- here_input("ac", d$array_type, d$array_realisation, d$gamma, "overlaps.rds")
  if (!file.exists(out_file)) {
    array <-
      arrays[[d$array_type]] |>
      filter(array_id == d$array_realisation) |>
      mutate(receiver_range = d$gamma) |>
      as.data.table()
    # Define detection containers
    containers <- acs_setup_detection_containers(grid, array)
    # Define overlaps & save
    overlaps   <- acs_setup_detection_overlaps(containers, array)
    saveRDS(overlaps, out_file)
    TRUE
  }
}) |> invisible()
toc()

#### Define detection containers (~45 mins [?] with 10 cl)
gc()
tic()
pbapply::pblapply(split(det_pars_all, seq_len(nrow(det_pars_all))), cl = 10L, function(d) {
  # Define array data
  out_file <- here_input("ac", d$array_type, d$array_realisation,
                         d$gamma, d$alpha, d$beta, "kernels.rds")
  if (!file.exists(out_file)) {
    array <-
      arrays[[d$array_type]] |>
      filter(array_id == d$array_realisation) |>
      mutate(receiver_range = d$gamma) |>
      as.data.table()
    # Read overlaps
    overlaps <- readRDS(here_input("ac", d$array_type, d$array_realisation,
                                   d$gamma, "overlaps.rds"))
    # Define detection containers & save
    kernels <- acs_setup_detection_kernels(array,
                                           .bathy = grid,
                                           .calc_detection_pr = acs_setup_detection_pr,
                                           .alpha = d$alpha, .beta = d$beta, .gamma = d$gamma,
                                           .verbose = FALSE)
    kernels <- wrapr(kernels)
    saveRDS(kernels, out_file)
  }
  TRUE
}) |> invisible()
toc()


#### End of code.
#########################
#########################
