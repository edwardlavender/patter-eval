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
alg_pars <- readRDS(here_input("alg_pars.rds"))


#########################
#########################
#### Prepare AC* inputs

#### Clean directories
if (FALSE) {
  tic()
  unlink(here_input("ac"), recursive = TRUE)
  toc()
}

#### Define detection parameters
# Collect all unique parameter combinations (alpha, beta, gamma)
det_pars_all <-
  alg_pars |>
  # purrr::flatten() |>
  rbindlist() |>
  dplyr::select(alpha, beta, gamma) |>
  distinct() |>
  as.data.table()
# Isolate unique gamma values
gammas <- sort(unique(det_pars_all$gamma))
sapply(gammas, \(.gamma) dir.create(here_input("ac", .gamma)))

#### Build detection containers, overlaps and containers (~7 mins)
# This code writes AC inputs to file:
# * ac/{gamma}/overlaps.rds
# * ac/{gamma}/{alpha}/{beta}/kernels.rds (required unwrapping)
dir.create(here_input("ac"))
tic()
ac_pars <-
  lapply(arrays, function(.arrays) {
    lapply(split(.arrays, .arrays$array_id), function(.array) {
      objs <-
        lapply(gammas, function(.gamma) {

          ##### Define array & gamma-specific objects (containers, overlaps)
          .array$receiver_range <- .gamma
          containers <- acs_setup_detection_containers(grid, .array)
          overlaps   <- acs_setup_detection_overlaps(containers, .array)
          saveRDS(overlaps, here_input("ac", .gamma, "overlaps.rds"))

          #### Define array and alpha, beta, gamma-specific objects (kernels)
          # Define unique alpha/beta/gamma combinations
          ab <-
            det_pars_all |>
            filter(gamma == .gamma) |>
            distinct() |>
            as.data.table()
          # Define a list of kernels, by alpha/beta (and gamma) combination
          kernels <-
            lapply(split(ab, seq_len(nrow(ab))), function(.ab) {
              dir.create(here_input("ac", .gamma, .ab$alpha, .ab$beta),
                         recursive = TRUE)
              k <- acs_setup_detection_kernels(.array,
                                               .bathy = grid,
                                               .calc_detection_pr = acs_setup_detection_pr,
                                               .alpha = .ab$alpha, .beta = .ab$beta, .gamma = .ab$gamma,
                                               .verbose = FALSE)
              k <- wrapr(k)
              saveRDS(k, here_input("ac", .gamma, .ab$alpha, .ab$beta, "kernels.rds"))
              NULL
            })
          names(kernels) <- paste(ab$alpha, ab$beta)

          #### Return outputs
          # list(overlaps = overlaps, kernels = kernels)
          NULL

        })
      names(objs) <- gammas
      objs
    })
  })
toc()


#### End of code.
#########################
#########################
