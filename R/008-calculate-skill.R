#########################
#########################
#### calculate-skill.R

#### Aims
# 1) Calculates skill metrics from simulation outputs

#### Prerequisites
# 1) Analyse simulated data using selected algorithms

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
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
sapply(list.files(here::here("src"), full.names = TRUE), source)


#########################
#########################
#### Calculate skill

#### Collate skill scores
# * Loop over sims
# * Read path UD & algorithm UD (if available)
# * Calculate skill for relevant algorithms

skills <- pbapply::pblapply(split(sims, seq_len(nrow(sims))), function(sim) {

  #### Read UDs
  path   <- terra::rast(here_alg(sim, "path", "ud.tif"))
  #
  # TO DO (add control, a blank)
  #
  if (sim$performance) {
    coa_1  <- read_rast(here_alg(sim, "coa", "30 mins", "ud.tif"))
    coa_2  <- read_rast(here_alg(sim, "coa", "120 mins", "ud.tif"))
  }
  acpf   <- read_rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud.tif"))
  acdcpf <- read_rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud.tif"))

  #### Define blank skill table
  skill <- data.table(
    id = sim$id,
    performance = sim$performance,
    alg = c("COA(30)", "COA(120)", "ACPF", "ACDCPF"),
    mb = NA_real_,
    me = NA_real_,
    rmse = NA_real_,
    R = NA_real_,
    d = NA_real_
  )

  #### Calculate skill metrics
  if (sim$performance) {
    # We will calculate model skill scores for each algorithm
    uds <- list(coa_1, coa_2, acpf, acdcpf)
  } else {
    # We will calculate model skill scores for patter algorithms only
    uds <- list(NULL, NULL, acpf, acdcpf)
  }
  skill$mb   <- skill_by_alg(uds, path, skill_mb)
  skill$me   <- skill_by_alg(uds, path, skill_me)
  skill$rmse <- skill_by_alg(uds, path, skill_rmse)
  skill$R    <- skill_by_alg(uds, path, skill_R)
  skill$d    <- skill_by_alg(uds, path, skill_d)

  #### Return outputs
  saveRDS(skill, here_output("skill", paste0(sim$id, ".rds")))
  skill

}) |> rbindlist()

saveRDS(skills, here_data("sims", "synthesis", "skill-raw.rds"))


#########################
#########################
#### Process skill scores

#### Drop simulations which failed for all algorithms
# skills <- skills[rowSums(!is.na(skills[, c("mb", "me", "rmse", "R", "d")])) > 0, ]

#### Add required information
# merge with sims()

#### Save processed data
saveRDS(skills, here_data("sims", "synthesis", "skill.rds"))


#### End of code.
#########################
#########################
