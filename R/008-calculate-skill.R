#########################
#########################
#### calculate-skill.R

#### Aims
# 1) Calculates skill metrics from simulation outputs

#### Prerequisites
# 1) Analyse simulated data using selected algorithms


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
dv::src()

#### Load data
# sims <- readRDS(here_input("sims.rds"))
sims <- readRDS(here_input("sims-performance.rds"))


#########################
#########################
#### Calculate skill

#### Collate skill scores
# * Loop over sims
# * Read path UD & algorithm UD (if available)
# * Calculate skill for relevant algorithms

# ~ 3.7 mins, 10 cl, 600 performance simulations

tic()
skills <- cl_lapply(.x = split(sims, seq_len(nrow(sims))),
                    .cl = 10L,
                    .chunk = TRUE,
                    .fun = function(sim) {

  #### Define outfile
  out_file <- here_output("skill", paste0(sim$id, ".rds"))
  if (file.exists(out_file)) {
    return(readRDS(out_file))
  }

  #### Read UDs
  # sim <- sims[1, ]
  print(sim$row)
  path   <- terra::rast(here_alg(sim, "path", "ud.tif"))
  null   <- terra::rast(here_input("blank.tif"))
  if (sim$performance) {
    coa_1  <- read_rast(here_alg(sim, "coa", "30 mins", "ud.tif"))
    coa_2  <- read_rast(here_alg(sim, "coa", "120 mins", "ud.tif"))
    # TO DO
    # * Add RSP
  }
  acpfk   <- read_rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-k.tif"))
  acdcpfk <- read_rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-k.tif"))
  # TO DO
  # * Add backward sampler results

  #### Define blank skill table
  skill <- data.table(
    id = sim$id,
    performance = sim$performance,
    alg = c("Null", "COA(30)", "COA(120)", "ACPF(K)", "ACDCPF(K)"),
    mb = NA_real_,
    me = NA_real_,
    rmse = NA_real_,
    R = NA_real_,
    d = NA_real_
  )

  #### Calculate skill metrics
  if (sim$performance) {
    # We will calculate model skill scores for each algorithm
    uds <- list(null, coa_1, coa_2, acpfk, acdcpfk)
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
  saveRDS(skill, out_file)
  skill

}) |> rbindlist()
toc()

saveRDS(skills, here_data("sims", "synthesis", "skill-raw.rds"))


#########################
#########################
#### Process skill scores

#### Drop simulations which failed for all algorithms
pos_success <- rowSums(!is.na(skills[, c("mb", "me", "rmse", "R", "d")])) > 0
skills[!pos_success, ]
skills <- skills[pos_success, ]

#### Add required information
skills <- merge(skills, sims, by = "id")

#### Save processed data
saveRDS(skills, here_data("sims", "synthesis", "skill.rds"))


#### End of code.
#########################
#########################
