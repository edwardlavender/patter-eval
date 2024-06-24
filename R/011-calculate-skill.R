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

#### (optional) Test
test <- FALSE
if (test) {
  sims  <- sims[1:2, ]
  cl    <- 1L
  chunk <- FALSE
} else {
  cl    <- 10L
  chunk <- TRUE
}
# (optional) Debug with parallel = FALSE
parallel <- TRUE
if (!parallel) {
  cl    <- NULL
  chunk <- FALSE
}

#### Collate skill scores
# * Loop over sims
# * Read path UD & algorithm UD (if available)
# * Calculate skill for relevant algorithms

# ~10.8 mins, 1181 performance simulations

tic()
gc()
skills <- cl_lapply(.x = split(sims, seq_len(nrow(sims))),
                    .cl = cl,
                    .chunk = chunk,
                    .fun = function(sim) {

  #### Define outfile
  # sim <- sims[1, ]
  out_file  <- here_output("skill", paste0(sim$id, ".rds"))
  overwrite <- FALSE
  if (!overwrite && file.exists(out_file)) {
    return(readRDS(out_file))
  }

  #### Read UDs
  print(sim$row)
  path   <- terra::rast(here_alg(sim, "path", "ud.tif"))
  null   <- terra::rast(here_input("blank.tif"))
  if (sim$performance) {
    coa_1  <- read_rast(here_alg(sim, "coa", "30 mins", "ud.tif"))
    coa_2  <- read_rast(here_alg(sim, "coa", "120 mins", "ud.tif"))
    rsp_1  <- read_rast(here_alg(sim, "rsp", "default", "ud.tif"))
    rsp_2  <- read_rast(here_alg(sim, "rsp", "custom", "ud.tif"))
  }
  acpff   <- read_rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-f.tif"))
  acdcpff <- read_rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-f.tif"))
  acpfs   <- read_rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-s.tif"))
  acdcpfs <- read_rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-s.tif"))

  #### Define blank skill table
  skill <- data.table(
    id = sim$id,
    performance = sim$performance,
    alg = c("Null", "COA(30)", "COA(120)", "RSP(1)", "RSP(2)", "ACPF(F)", "ACDCPF(F)", "ACPF(S)", "ACDCPF(S)"),
    mb = NA_real_,
    me = NA_real_,
    rmse = NA_real_,
    R = NA_real_,
    d = NA_real_
  )

  #### Calculate skill metrics
  if (sim$performance) {
    # We will calculate model skill scores for each algorithm
    uds <- list(null, coa_1, coa_2, rsp_1, rsp_2, acpff, acdcpff, acpfs, acdcpfs)
  } else {
    # We will calculate model skill scores for patter algorithms only
    uds <- list(null, NULL, NULL, NULL, NULL, acpff, acdcpff, acpfs, acdcpfs)
  }
  skill[, mb := skill_by_alg(uds, path, skill_mb)]
  skill[, me   := skill_by_alg(uds, path, skill_me)]
  skill[, rmse := skill_by_alg(uds, path, skill_rmse)]
  skill[, R := skill_by_alg(uds, path, skill_R)]
  skill[, d := skill_by_alg(uds, path, skill_d)]

  #### Return outputs
  saveRDS(skill, out_file)
  skill

}) |> rbindlist()
toc()

saveRDS(skills, here_output("synthesis", "skill-raw.rds"))


#########################
#########################
#### Process skill scores

#### Drop simulations which failed for all algorithms
nrow(skills)
pos_success <- rowSums(!is.na(skills[, c("mb", "me", "rmse", "R", "d")])) > 0
skills[!pos_success, ]
skills <- skills[pos_success, ]
nrow(skills)

#### Add required information
skills <- merge(skills, sims, by = "id")

#### Save processed data
saveRDS(skills, here_output("synthesis", "skill.rds"))


#### End of code.
#########################
#########################
