#########################
#########################
#### run-sims.R

#### Aims
# 1) Simulate data

#### Prerequisites
# 1) NA

#### TO DO
# * Run time trials for 005-workflow.R
# * Discuss cluster options with SD


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
library(prettyGraphics)
library(tictoc)

#### Load data
source(here_r("001-define-global-param.R"))
source(here_r("002-define-helpers.R"))
grid       <- terra::rast(here_input("grid.tif"))
arrays     <- readRDS(here_input("arrays.rds"))
paths      <- readRDS(here_input("paths.rds"))
detections <- readRDS(here_input("detections.rds"))
sims       <- readRDS(here_input("sims.rds"))


#########################
#########################
#### Implement simulations

#### (optional) Use a subset of data for testing
sims <- sims[which(sims$performance)[1:3], ]

#### Set up cluster
cl <- NULL
workers <- 2L
cl <- workers
cl <- parallel::makeCluster(workers)
parallel::clusterEvalQ(cl, {
  library(patter)
  library(data.table)
  library(dtplyr)
  library(dplyr, warn.conflicts = TRUE)
})
parallel::clusterExport(cl,
                        varlist = c("paths", "detections", "seed",
                                    "unwrapr", "mins",
                                    "here_input", "here_output",, "here_alg",
                                    "pf_coords", "skill_by_alg"))

#### Run simulations
# (optional) TO DO
# * Define chunks to loop over in parallel & create a log file for each chunk
pbapply::pblapply(split(sims, seq_len(nrow(sims))), cl = cl, function(sim) {

  # Run sim function
  workflow(
    sim = sim,
    paths = paths,
    detections = detections,
    seed = seed,
    # algorithm parameters
    sigma = spatstat.explore::bw.diggle,
    # utility functions, passed as arguments
    unwrapr = unwrapr,
    mins = mins,
    here_input = here_input,
    here_output = here_output,
    here_alg = here_alg,
    pf_coords = pf_coords,
    skill_by_alg = skill_by_alg
  )

}) |> invisible()


#### End of code.
#########################
#########################
