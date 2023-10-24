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
source(here_r("005-workflow.R"))
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
grid <- terra::wrap(grid)
cl <- NULL
workers <- 2L
if (.Platform$OS.type == "unix") {
  cl <- workers
} else if (.Platform$OS.type == "windows") {
  cl <- parallel::makeCluster(workers)
  parallel::clusterEvalQ(cl, {
    library(patter)
    library(data.table)
    library(dtplyr)
    library(dplyr, warn.conflicts = TRUE)
  })
  parallel::clusterExport(cl,
                          varlist = c("workflow",
                                      "grid", "paths", "detections", "seed",
                                      "unwrapr", "mins",
                                      "here_input", "here_output", "here_alg",
                                      "pf_coords", "skill_by_alg"))
}

#### Run simulations
# (optional) TO DO
# * Define chunks to loop over in parallel & create a log file for each chunk
pbapply::pblapply(split(sims, seq_len(nrow(sims))), cl = cl, function(sim) {

  # Run sim function
  # sim <- sims[1, ]
  success <-
    workflow(
    sim = sim,
    grid = grid,
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
    skill_by_alg = skill_by_alg,
    test = TRUE
  )

  '
  # Check outputs
  rstudioapi::navigateToFile(here_output("algorithm", sim$id, "log.txt"))
  readRDS(here_output("skill", paste0(sim$id, ".rds")))
  readRDS(here_output("time", paste0(sim$id, ".rds")))
  '

  success

}) |> invisible()


#### End of code.
#########################
#########################
