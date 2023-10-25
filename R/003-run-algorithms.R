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

#### Path UDs
# * ~3 mins, cl = 1L
# * ~40 s, cl = 10L forks
gc()
pbapply::pblapply(sims_for_performance_ls, cl = 10L, function(sim) {
  # sim <- sims_for_performance_ls[[1]]
  grid <- terra::unwrap(grid)
  workflow_path(sim, grid)
}) |> invisible()

#### COA UDs
# * ~3 mins, cl = 1L
# * ~40 s, cl = 10L forks
gc()
pbapply::pblapply(sims_for_performance_ls, cl = 10L, function(sim) {
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
m <- read_array(s)
points(m$receiver_easting, m$receiver_northing)


#########################
#########################
#### Patter simulations

#### Overview
# Here, {patter} workflows are implemented (for both performance & sensitivity simulations)

#### Set up cluster
# Number of forks
cl <- 2L
# Define chunks to iterate over in parallel
chunks <- cl_chunks(cl, nrow(sims))
# Creare a log file for each chunk
lapply(seq_len(cl), function(i) {
  file.create(here_output("logs", paste0("log-", i, ".txt")))
}) |> invisible()

#### Time trials
# Estimated duration (days) on {cl} CPUs, assuming simulations take {guess} seconds
guess <- 30 # 30 s
(nrow(sims) * 30)/60/60/24/cl

#### Implementation

# Performance
# * All performance simulations should converge!
# sims <- sims_for_performance

gc()
success <-
  pbapply::pblapply(1, cl = NULL, function(i) {
    sink(here_output("logs", paste0("log-", i, ".txt")))
    cat(paste("CHUNK" , i, "\n"))
    t1_chunk <- Sys.time()
    cat(paste("Start:", as.character(t1_chunk), "\n"))
    grid <- terra::unwrap(grid)
    sims_for_chunk <- sims[chunks[[i]], ]
    success <-
      lapply(split(sims_for_chunk, seq_len(nrow(sims_for_chunk)))[1:5], function(sim) {
        cat(paste0("\n", sim$id, ":\n"))
        t1 <- Sys.time()
        workflow_patter(sim, grid)
        t2 <- Sys.time()
        cat(as.numeric(difftime(t2, t1, units = "secs")))
      }) |> unlist()
    t2_chunk <- Sys.time()
    cat("\n\n")
    cat(paste("End:", as.character(t2_chunk), "\n"))
    cat(paste("Duration:", round(difftime(t2_chunk, t1_chunk, units = "hours"), 2), "hours."))
    sink()
    success
  }) |> unlist()


#########################
#########################
#### Quick visuals

qplot <- function(file) {
  terra_qplot(file)
  p <- read_path(sim)
  prettyGraphics::add_sp_path(p$x, p$y, length = 0)
  m <- read_array(sim)
  points(m$receiver_easting, m$receiver_northing)
}

sim <- sims[4, ]

pp <- par(mfrow = c(2, 2))
qplot(here_alg(sim, "path", "ud.tif"))
qplot(here_alg(sim, "coa", "120 mins", "ud.tif"))
qplot(here_alg(sim, "patter", "acpf", sim$alg_par, "ud.tif"))
qplot(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud.tif"))
par(pp)


#### End of code.
#########################
#########################
