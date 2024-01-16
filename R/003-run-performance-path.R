#########################
#########################
#### run-performance-path.R

#### Aims
# 1) Generate UDs for simulated paths

#### Prerequisites
# 1) Simulate data & prepare for algorithm implementation


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
gridw      <- readRDS(here_input("gridw.rds"))
im         <- qs::qread(here_input("im.qs"))
win        <- qs::qread(here_input("win.qs"))
sims_for_performance    <- readRDS(here_input("sims-performance.rds"))
sims_for_performance_ls <- split(sims_for_performance, sims_for_performance$id)


#########################
#########################
#### Estimate UDs

# * ~5 mins, cl = 1L
# * ~40 s, cl = 10L forks

gc()
tic()
cl_lapply(sims_for_performance_ls,
          .fun = function(sim, .chunkargs) {
            # sim <- sims_for_performance_ls[[568]]
            print(sim$row)
            workflow_path(sim,
                          spat = .chunkargs$spat,
                          im = im, win = win)
          },
          .chunk = TRUE,
          .chunk_fun = function(sim) {
            list(spat = terra::unwrap(gridw))
          },
          .cl = 10L)
toc()


#### End of code.
#########################
#########################
