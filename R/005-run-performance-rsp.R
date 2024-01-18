#########################
#########################
#### run-performance-rsp.R

#### Aims
# 1) Estimate UDs using COAs()

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
spatw        <- readRDS(here_input("spatw.rds"))
spat_ll_dbbw <- terra::wrap(terra::rast(here_input("spat_ll_dbb.tif")))
tm           <- qs::qread(here_input("actel", "tm.qs"))
sims_for_performance    <- readRDS(here_input("sims-performance.rds"))
sims_for_performance_ls <- split(sims_for_performance, sims_for_performance$id)


#########################
#########################
#### Estimate UDs

gc()
tic()
success <-
  cl_lapply(sims_for_performance_ls,
          .fun = function(sim, .chunkargs) {
            # sim <- sims_for_performance_ls[[1]]
            print(sim$row)
            tic()
            workflow_rsp(sim = sim,
                         spat = .chunkargs$spat,
                         spat_ll_dbb = .chunkargs$spat_ll_dbb,
                         tm = tm)
            toc()
          },
          .chunk = TRUE,
          .chunk_fun = function(sim) {
            .chunkargs <- list(spat = terra::unwrap(spatw),
                 spat_ll_dbb = terra::unwrap(spat_ll_dbbw))
            .chunkargs
          },
          .cl = NULL)
toc()

# Check success
rbindlist(success)
