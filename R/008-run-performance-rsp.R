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

# Validate gdistance for RSP
library(gdistance)
Sys.sleep(60)
library(gdistance)

#### Load data
spatw        <- readRDS(here_input("spatw.rds"))
spat_ll_dbbw <- terra::wrap(terra::rast(here_input("spat_ll_dbb.tif")))
tm           <- qs::qread(here_input("actel", "tm.qs"))
sims_for_performance    <- readRDS(here_input("sims-performance.rds"))
sims_for_performance_ls <- split(sims_for_performance, sims_for_performance$id)


#########################
#########################
#### Estimate UDs

#### (optional) Test
test <- FALSE
if (test) {
  cl <- 1L
  sims_for_performance_ls <- sims_for_performance_ls[1:2]
} else {
  cl <- 50L
}

#### Run workflow (~13 h on one cl)
gc()
tic()
success <-
  cl_lapply(sims_for_performance_ls,
          .fun = function(sim, .chunkargs) {
            # sim <- sims_for_performance_ls[[1]]
            print(sim$row)
            tic()
            success <- workflow_rsp(sim = sim,
                                    spat = .chunkargs$spat,
                                    spat_ll_dbb = .chunkargs$spat_ll_dbb,
                                    tm = tm)
            toc()
            success
          },
          .chunk = TRUE,
          .chunk_fun = function(sim) {
            .chunkargs <- list(spat = terra::unwrap(spatw),
                 spat_ll_dbb = terra::unwrap(spat_ll_dbbw))
            .chunkargs
          },
          .cl = cl)
toc()
# beepr::beep(10L)

#### Record success
sdt <- rbindlist(success)
saveRDS(sdt, here_data("sims", "output", "success", "rsp.rds"))

#### Quick checks
if (interactive()) {
  s <- sims_for_performance_ls[[1]]
  pp <- par(mfrow = c(2, 3))
  here_alg(s, "path", "ud.tif")  |> terra_qplot()
  here_alg(s, "coa", "30 mins", "ud.tif") |> terra_qplot()
  here_alg(s, "coa", "120 mins", "ud.tif") |> terra_qplot()
  here_alg(s, "rsp", "default", "ud.tif") |> terra_qplot()
  here_alg(s, "rsp", "custom", "ud.tif") |> terra_qplot()
  m <- read_array(s)
  points(m$receiver_easting, m$receiver_northing)
  par(pp)
}



#### End of code.
#########################
#########################
