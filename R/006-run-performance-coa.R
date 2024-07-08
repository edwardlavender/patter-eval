#########################
#########################
#### run-performance-coa.R

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

#### Essential packages
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tictoc)
dv::src()

#### Load data
spatw      <- readRDS(here_input("spatw.rds"))
win        <- qs::qread(here_input("win.qs"))
sims_for_performance    <- readRDS(here_input("sims-performance.rds"))
sims_for_performance_ls <- split(sims_for_performance, sims_for_performance$id)


#########################
#########################
#### Estimate UDs

#### (optional) Test
test <- FALSE
if (test) {
  cl <- NULL
  sims_for_performance_ls <- sims_for_performance_ls[1:2]
} else {
  cl <- 10L
}

#### Run workflow (~1.5 mins)
gc()
tic()
success <- cl_lapply(sims_for_performance_ls,
                     .fun = function(sim, .chunkargs) {
                       # sim <- sims_for_performance_ls[[1]]
                       print(sim$row)
                       workflow_coa(sim = sim,
                                    spat = .chunkargs$spat,
                                    win = win)
                     },
                     .chunk = TRUE,
                     .chunk_fun = function(sim) {
                       list(spat = terra::unwrap(spatw))
                     },
                     .cl = cl)
toc()

#### Record success
sdt <- rbindlist(success)
sdt[coa_1 == FALSE, ]
sdt[coa_2 == FALSE, ]
saveRDS(sdt, here_data("sims", "output", "success", "coa.rds"))

#### Quick checks
if (interactive()) {
  s <- sims_for_performance_ls[[10]]
  pp <- par(mfrow = c(1, 3))
  here_alg(s, "path", "ud.tif")  |> terra_qplot()
  here_alg(s, "coa", "30 mins", "ud.tif") |> terra_qplot()
  here_alg(s, "coa", "120 mins", "ud.tif") |> terra_qplot()
  m <- read_array(s)
  points(m$receiver_x, m$receiver_y)
  par(pp)
}


#### End of code.
#########################
#########################
