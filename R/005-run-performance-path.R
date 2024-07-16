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

# spat <- terra::unwrap(spatw)

#########################
#########################
#### Estimate UDs

#### Run workflow (~62 s, 10 cl)
gc()
tic()
cl_lapply(sims_for_performance_ls,
          .fun = function(sim, .chunkargs) {
            # sim <- sims_for_performance_ls[[568]]
            print(sim$row)
            workflow_path(sim,
                          spat = .chunkargs$spat,
                          win = win)
          },
          .chunk = TRUE,
          .chunk_fun = function(sim) {
            list(spat = terra::unwrap(spatw))
          },
          .cl = 10L)
toc()

#### Quick checks
if (interactive()) {
  pp <- par(mfrow = c(3, 3))
  nms <- sample(names(sims_for_performance_ls), 9)
  for (nm in nms) {
    s <- sims_for_performance_ls[[nm]]
    here_alg(s, "path", "ud.tif")  |> terra_qplot()
    read_path(s) |> select(x, y) |> as.matrix() |> points(pch = ".")
  }
  par(pp)
}


#### End of code.
#########################
#########################
