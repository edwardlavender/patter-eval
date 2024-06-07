#########################
#########################
#### process-data.R

#### Aims:
# (1) Prepare for algorithm implementations
# * This code is implemented prior to simulations for improved speed

#### Prerequisites
# 1) Simulate data (sim-data.R)


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
spat       <- terra::rast(here_input("spat.tif"))
arrays     <- readRDS(here_input("arrays.rds"))
paths      <- readRDS(here_input("paths.rds"))
acoustics  <- qs::qread(here_input("acoustics.qs"))
detections <- qs::qread(here_input("detections.qs"))
sims       <- readRDS(here_input("sims.rds"))


#########################
#########################
#### Simulations

#### Identify performance simulations
# Isolate relevant simulations
sims_for_performance <-
  sims |>
  filter(performance) |>
  mutate(row = row_number()) |>
  as.data.table()
# Save
nrow(sims_for_performance) # 599
saveRDS(sims_for_performance, here_input("sims-performance.rds"))


#########################
#########################
#### Parameters

#### Isolate detection parameters
det_pars_all <-
  sims |>
  select(array_type, array_realisation, alpha, beta, gamma) |>
  distinct() |>
  mutate(index = row_number()) |>
  as.data.table()

#### Define gamma parameters
gammas <-
  det_pars_all |>
  select(array_type, array_realisation, gamma) |>
  distinct() |>
  as.data.table()


#########################
#########################
#### Prepare datasets

#
#
# TO DO
# REVIEW THIS CODE WITH RESPECT TO WHAT WE NEED
#
#

#### Arrays (~5 s)
# Define relevant simulation information
# * /array_type/array_realisation/
sims_for_array <-
  sims |>
  group_by(array_type, array_realisation) |>
  slice(1L) |>
  as.data.table()
# Define data
sims_for_arrays_ls <- split(sims_for_array, seq_len(nrow(sims_for_array)))
pbapply::pblapply(sims_for_arrays_ls, function(sim) {
  folder <- here_input("arrays", sim$array_type, sim$array_realisation)
  dir.create(folder, recursive = TRUE)
  out <- get_array(sim, arrays)
  qs::qsave(out, file.path(folder, "array.qs"))
  TRUE
}) |> invisible()

#### Acoustics (0, 1) & detections (1) (~1.5 mins, 10 cl)
# * /combination {system type, path type}/array_type/array_realisation/path_realisation
# * acoustics (0, 1) is required for particle algorithms
# * detections (1) is required for coa()/rsp()
sims_for_realisations <-
  sims |>
  group_by(combination,
           array_type, array_realisation,
           path_realisation) |>
  slice(1L) |>
  as.data.table()
sims_for_realisations_ls <-
  split(sims_for_realisations, seq_len(nrow(sims_for_realisations)))
pbapply::pblapply(sims_for_realisations_ls, cl = 10L, function(sim) {
  folder <- here_input("acoustics",
                       sim$combination,
                       sim$array_type, sim$array_realisation,
                       sim$path_realisation)
  dir.create(folder, recursive = TRUE)
  acc <- get_acoustics(sim, acoustics)
  det <- get_detections(sim, detections)
  qs::qsave(acc, file.path(folder, "acoustics.qs"))
  qs::qsave(det, file.path(folder, "detections.qs"))
  TRUE
}) |> invisible()

#### Paths (~11 * 2 s)
# * /combination {system type, path type}/array_type/array_realisation/path_realisation
pbapply::pblapply(sims_for_realisations_ls, cl = 10L, function(sim) {
  folder <- function(top) {
    here_input(top,
               sim$combination,
               sim$array_type, sim$array_realisation,
               sim$path_realisation)
  }
  dir.create(folder("paths"), recursive = TRUE)
  det <- get_detections(sim, detections)
  out <- get_path(sim, paths, det)
  qs::qsave(out, file.path(folder("paths"), "path.qs"))
  TRUE
}) |> invisible()


#########################
#########################
#### Package objects

#### Actel objects (~41 s)
# This code has to be run non interactively to avoid prompts.
tic()
dir.create(here_data("sims", "input", "actel"))
system("R CMD BATCH --no-save --no-restore ./R/004-process-data-rsp.R ./data/sims/input/actel/log.txt")
toc()

#### Patter objects
# Yobs
# TO DO


#########################
#########################
#### Prepare spatspat inputs

im  <- as.im.SpatRaster(spat)
win <- as.owin.SpatRaster(spat, .im = im)

plot(im)
plot(win)

qs::qsave(win, here_input("win.qs"))


#########################
#########################
#### Prepare output directories

#### Build folders for algorithm outputs (~25 s)
# Define data
sims_by_realisation <-
  sims |>
  group_by(combination, array_type, array_realisation, path_realisation) |>
  slice(1L) |>
  as.data.table()
# Build folders
# * output/{combination}/{array_type}/{array_realisation}/{path_realisation}/{algorithm}
pbapply::pblapply(split(sims_by_realisation, seq_len(nrow(sims_by_realisation))), function(sim) {
  top <- here_alg(sim)
  dir.create(top, recursive = TRUE)
  dir.create(file.path(top, "path"), recursive = TRUE)
  dir.create(file.path(top, "coa", "30 mins"), recursive = TRUE)
  dir.create(file.path(top, "coa", "120 mins"), recursive = TRUE)
  dir.create(file.path(top, "rsp", "default"), recursive = TRUE)
  dir.create(file.path(top, "rsp", "custom"), recursive = TRUE)
  dir.create(file.path(top, "patter", "acpf"), recursive = TRUE)
  dir.create(file.path(top, "patter", "acdcpf"), recursive = TRUE)
}) |> invisible()

#### Build folders for patter outputs (~25 s)
# Validate that each simulation is uniquely defined by the realisations & the algorithm parameters
sims |>
  group_by(combination, array_type, array_realisation, path_realisation, alg_par) |>
  summarise(n = n()) |>
  pull(n) |>
  table()
pbapply::pblapply(split(sims, seq_len(nrow(sims))), function(sim) {
  top <- file.path(here_alg(sim), "patter")
  dir.create(file.path(top, "acpf", sim$alg_par), recursive = TRUE)
  dir.create(file.path(top, "acdcpf", sim$alg_par), recursive = TRUE)
}) |> invisible()


#### End of code.
#########################
#########################
