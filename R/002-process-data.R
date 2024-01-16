#########################
#########################
#### run-algorithms-prep.R

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
detections <- readRDS(here_input("detections.rds"))
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
#### Prep datasets

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

#### Acoustics (~14 s)
# * /combination {system type, path type}/array_type/array_realisation/path_realisation
sims_for_realisations <-
  sims |>
  group_by(combination,
           array_type, array_realisation,
           path_realisation) |>
  slice(1L) |>
  as.data.table()
sims_for_realisations_ls <- split(sims_for_realisations, seq_len(nrow(sims_for_realisations)))
pbapply::pblapply(sims_for_realisations_ls, function(sim) {
  folder <- here_input("acoustics",
                       sim$combination,
                       sim$array_type, sim$array_realisation,
                       sim$path_realisation)
  dir.create(folder, recursive = TRUE)
  out <- get_acoustics(sim, detections)
  qs::qsave(out, file.path(folder, "acoustics.qs"))
  TRUE
}) |> invisible()

#### Paths (~11 s)
# * /combination {system type, path type}/array_type/array_realisation/path_realisation
pbapply::pblapply(sims_for_realisations_ls, function(sim) {
  folder <- function(top) {
    here_input(top,
               sim$combination,
               sim$array_type, sim$array_realisation,
               sim$path_realisation)
  }
  dir.create(folder("paths"), recursive = TRUE)
  acoustics <- qs::qread(file.path(folder("acoustics"), "acoustics.qs"))
  out <- get_path(sim, paths, acoustics)
  qs::qsave(out, file.path(folder("paths"), "path.qs"))
  TRUE
}) |> invisible()

#### Data lists (~50 s)
# * /combination {system type, path type}/array_type/array_realisation/path_realisation/gamma
pbapply::pblapply(sims_for_realisations_ls, function(sim) {
  # sim <- sims_for_realisations_ls[[1]]
  folder <- here_input("dlist",
                       sim$combination,
                       sim$array_type, sim$array_realisation,
                       sim$path_realisation,
                       sim$gamma)
  dir.create(folder, recursive = TRUE)
  acc      <- get_acoustics(sim, detections)
  arc      <- get_path(sim, paths, acc)
  moorings <- get_array(sim, arrays)
  moorings$receiver_range <- sim$gamma
  out <- pat_setup_data(.acoustics = acc,
                        .archival = arc,
                        .moorings = moorings,
                        .bathy = spat,
                        .lonlat = FALSE)
  out$spatial$bathy <- NULL
  qs::qsave(out, file.path(folder, "dlist.qs"))
  TRUE
}) |> invisible()


#########################
#########################
#### Prepare AC* inputs

#### Define directories (~7 s)
# Define folders
# * ac/{array_type}/{array_realisation}/{gamma}/overlaps.qs
# * ac/{array_type}/{array_realisation}/{gamma}/{alpha}/{beta}/kernels.qs
pbapply::pbsapply(split(det_pars_all, seq_len(nrow(det_pars_all))), function(d) {
  dir.create(here_input("ac", d$array_type, d$array_realisation,
                        d$gamma, d$alpha, d$beta),
             recursive = TRUE)
}) |> invisible()

#### Define detection overlaps (6 s, 10 cl; 28 s, 0 cl)
gc()
tic()
pbapply::pblapply(split(gammas, seq_len(nrow(gammas))), cl = 10L, function(d) {
  # Isolate array data
  # d = split(gammas, seq_len(nrow(gammas)))[[1]]
  out_file <- here_input("ac", d$array_type, d$array_realisation, d$gamma, "overlaps.qs")
  if (!file.exists(out_file)) {
    # Define array
    array <-
      arrays[[d$array_type]] |>
      filter(array_id == d$array_realisation) |>
      mutate(receiver_range = d$gamma) |>
      as.data.table()
    # Define data list
    dlist <- pat_setup_data(.moorings = array,
                            .bathy = spat,
                            .lonlat = FALSE)
    # Define detection containers
    overlaps   <- acs_setup_detection_overlaps(dlist)
    qs::qsave(overlaps, out_file)
    TRUE
  }
}) |> invisible()
toc()

#### Define detection kernels (~5 mins with 10 cl)
gc()
tic()
pbapply::pblapply(split(det_pars_all, seq_len(nrow(det_pars_all))), cl = 10L, function(d) {
  # Define array data
  print(paste(d$index, "/", max(det_pars_all$index)))
  out_file <- here_input("ac", d$array_type, d$array_realisation,
                         d$gamma, d$alpha, d$beta, "kernels.qs")
  if (!file.exists(out_file)) {
    # Define array
    array <-
      arrays[[d$array_type]] |>
      filter(array_id == d$array_realisation) |>
      mutate(receiver_range = d$gamma) |>
      as.data.table()
    # Define data list
    dlist <- pat_setup_data(.moorings = array,
                            .bathy = spat,
                            .lonlat = FALSE)
    # Define detection kernels
    kernels <- acs_setup_detection_kernels(dlist,
                                           .alpha = d$alpha, .beta = d$beta,
                                           .verbose = FALSE)
    kernels <- wrapr(kernels)
    qs::qsave(kernels, out_file)
  }
  TRUE
}) |> invisible()
toc()


#########################
#########################
#### Prepare spatspat inputs

im  <- as.im.SpatRaster(spat)
win <- as.owin.SpatRaster(spat, .im = im)

plot(im)
plot(win)

qs::qsave(im, here_input("im.qs"))
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
  dir.create(file.path(top, "rsp"), recursive = TRUE)
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
