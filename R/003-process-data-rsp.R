#########################
#########################
#### process-data-rsp.R

#### Aims:
# (1) Prepare for algorithm implementations (RSP)

#### Prerequisites
# 1) This code is implemented non interactively via process-data.R


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
detections <- readRDS(here_input("detections.rds"))
sims       <- readRDS(here_input("sims.rds"))


#########################
#########################
#### Actel objects

#### Define simulations (as in process-data.R)
sims_for_realisations <-
  sims |>
  group_by(combination,
           array_type, array_realisation,
           path_realisation) |>
  slice(1L) |>
  as.data.table()
sims_for_realisations_ls <-
  split(sims_for_realisations, seq_len(nrow(sims_for_realisations)))

#### Build explore objects (~41 s)
# * /combination {system type, path type}/array_type/array_realisation/path_realisation

actel <-
  pbapply::pblapply(sims_for_realisations_ls, function(sim) {

    #### Standard set up
    # sim = sims_for_realisations_ls[[1]]
    folder <- here_input("actel",
                         sim$combination,
                         sim$array_type, sim$array_realisation,
                         sim$path_realisation)
    dir.create(folder, recursive = TRUE)

    #### Get data
    acoustics <- get_acoustics(sim, detections) |> as.data.frame()
    moorings  <- get_array(sim, arrays) |> as.data.frame()

    #### Prepare actel dataserts
    # Biometrics
    act_bio <-
      data.frame(Release.date = acoustics$timestamp[1] - 1,
                 Release.site = "unspecified",
                 Group = 1,
                 Signal = 1L)
    # Deployments
    act_deployments <-
      moorings |>
      mutate(Start = paste(receiver_start, "00:00"),
             Stop = paste(receiver_end + 1, "00:00")) |>
      select(Receiver = receiver_id,
             Station.name = receiver_id,
             Start, Stop)
    # Spatial datasets
    moorings_ll <-
      moorings |>
      select(receiver_easting, receiver_northing) |>
      as.matrix() |>
      terra::vect(crs = terra::crs(spat)) |>
      terra::project("EPSG: 4326") |>
      terra::crds()
    moorings$lon <- moorings_ll[, 1]
    moorings$lat <- moorings_ll[, 2]
    act_spatial <-
      moorings |>
      mutate(
        Longitude = lon, Latitude = lat,
        x = receiver_easting, y = receiver_northing,
        Array = "A0", Section = "unspecified",
        Type = "Hydrophone", Range = sim$gamma) |>
      select(Station.name = receiver_id,
             Longitude, Latitude, x, y,
             Array, Section, Type, Range)
    # Detections
    act_detections <-
      acoustics |>
      mutate(Receiver = receiver_id,
             Timestamp = timestamp,
             CodeSpace = "unspecified",
             Signal = 1L) |>
      select(Receiver, Timestamp, CodeSpace, Signal)

    #### Build explore object
    act <- actel::preload(act_bio, act_spatial,
                          act_deployments, act_detections,
                          tz = "UTC")
    act <- actel::explore(act, tz = "UTC", GUI = "never")
    qs::qsave(act, file.path(folder, "actel.qs"))

    TRUE
  }) |> invisible()

#### Build transitionMatrix (~1 s)
tic()
spat_ll <- terra::project(spat, "EPSG:4326")
tm <- actel::transitionLayer(raster::raster(spat_ll))
qs::qsave(tm, here_input("actel", "tm.qs"))
toc()


#### End of code.
#########################
#########################
