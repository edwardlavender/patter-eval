#########################
#########################
#### process-data-density.R

#### Aims
# 1) Improve speed of the backward sampler via precalculation of movement densities:
# * The backward sampler is expensive.
# * If we calculate densities on the fly, each operation takes two minutes.
# * For ~40000 simulations, ETA is: 40000 * 2/60/24 = 55 days on one core.
# * An alternative is to pre-calculate movement densities.
# * This is an expensive but one-off operation:
# * ETA: 8 hours on one core (feasible);
# * Output: one million density files

#### Prerequisites
# 1) NA


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
library(prettyGraphics)
library(tictoc)
dv::src()

#### Load data
spat      <- terra::rast(here_input("spat.tif"))
path_pars <- readRDS(here_input("path_pars.rds"))


#########################
#########################
#### File size estimation

if (interactive()) {

  # Here we check whether pre-calculation is feasible, given the size of the grid.
  # (For very large grids this is not possible b/c the disk space required
  # to store densities becomes prohibitive.)

  #### Define an example location
  # terra::plot(spat)
  xy <- cbind(5473.568, 5436.858)
  pt_cell <- terra::cellFromXY(spat, xy)
  pt_xy   <- terra::xyFromCell(spat, pt_cell)

  #### Define reachable zone around point
  tic()
  i <- 1
  zone <-
    xy |>
    terra::vect() |>
    terra::buffer(width = path_pars$mobility[i]) |>
    sf::st_as_sf()

  #### Define the probability density of movements into all reachable cells
  reachable <-
    # Identify reachable cells
    exactextractr::exact_extract(spat, zone,
                                 include_cell = TRUE,
                                 include_cols = NULL,
                                 include_xy = TRUE)[[1]] |>
    as.data.table() |>
    # Compute densities between cells
    mutate(cell = as.integer(cell),
           pt_cell = pt_cell,
           pt_x = pt_xy[, 1],
           pt_y = pt_xy[, 2],
           len = patter:::dist_2d(pt_x, pt_y, x, y),
           dens = dtruncgamma(len,
                              .shape = path_pars$shape[i],
                              .scale = path_pars$scale[i],
                              .mobility = path_pars$mobility[i])) |>
    # Retain only essential columns
    select(cell, dens) |>
    as.data.table()
  toc()

  #### Write to file and check size
  # Compare multiple formats
  mapply(FUN = function(extension, write) {
    ftmp <- tempfile(fileext = extension)
    write(reachable, ftmp)
    # Get file size (MB)
    file.size(ftmp) / 1e6L
  },
  c(".rds", ".parquet", ".feather", ".fst", ".qs"),
  list(saveRDS, arrow::write_parquet, feather::write_feather, fst::write_fst, qs::qsave)) |>
    sort()
  # qs::qsave() produces the smallest files (followed by saveRDS)
  fmb <- 0.028910 # MB

  #### Estimate approx total size of all files required to store densities (GB):
  # We will have one file for every grid cell and every set of movement parameters
  (fmb * terra::ncell(spat) * nrow(path_pars)) / 1e3L

}


#########################
#########################
#### Implement precalculation

#### Create directories
dir.create(here_input("density", "1"), recursive = TRUE)
dir.create(here_input("density", "2"), recursive = TRUE)

#### Define spat coordinates (~6 s)
tic()
spat_dt <-
  terra::as.data.frame(spat, xy = TRUE, cell = TRUE) |>
  select(cell, x, y) |>
  mutate(file_1 = here_input("density", "1", paste0(cell, ".qs")),
         file_2 = here_input("density", "2", paste(cell, ".qs"))) |>
  as.data.table()
# Define coordinate matrix
spat_xy <- cbind(spat_dt$x, spat_dt$y)
toc()

#### Define test
test <- FALSE
if (test) {
  spat_xy <- spat_xy[1:1000, ]
  cl <- 10L
} else {
  cl <- 50L
}
nrow(spat_xy)

#### Precalculate densities (1): ~17 mins (50 cl)
# Note that mobility is increased by one grid cell
tic()
if (FALSE) {
  gc()
  log.txt <- here_data("sims", "output", "log", "server",
                       "004-process-data-density-progress-1.txt")
  log.con <- log_open(log.txt)
  cl_lapply(seq_len(nrow(spat_xy)), .fun = function(i) {
    spatDens(.spat = spat,
             .xy = spat_xy[i, , drop = FALSE],
             .shape = path_pars$shape[1],
             .scale = path_pars$scale[1],
             .mobility = path_pars$mobility[1] + sr,
             .file = spat_dt$file_1[i])
  },
  .chunk = TRUE,
  .cl = cl)
  log_close(log.con)
  # rstudioapi::navigateToFile(log.txt)
}
toc()

#### Precalculate densities (2): ~19 mins (50 cl)
gc()
tic()
if (TRUE) {
  log.txt <- here_data("sims", "output", "log", "server",
                       "004-process-data-density-progress-2.txt")
  log.con <- log_open(log.txt)
  cl_lapply(seq_len(nrow(spat_xy)), .fun = function(i) {
    spatDens(.spat = spat,
             .xy = spat_xy[i, , drop = FALSE],
             .shape = path_pars$shape[2],
             .scale = path_pars$scale[2],
             .mobility = path_pars$mobility[2] + sr,
             .file = spat_dt$file_2[i])
  },
  .chunk = TRUE,
  .cl = cl)
  log_close(log.con)
  # rstudioapi::navigateToFile(log.txt)
}
toc()


#### End of code.
#########################
#########################
