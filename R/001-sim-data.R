#########################
#########################
#### sim-data.R

#### Aims
# 1) Simulate datasets for analysis

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
library(collapse)
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(glue)
library(JuliaCall)
library(prettyGraphics)
library(tictoc)
dv::src()
tic()

#### patter
julia_connect()
set_seed()


#########################
#########################
#### High-level directory structure

# Data directories
dir.create(here_data())
dir.create(here_data("rsync"))
dir.create(here_data("sims"))
dir.create(here_data("sims", "input"))
dir.create(here_data("sims", "output"))
dir.create(here_data("sims", "output", "log"))
dir.create(here_data("sims", "output", "log", "patter"))
dir.create(here_data("sims", "output", "log", "patter", "performance"))
dir.create(here_data("sims", "output", "log", "patter", "sensitivity"))
dir.create(here_data("sims", "output", "log", "server"))
dir.create(here_data("sims", "output", "run"))
dir.create(here_data("sims", "output", "skill"))
dir.create(here_data("sims", "output", "success"))
dir.create(here_data("sims", "output", "synthesis"))


#########################
#########################
#### Simulate study

#### Define study site
# We define a simple rectangular study site
# * Grid resolution affects the speed for some routines
# * And required disk space
spat <- terra::rast(vals = 1,
                    res = sr,
                    xmin = 0, xmax = 1e4,
                    ymin = 0, ymax = 1e4,
                    crs = crs)
terra::ncell(spat)
# Define uniform ('blank') SpatRaster
# (used as a NULL model)
spat <- spat / terra::global(spat, "sum")[1, 1]
terra::writeRaster(spat, here_input("blank.tif"), overwrite = TRUE)
# Generate hypothetical bathymetry values (deterministically)
g       <- terra::as.data.frame(spat, xy = TRUE)
g$depth <- gen_depth(g)
spat    <- terra::rasterize(as.matrix(g[, c("x", "y")]),
                            spat,
                            values = g$depth)
names(spat) <- "bathy"
terra::global(terra::cellSize(spat, unit = "km"), "sum")
# Export map
set_map(spat)
# Write SpatRaster(s)
# * Include wrapped version for parallelised routines
saveRDS(terra::wrap(spat), here_input("spatw.rds"))
terra::writeRaster(spat, here_input("spat.tif"), overwrite = TRUE)
ext <- terra::ext(spat)
terra::ncell(spat)

#### Define study period (`timeline`)
# Define number of days
n_days  <- 2
# Define number of two minute time steps in n_days
step    <- "2 mins"
n_steps <- 60/2*24 * n_days
period  <- seq.POSIXt(as.POSIXct("2023-01-01", tz = "UTC"),
                      length.out = n_steps,
                      by = step)
range(period)


#########################
#########################
#### Outline simulations

# We will simulate:
# * one 'study'          (one area, one study period)
# * n_systems            (different detection probability parameters; n = 2)
# * n_arrays             (different arrays; n = 20)
# * n_array_realisations (different realisations; n = 1)
# * n_paths              (different movement parameters, one type per system; n = 2)
# * n_path_realisations  (different realisations; n = 30)
# * simulated dataset for each array/path realisation

# We then generate:
# * modelled outcomes for each simulated dataset & parameter values


#########################
#########################
#### Simulate arrays

#### Overview
# We simulate arrays that vary depending on:
# * Detection probability parameters
# * Receiver placement (regular, random)
# * Receiver number (5, 10, 15, ...)

#### Define detection pars
detection_pars <- data.table(alpha = c(4, 5),
                             beta = c(-0.01, -0.025),
                             gamma = c(750, 500))
saveRDS(detection_pars, here_input("detection_pars.rds"))
# Visualise detection probability function
i <- 2
x <- seq(0, 1000, by = 1)
y <- ddetlogistic(x,
                  .alpha = detection_pars$alpha[i],
                  .beta = detection_pars$beta[i],
                  .gamma = detection_pars$gamma[i])
plot(x, y, type = "l")
# Define the number of systems
n_systems <- nrow(detection_pars)

#### Define array parameters (receiver number, arrangement)
# The number of receivers is chosen such that we have low to close to ~100 % coverage
# (at least within receiver detection ranges)
array_pars <-
  CJ(arrangement = c("regular", "random"),
              number = seq(10, 100, by = 10)) |>
  arrange(arrangement, number) |>
  mutate(index = row_number()) |>
  as.data.table()

#### Simulate arrays
# We simulate a list of n_array designs
# ... as defined by the parameters above (one realisation of each)
# We hold array designs constant across different detection probability parameters
n_array              <- nrow(array_pars)
n_array_realisations <- 1L
tic()
arrays <-
  pbapply::pblapply(split(array_pars, seq_len(nrow(array_pars))), function(d) {
    a <- sim_array(.map = spat,
                   .timeline = period,
                   .arrangement = d$arrangement,
                   .n_array = n_array_realisations,
                   .n_receiver = d$number,
                   .receiver_alpha = NULL,
                   .receiver_beta = NULL,
                   .receiver_gamma = NULL,
                   .plot = FALSE)
    a$n_receiver  <- d$number
    a$arrangement <- d$arrangement
    a
})
head(arrays[[1]])
toc()
saveRDS(arrays, here_input("arrays.rds"))

#### Examine range in array coverage (~5 s)
area_total <- terra::expanse(spat, transform = TRUE)$area
coverage <-
  lapply(split(detection_pars, seq_len(nrow(detection_pars))), function(d) {
    pbapply::pblapply(seq_len(length(arrays)), function(i) {
      # Define detection containers
      containers <- terra::setValues(spat, NA)
      array      <- arrays[[i]]
      array[, cell_id := terra::cellFromXY(containers, cbind(receiver_x, receiver_y))]
      containers[array$cell_id] <- 1
      containers <- terra::buffer(containers, d$gamma)
      containers <- terra::classify(containers, cbind(0, NA))
      # Calculate area in containers
      area_in_containers <- terra::cellSize(containers)
      area_in_containers <- terra::mask(area_in_containers, containers)
      area_in_containers <- terra::global(area_in_containers, "sum", na.rm = TRUE)$sum
      # Calculate % area in containers
      data.table(array_type = i,
                 n_receiver = array$n_receiver[1],
                 arrangement = array$arrangement[1],
                 pc = (area_in_containers / area_total) * 100)
    }) |> rbindlist()
  })
coverage[[1]]; coverage[[2]]
range(coverage[[1]]$pc); range(coverage[[2]]$pc)
# Create tidy table
array_table <- coverage[[1]]
array_table$`Coverage (750)` <- tidy_numbers(coverage[[1]][, .(pc)], 2)
array_table$`Coverage (500)` <- tidy_numbers(coverage[[2]][, .(pc)], 2)
array_table |>
  arrange(arrangement) |>
  mutate(arrangement = stringr::str_to_sentence(arrangement)) |>
  select(Array = array_type,
         Arrangement = arrangement,
         `Receiver` = n_receiver,
         `Coverage (750)`,
         `Coverage (500)`) |>
  tidy_write(here_fig("arrays.txt"))


#########################
#########################
#### Simulate movement

#### Define path parameters
# For each system (set of detection_pars), we consider one path type.
# For convenience, we assume:
# * Gamma distribution of step lengths (with varying parameters)
# * Default settings for turning angles
# (We examine algorithm performance for each pair of system/path types
# ... since we only simulate multiple systems/paths to check the validity
# ... of conclusions, we do not consider all combinations of systems/paths,
# ... to reduce the number of simulations required & to improve speed)
path_pars <- data.table(mobility = c(500, 750),
                        shape = c(1, 15),
                        scale = c(250, 15))
saveRDS(path_pars, here_input("path_pars.rds"))
# Examine simulated paths
i <- 1L
hist(rtruncgamma(.n = 1e5,
                 .shape = path_pars$shape[i],
                 .scale = path_pars$scale[i],
                 .mobility = path_pars$mobility[i]),
     xlim = c(0, path_pars$mobility[i]),
     breaks = 100,
     xlab = "Step length",
     main = i)

#### Simulate paths & associated observations (e.g., depths)
# We simulate n_path types
n_path              <- nrow(path_pars)
stopifnot(n_systems == n_path)
n_path_realisations <- 30L
origin <- matrix(mean(ext[1:2], mean(ext[3:4])), ncol = 2)
origin <- data.table(map_value = terra::extract(spat, origin)[1, 1],
                     x = origin[1], y = origin[2])
# Simulate paths (~3 s)
# * One element for each set of parameters
#   * One data.table for all realisations
tic()
paths <- lapply(split(path_pars, seq_len(nrow(path_pars))), function(d) {
  # Define movement model
  model_move <-
    move_xy(dbn_length = glue("truncated(Gamma({d$shape}, {d$scale}),
                                         upper = {d$mobility})"),
            dbn_angle = "Uniform(-pi, pi)")
  # Generate paths
  sim_path_walk(.map = spat,
                .timeline = period,
                .state = "StateXY",
                .xinit = origin,
                .model_move = model_move,
                .n_path = n_path_realisations,
                .plot = FALSE) |>
    # Simulate depths
    mutate(depth = sim_depth(map_value)) |>
    as.data.table()
})
toc()
stopifnot(length(unique(paths[[1]]$path_id)) == n_path_realisations)
saveRDS(paths, here_input("paths.rds"))


#########################
#########################
#### Simulate acoustic observations (0, 1)

#### Overview
# For each array/path realisation, we simulate acoustic observations
# This code returns a list:
# * One element for each set of detection pars & path type
#   * One element for each array design
#       * One data.table with acoustic observations for all array/path realisations

#### Simulate acoustic observations (~15 mins)
# (This is very slow for backward compatibility/historical reasons)
tic()
pairs <- CJ(a = seq_len(n_array_realisations),
            p = seq_len(n_path_realisations))
pairs$key <- paste(pairs$a, pairs$p)
acoustics <-
  pbapply::pblapply(seq_len(n_systems), function(i) {

    # Define detection probability parameters
    d <- detection_pars[i, , drop = FALSE]

    # Define path realisations in Julia
    julia_assign("paths", paths[[i]])
    julia_code(
      '
      # Format paths as matrix
      np = length(unique(paths.path_id))
      nt = length(unique(paths.timestep))
      xout = Matrix{StateXY}(undef, np, nt);
      for i in 1:nrow(paths)
        xout[paths.path_id[i], paths.timestep[i]] =
          StateXY(paths.map_value[i], paths.x[i], paths.y[i])
      end
      # Define paths object
      paths = xout
      '
    )

    # Simulate observations for each array
    pbapply::pblapply(seq_len(n_array), function(j) {

      print(glue("{i} / {n_systems}; {j} / {n_array}"))

      # Define observation model parameter(s)
      array <-
        arrays[[j]] |>
        select("sensor_id" = "receiver_id",
               "receiver_x", "receiver_y") |>
        mutate(receiver_alpha = d$alpha,
               receiver_beta = d$beta,
               receiver_gamma = d$gamma) |>
        as.data.table()

      # Simulate observations
      # > This returns a list with two elements:
      # * ModelObsDepthUniform (used for checking only)
      # * ModelObsAcousticLogisTrunc element
      # > That itself is a list, with one element for each path realisation
      sim <- sim_observations(.timeline = period,
                              .model_obs = c("ModelObsAcousticLogisTrunc", "ModelObsDepthUniform"),
                              .model_obs_pars = list(array,
                                                     data.table(sensor_id = 1L,
                                                                depth_shallow_eps = 0L,
                                                                depth_deep_eps = 0L))
                              )

      # Check that the list elements are correctly ordered by path ID (yes)
      check <- sim$ModelObsDepthUniform
      for (k in seq_len(length(check))) {
        check[[k]][, path_id := k]
      }
      check <- rbindlist(check)
      stopifnot(all.equal(paths[[i]]$map_value, check$obs))

      # Collate simulated acoustic observations (0, 1)
      # > We want one data.table with observations for all array/path realisations
      sim <- sim$ModelObsAcousticLogisTrunc
      for (k in seq_len(length(sim))) {
        sim[[k]][, array_id := arrays[[j]]$array_id[1]]
        sim[[k]][, path_id := k]
      }
      sim <- rbindlist(sim)
      # Check for simulations without detections
      sim$key <- paste(sim$array_id, sim$path_id)
      if (!all(pairs$key %in% sim$key[sim$obs > 0L])) {
        print(pairs[!(pairs$key %in% sim$key), ])
        warning("Some keys did not generate detections.",
                call. = FALSE, immediate. = TRUE)
      }
      sim
    })
  })
toc()

#### Drop realisations without sufficient detections
# Not all array/path realisations may generate detections
# We exclude simulations without sufficient observations
# (see below).

#### Save acoustic observations (0, 1) and detections (0, 1)
# Save acoustics
head(acoustics[[1]][[1]])
qs::qsave(acoustics, here_input("acoustics.qs"))
# Save detections
detections <- lapply(acoustics, function(elm) {
  lapply(elm, function(d) {
    d[obs == 1L, ]
  })
})
head(detections[[1]][[1]])
qs::qsave(detections, here_input("detections.qs"))


#########################
#########################
#### Define parameter sets

#### Overview
# For each array/path realisation, we will implement the algorithms
# ... using correct & incorrect parameter values for:
# * detection parameters (alpha, beta, gamma)
# * movement parameters  (mobility, shape, scale)
# * & sufficient/approximate values for algorithm controls (time step lengths, n particles)
# To do this, we need to build a dataframe with correct & incorrect (or approx) values
# * We hold each parameter at the right value & change others
# * We change some combinations of parameters simultaneously
# Specific research questions:
# (1) What happens when we under/overestimate gamma/mobility?
# ... I.e., can we confirm the results from Lavender et al. (2023)
# ... and examine what happens when this happens in combination?
# (2) What happens when we under/overestimate beta/shape?
# Note that we need to minimise the number of simulations.

#### Define algorithm parameters
# This returns a list
# * One element for each set of detection parameters
#   * One element for each set of path parameters
#     * One data.table with all parameter sets
alg_controls <- data.table(step = 2, n_particles = 10000)
alg_pars <-
  lapply(seq_len(n_systems), function(i){

      #### Define true parameters
      dpar <- detection_pars[i, , drop = FALSE]
      ppar <- path_pars[i, , drop = FALSE]
      pars <- cbind(dpar, ppar[])
      # Add algorithm control defaults
      pars$step        <- alg_controls$step # in mins
      pars$n_particles <- alg_controls$n_particles
      pars$flag        <- "defaults"

      #### Change selected parameters while holding others constant
      # * To minimise the number of simulations,
      # ... we (optionally) select a subset of parameters
      selected_pars <- c("alpha", "beta", "gamma",
                         "mobility", "shape", "scale",
                         "step", "n_particles")
      constants <-
        lapply(selected_pars, function(x) {
          # x <- selected_pars[1]

          # Define a sequence of parameter values
          vals <- pars[[x]] * c(0.1, 0.5, 1, 1.5, 2)
          # For delta t & n_particles, we use selected hard-coded defaults
          if (x == "step") {
            vals <- c(2, 4, 8)
          }
          if (x == "n_particles") {
            vals <- c(5000, 10000, 20000)
          }
          dp <- data.table(x = vals)
          colnames(dp) <- x
          keep <- colnames(pars)[!(colnames(pars) %in% x)]
          dp <-
            dp |>
            cbind(pars[, ..keep]) |>
            distinct() |>
            dplyr::select(colnames(pars)) |>
            as.data.table()

          # For step, we need to update the 'correct' step length parameters
          # (to account for the longer duration)
          if (x == "step") {
            for (v in seq_len(length(vals))) {
              dp$mobility[v] <- dp$mobility[v] * (dp$step[v] / alg_controls$step)
              dp$shape[v]    <- dp$shape[v] * (dp$step[v] / alg_controls$step)
            }
          }

          # Return dp
          dp

        }) |> rbindlist()
      constants$flag <- "constants"

      #### Change gamma and mobility simultaneously
      # (while holding other parameters constant)
      gm <- CJ(gamma = unique(constants$gamma),
               mobility = unique(constants$mobility))
      keep <- colnames(pars)[!(colnames(pars) %in% c("gamma", "mobility"))]
      gm <- cbind(gm, pars[, ..keep])
      gm$flag <- "gm"

      #### Change beta and shape simultaneously
      # (while holding other parameters constant)
      bs <- CJ(beta = unique(constants$beta),
                        shape = unique(constants$shape))
      keep <- colnames(pars)[!(colnames(pars) %in% c("beta", "shape"))]
      bs <- cbind(bs, pars[, ..keep])
      bs$flag  <- "bs"

      #### Define all unique parameter sets
      sets <-
        rbind(pars, constants, gm, bs) |>
        distinct() |>
        mutate(alg_par = row_number(),
               system_type = i,
               path_type = i) |>
        data.table()

      #### Return parameter sets
      attr(sets, "true_pars") <- pars
      sets
  })
# Save algorithm parameters
n_alg_pars <- nrow(alg_pars[[1]])
saveRDS(alg_pars, here_input("alg_pars.rds"))


#########################
#########################
#### Collate simulations

#### Collect simulations
sims <-
  CJ(
  combination = seq_len(n_systems),
  array_type = seq_len(n_array),
  array_realisation = seq_len(n_array_realisations),
  path_type = NA_integer_,
  path_realisation = seq_len(n_path_realisations),
  alg_par = seq_len(n_alg_pars)
  ) |>
  arrange(combination, array_type, array_realisation,
          path_type, path_realisation, alg_par) |>
  as.data.table()
sims$system_type <- sims$combination
sims$path_type   <- sims$combination

#### Add array/path/parameter information
# Add array information
array_pars_dt <- lapply(seq_len(length(arrays)), \(i) {
  arrays[[i]] |>
    slice(1L) |>
    mutate(array_type = i,
           array_realisation = array_id) |>
    select(array_type, array_realisation, arrangement, n_receiver) |>
    as.data.table()
}) |>
  rbindlist()
nrow(sims)
sims <- merge(sims, array_pars_dt, by = c("array_type", "array_realisation"))
nrow(sims)
# Add algorithm parameters
alg_pars_dt <-
  alg_pars |>
  rbindlist() |>
  select(system_type, path_type, alg_par,
         alpha, beta, gamma,
         mobility, shape, scale,
         step, n_particles,
         flag) |>
  as.data.table()
nrow(sims)
sims <- merge(sims, alg_pars_dt,
              by = c("system_type", "path_type", "alg_par"))
nrow(sims)
# Drop duplicate parameter combinations
nrow(sims)
sims <-
  sims |>
  group_by(path_type, array_type, array_realisation, path_realisation,
           alpha, beta, gamma,
           mobility, shape, scale,
           step, n_particles) |>
  slice(1L) |>
  as.data.table()
nrow(sims)
# Define 'performance' (TRUE/FALSE)
# * 'performance' indicates a simulation using the correct parameters
# * (these simulations are used to compare algorithm performance)
sims$performance <- FALSE
true_pars <- cbind(detection_pars, path_pars, alg_controls)
for (i in seq_len(nrow(true_pars))) {
  p <- true_pars[i, ]
  pos <- which(
    sims$alpha == p$alpha &
    sims$beta == p$beta &
    sims$gamma == p$gamma &
    sims$mobility == p$mobility &
    sims$shape == p$shape &
    sims$scale == p$scale &
    sims$step == p$step &
    sims$n_particles == p$n_particles
  )
  sims$performance[pos] <- TRUE
}
table(sims$performance)
stopifnot(length(which(sims$performance)) ==
            n_array * n_array_realisations * n_path * n_path_realisations)

#### Subset sensitivity combinations analysis
# * Select 10 array designs for this analysis
# * This helps to minimise the number of simulations
# * (And we can't plot more results easily anyway)
nr <- sort(unique(array_pars$number))
nr <- nr[seq(1, length(nr), by = 10)]
pos <- which(!sims$performance & sims$flag == "gm" & !(sims$n_receiver %in% nr))
if (length(pos) > 0) {
  sims <- sims[-pos, ]
}
pos <- which(!sims$performance & sims$flag == "bs" & !(sims$n_receiver %in% nr))
if (length(pos) > 0) {
  sims <- sims[-pos, ]
}

#### Exclude simulations with insufficient detections (~37 s)
# Count the number of detections for each simulation
tic()
n <- nrow(sims)
sims$count <- 0L
for (i in seq_len(n)) {
  svMisc::progress(i, n)
  sims$count[i] <-
    detections[[sims$system_type[i]]][[sims$array_type[i]]][
    array_id == sims$array_realisation[i] & path_id == sims$path_realisation[i] & obs == 1L, ] |>
    fnrow()
}
toc()
hist(sims$count, breaks = 100)
# Drop simulations with insufficient detections
n; nrow(sims[sims$count > 5L, ])
sims <- sims[sims$count > 5L, ]
range(sims$count)

#### Define IDs
sims$id <- seq_len(nrow(sims))
head(sims)
nrow(sims)

#### Save sims
length(which(sims$performance)) # 1193
nrow(sims)                      # 39516
saveRDS(sims, here_input("sims.rds"))


toc()
# 18 mins

#### End of code.
#########################
#########################
