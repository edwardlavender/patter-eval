#########################
#########################
#### run-patter-debug.R

#### Aims
# 1) Resolve convergence failures for ACDCPF algorithm
# * All performance simulations in run-patter.R should convergence
# * Some ACDCPF algorithm implementations did not converge
# * This script examines convergence failures & confirms that convergence
# * ... is achieved with larger numbers of particles

#### Prerequisites
# 1) run-patter.R


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
library(prettyGraphics)
library(tictoc)
dv::src()

#### Load data
sims   <- readRDS(here_input("sims-performance.rds"))
spat    <- terra::rast(here_data("sims", "input", "spat.tif"))
success <- readRDS(here_data("sims", "output", "success", "patter-performance.rds"))


#########################
#########################
#### Examine success

#### Count convergence failures
table(success$acpf)
table(success$acdcpf)

#### Identify convergence failures
failures <- success[acdcpf == FALSE, ]
sims <- sims[row %in% failures$row, ]
nrow(sims)


#########################
#########################
#### Re-run forward run

spat <- terra::unwrap(terra::wrap(spat))
terra::inMemory(spat)

tic()
convergence <-
  pbapply::pblapply(
    X = split(sims, seq_len(nrow(sims))),
    cl = NULL, # floor(nrow(sims) / 2L),
    FUN = function(sim) {

      #### Define simulation
      # sim   <- sims[1, ]
      print(sim$row)
      out_file <- here_data("sims", "output", "debug", "patter", "convergence", paste0(sim$row, ".rds"))
      if (file.exists(out_file)) {
        return(readRDS(out_file))
      }

      #### Set up data list and observations
      # (This code is copied from workflow_patter())
      dlist <- read_dlist(sim)
      dlist$spatial$bathy                <- spat
      dlist$algorithm$detection_overlaps <- read_overlaps(sim)
      dlist$algorithm$detection_kernels  <- read_kernels(sim)
      # Define observation timeline
      obs <- pf_setup_obs(.dlist = dlist,
                          .trim = FALSE,
                          .step = paste(sim$step, "mins"),
                          .mobility = sim$mobility,
                          .receiver_range = dlist$data$moorings$receiver_range[1])
      # Adjust mobility
      # * obs$mobility is used in pf_rpropose_reachable()
      # * We adjust mobility here to account for the discretisation error
      obs[, mobility := mobility + sr]
      # Include depths
      obs[, depth_shallow := obs$depth - 5]
      obs[, depth_deep := obs$depth + 5]

      #### Define arguments for forward run
      # (This code is copied from get_ud_patter())
      rargs <- list(.shape = sim$shape, .scale = sim$scale,
                    .mobility = sim$mobility)
      dargs <- list(.shape = sim$shape, .scale = sim$scale,
                    .mobility = sim$mobility + sr)
      # Likelihood functions
      lik <-  list(pf_lik_dc = pf_lik_dc,
                   acs_filter_container = acs_filter_container,
                   pf_lik_ac = pf_lik_ac)
      # Record opts
      record <- pf_opt_record(.save = TRUE)

      #### Run forward simulation
      # (We will try multiple values for n to achieve convergence)
      # Define output list
      out <- list()
      # Define starting 'success'
      s     <- FALSE
      ns    <- c(1e3L, 1e4L, 1e5L)
      count <- 1L
      # Run forward simulation until convergence is achieved
      while (isFALSE(s) && count <= length(ns)) {
        # Run simulation
        print(count)
        n <- ns[count]
        print(n)
        tic()
        ssf()
        out_pff  <- pf_forward(.obs = obs,
                               .dlist = dlist,
                               .rargs = rargs,
                               .dargs = dargs,
                               .likelihood = lik,
                               .n = n,
                               .record = record,
                               .trial = pf_opt_trial(.trial_sampler = FALSE),
                               # .control = pf_opt_control(.sampler_batch_size = 1000L),
                               .verbose = here_data("sims", "output", "log", "patter", "performance-debug",
                                                    paste0(sim$row, "-", n, ".txt")))
        toc()
        # Record success
        s <- out_pff$convergence
        # Update output list
        out[[count]] <- data.frame(row = sim$row,
                                   n = n,
                                   success = s)
        # Update loop controls
        count <- count + 1L
      }

      #### Return outputs
      out <- rbindlist(out)
      saveRDS(out, out_file)
      out
    })
toc()


#########################
#########################
#### Examine results

cdt <- rbindlist(convergence)
# View(cdt)

cdt |>
  group_by(row) |>
  arrange(n, .by_group = TRUE) |>
  slice_tail(n = 1L) |>
  as.data.table()


#########################
#########################
#### Debugging

#### Define simulation
# * row 602 struggles to converge (position ~649)
sim <- sims[row == 602, ]

#### Read data
path     <- read_path(sim)
moorings <- read_array(sim)
acc      <- read_acoustics(sim)
dlist    <- read_dlist(sim)

#### Plot area & path
terra::plot(spat)
add_sp_path(path$x, path$y, length = 0.01)
text(moorings$receiver_easting, moorings$receiver_northing, moorings$receiver_id, col = "red", font = 2)

#### Plot acoustic and archival observations
acc
# plot(acc$timestep, rep(1, nrow(acc)))
plot(dlist$data$archival$timestep,
     dlist$data$archival$depth * -1, type = "l",
     ylim = c(-350, 0))

points(dlist$data$acoustics$timestep, rep(1, nrow(acc)))

#### Run the forward simulation for this individual (~52 s)
# Define arguments (as above)
dlist <- read_dlist(sim)
dlist$spatial$bathy                <- spat
dlist$algorithm$detection_overlaps <- read_overlaps(sim)
dlist$algorithm$detection_kernels  <- read_kernels(sim)
obs <- pf_setup_obs(.dlist = dlist,
                    .trim = FALSE,
                    .step = paste(sim$step, "mins"),
                    .mobility = sim$mobility,
                    .receiver_range = dlist$data$moorings$receiver_range[1])
obs[, mobility := mobility + sr]
obs[, depth_shallow := obs$depth - 5]
obs[, depth_deep := obs$depth + 5]
rargs <- list(.shape = sim$shape, .scale = sim$scale,
              .mobility = sim$mobility)
dargs <- list(.shape = sim$shape, .scale = sim$scale,
              .mobility = sim$mobility + sr)
lik <-  list(pf_lik_dc = pf_lik_dc,
             acs_filter_container = acs_filter_container,
             pf_lik_ac = pf_lik_ac)
record <- pf_opt_record(.save = TRUE)
# Run forward simulation
# * 1e4 particles: ~41 s (fail)
# * 1e6 particles: ~1.5 mins (fail)
# * 1e6 particles: ~5.5 mins (fail)
tic()
ssf()
out_pff  <- pf_forward(.obs = obs,
                       .dlist = dlist,
                       .rargs = rargs,
                       .dargs = dargs,
                       .likelihood = lik,
                       .n = 1e6L,
                       .record = pf_opt_record(.sink = tempdir()),
                       .trial = pf_opt_trial(.trial_sampler = FALSE),
                       .verbose = TRUE)

toc()

#### Examine diagnostics
# TO DO

#### Plot particle samples
# (1e3 particles: ~55s, 20 cl)
# (1e6 particles: 3 mins, 50 cl)
# At time step 554, the last particles in the vicinity of the true position die
# This may be partly related to the true position being near the boundaries of the study area
# Importantly, once the right patch of particles dies, it is very difficult for the algorithm to recover
# B/c there are no particles in the right vicinity
# in this instance the only patch is on the other side of the study area
tic()
pf_plot_history(.dlist = dlist,
                .steps = NULL,
                .forward = file.path(tempdir(), "history"),
                .png = list(filename = here_fig("debug", "602")),
                .add_forward = list(pch = 21, col = "red", bg = "red", cex = 0.25),
                .add_layer = function(t) {
                  # Plot true path
                  add_sp_path(path$x, path$y, length = 0.01, lwd = 0.75, col = scales::alpha("dimgrey", 0.75))
                  # Add moorings
                  points(moorings$receiver_easting, moorings$receiver_northing,
                         col = "blue", lwd = 2)
                  # Plot true position at time t
                  points(path$x[t], path$y[t], col = "orange", lwd = 3)
                  terra::sbar(d = 500, xy = "bottomleft")
                },
                .cl = 50L)
toc()


#### End of code.
#########################
#########################
