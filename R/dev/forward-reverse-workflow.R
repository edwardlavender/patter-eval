###########################
###########################
#### forward-reverse-workflow.R

#### Aims
# 1) Test forward/reverse workflow for PF

#### Prerequisites
# 1) Define datasets


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
# rstudioapi::restartSession()
dv::clear()
op <- options(error = NULL)
# op <- options(error = recover)
# op <- options(error = function(...) beepr::beep(7))

#### Essential packages
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr)
library(prettyGraphics)
library(tictoc)

#### Load data
dv::src()
spat <- terra::rast(here_data("sims", "input", "spat.tif"))


###########################
###########################
#### Preparation

#### Define simulation
sims   <- readRDS(here_input("sims-performance.rds"))
sim    <- sims[row == 602, ]
# sim    <- sims[row == 897, ]
sink   <- here_data("sims", "output", "debug", "patter", "experiments", "forward-backward", sim$row)
dir.create(sink, recursive = TRUE)
dir.create(file.path(sink, "forward", "fwd"), recursive = TRUE)
dir.create(file.path(sink, "forward", "bwd"), recursive = TRUE)
dir.create(file.path(sink, "forward", "full"), recursive = TRUE)
dir.create(file.path(sink, "backward"), recursive = TRUE)

#### Define data list
dlist <- read_dlist(sim)
dlist$spatial$bathy <- spat
dlist$algorithm$detection_overlaps <- read_overlaps(sim)
dlist$algorithm$detection_kernels  <- read_kernels(sim)
# (optional) Define origin
# * Not currently implemented

#### Define observation timeline
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

#### Proposal functions
rargs <- list(.shape = sim$shape, .scale = sim$scale,
              .mobility = sim$mobility)
dargs <- list(.shape = sim$shape, .scale = sim$scale,
              .mobility = sim$mobility + sr)

#### Likelihood functions
lik <-  list(pf_lik_dc = pf_lik_dc,
             acs_filter_container = acs_filter_container,
             pf_lik_ac = pf_lik_ac)

#### Control options
n       <- 1000L
trial   <- pf_opt_trial(.trial_sampler = FALSE,
                        .trial_revert = 0L,
                        .trial_resample_crit = 0L)
record  <- pf_opt_record(.save = TRUE,
                         .sink = file.path(sink, "forward", "fwd"))
verbose <- file.path(sink, "forward", "fwd", "log-1.txt")


###########################
###########################
#### Run forwards

#### Define argument list
args <- list(.obs = copy(obs),
             .dlist = dlist,
             .rargs = rargs,
             .dargs = dargs,
             .likelihood = lik,
             .n = n,
             .trial = trial,
             .record = record,
             .verbose = verbose)

#### Initiate parameters
pos_detections <- which(obs$detection == 1L)
continue       <- TRUE
breaks         <- list()
count          <- 1L

tic()
while (continue) {

  #### User messages
  cat("\n")
  cat(paste0(rep("-", 25), collapse = "-"))
  message(args$.obs$timestep[1])

  #### Run simulation
  tic()
  ssf()
  out_pff <- do.call(patter::pf_forward, args)
  toc()

  #### Identify the time step we have reached
  # (This code assumes outputs are saved in memory)
  convergence <- out_pff$convergence
  timestep    <- out_pff$history[[length(out_pff$history)]]$timestep[1]
  # Record the time step at which the algorithm failed
  # * We will use this later to identify surrounding anchor points
  # * & rerun the algorithm backwards
  if (isFALSE(convergence)) {
    breaks[[count]] <- timestep
  }

  #### Define whether or not to continue the algorithm
  # * We continue if convergence = FALSE & we have not passed the last anchor point
  continue <- isFALSE(convergence) && timestep < max(pos_detections)
  if (continue) {

    count <- count + 1L

    #### Jump to the next anchor point
    jump <- pos_detections[which.min(abs(pos_detections - timestep))]
    cat(glue::glue("Jumping to {jump}...\n"))

    #### Update arguments for next run
    # Update observations
    args$.obs <-
      obs |>
      copy() |>
      slice(jump:n()) |>
      as.data.table()
    stopifnot(args$.obs$detection[1] == 1L)
    # Update origin (optional)
    args$.dlist$spatial$origin <- NULL
    # Update verbose
    args$.verbose     <- file.path(sink, "forward", "fwd", paste0("log-", count, ".txt"))

  }

}
toc()


###########################
###########################
#### Run backwards

#### Update output sink
record  <- pf_opt_record(.save = TRUE,
                         .sink = file.path(sink, "forward", "bwd"))
args$.record <- record

#### Identify anchor points
# Define function
x <- 5
v <- c(1, 2, 4, 6)
interval <- function(x, v) {
  stopifnot(!is.unsorted(v))
  index <- which(v > x)[1]
  v[c(index - 1L, index)]
}
interval(x, v)
# Identify points
index <- c(1, pos_detections,  nrow(obs)) |> unique() |> sort()
anchors <-
  lapply(breaks, \(b) {
    interval(b, index)
  })

#### Run algorithm backwards between anchor points
tic()
pbapply::pblapply(anchors, function(anchor) {

  cat(paste0(rep("-", 25), collapse = "-"))
  # anchor = anchors[[1]]
  # anchor = c(1, 655)
  message(anchor[1])

  # Identify observations
  obs_bwd <-
    obs |>
    copy() |>
    slice(min(anchor):max(anchor))

  # Reverse observations
  obs_bwd <-
    obs_bwd |>
    arrange(desc(timestep)) |>
    mutate(detection_id = as.integer(factor(detection_id, unique(detection_id))),
           receiver_id_next = patter:::.pf_setup_obs_receiver_id_next(receiver_id),
           receiver_id_next_key = lapply(receiver_id_next, acs_setup_receiver_key) |> unlist()) |>
    mutate(
      step_forwards = row_number(),
      step_backwards = rev(.data$step_forwards),
      buffer_future = mobility * .data$step_backwards - 500,
      buffer_future_incl_gamma = .data$buffer_future + 750) |>
    as.data.table()

  # Update origin
  # (optional)
  dlist$spatial$origin <- NULL

  # Update argument list
  args$.obs         <- obs_bwd
  args$.dlist       <- dlist
  args$.verbose     <- file.path(sink, "forward", "bwd", paste0("log-", anchor[1], ".txt"))

  # Run algorithm backwards between anchor points
  tic()
  ssf()
  out_pff <- do.call(patter::pf_forward, args)
  toc()

  # Return outputs
  out_pff$convergence

})
toc()


###########################
###########################
#### Assembly

#### Join outputs from forward and backward run
# * Loop over time steps
# * Read fwd & bwd files (if exist)
# * rbind
# * Write to full/ directory
np <-
  lapply(obs$timestep, function(t) {

  print(t)
  fdt <- bdt <- NULL
  fwd <- file.path(sink, "forward", "fwd", "history", paste0(t, ".parquet"))
  bwd <- file.path(sink, "forward", "bwd", "history", paste0(t, ".parquet"))

  if (file.exists(fwd)) {
    fdt <- arrow::read_parquet(fwd)
    stopifnot(all(fdt$timestep %in% t))
  }

  if (file.exists(bwd)) {
    bdt <- arrow::read_parquet(bwd)
    stopifnot(all(bdt$timestep %in% t))
  }

  dt <-
    list(fdt, bdt) |>
    plyr::compact() |>
    rbindlist()
  stopifnot(nrow(dt) > 0L)

  arrow::write_parquet(dt, file.path(sink, "forward", "full", paste0(t, ".parquet")))

  nrow(dt)

}) |> invisible()
# Samples contain different numbers of particles
sort(table(unlist(np)))

#### Use backward sampler for assembly
# * Read files from full/
# * Run backward sampler in the usual way
# * (UPDATE TO HANDLE VARYING N and weights)
tic()
out_pfb <- pf_backward_sampler_v(.history = file.path(sink, "forward", "full"),
                                 .dlist = dlist,
                                 .dargs = dargs,
                                 .record = pf_opt_record(.save = TRUE))
toc()

table(sapply(out_pfb$history, nrow))

#### Plot particle samples
tic()
path     <- read_path(sim)
moorings <- read_array(sim)
sink_fig <- here_fig("debug", "experiments", sim$row)
dir.create(sink_fig)
dlist$spatial$bathy <- terra::unwrap(terra::wrap(spat))
pf_plot_history(.dlist = dlist,
                .steps = NULL,
                .forward = out_pfb,
                .png = list(filename = sink_fig),
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
###########################
###########################
