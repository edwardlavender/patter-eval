#########################
#########################
#### run-patter-trials.R

#### Aims
# 1) Run patter time trials
# * The backward sampler is an expensive algorithm
# * Here, we compare the speeds of alternative implementations of the algorithm
# * This code supports, in particular, the large number of sensitivity
# * analyses we need to implement

#### Prerequisites
# 1) Develop {patter} workflow


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

#### Load data
spat <- terra::rast(here_input("spat.tif"))
sims <- readRDS(here_input("sims-performance.rds"))


#########################
#########################
#### Forward run

# ~54 s (siam-linux20)
tic()

#### Define simulation
sim <- sims[1, ]

#### Define data list
dlist <- read_dlist(sim)
dlist$spatial$bathy <- spat
dlist$algorithm$detection_overlaps <- read_overlaps(sim)
dlist$algorithm$detection_kernels  <- read_kernels(sim)

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

#### Define proposal functions
rargs <- list(.shape = sim$shape, .scale = sim$scale,
              .mobility = sim$mobility)
dargs <- list(.shape = sim$shape, .scale = sim$scale,
              .mobility = sim$mobility + sr)

#### Define likelihood functions
lik <-  list(pf_lik_dc = pf_lik_dc,
             acs_filter_container = acs_filter_container,
             pf_lik_ac = pf_lik_ac)

#### Define control arguments
record  <- pf_opt_record(.save = TRUE)
control <- pf_opt_control(.sampler_batch_size = 1000L)

#### Forward simulation
tic()
ssf()
out_pff  <- pf_forward(.obs = obs,
                       .dlist = dlist,
                       .rargs = rargs,
                       .dargs = dargs,
                       .likelihood = lik,
                       .n = sim$n_particles,
                       .control = control,
                       .record = record,
                       .verbose = FALSE)
toc()

toc()


#########################
#########################
#### Backward run comparisons

#### pf_backward_sampler_v(): 4m 21s (siam-linux20)
if (FALSE) {
  out_pfbs <-
    pf_backward_sampler_v(.history = out_pff$history,
                          .dpropose = pf_dpropose,
                          .obs = obs,
                          .dlist = dlist,
                          .dargs = dargs,
                          .record = record,
                          .verbose = FALSE)
  out_pfbs$time$duration
}

#### pf_backward_sampler_p():
if (FALSE) {
  dlist$algorithm$sim <- sim
  out_pfbs <- pf_backward_sampler_p(.history = out_pff$history,
                                    .dpropose = pf_dpropose_read,
                                    .obs = obs,
                                    .dlist = dlist,
                                    .dargs = list(),
                                    .record = record,
                                    .verbose = TRUE)
}

#### pf_backward_sampler_p_fst() with .read = TRUE:
if (FALSE) {
  out_pfbs <- pf_backward_sampler_p_fst(.history = out_pff$history,
                                        .dlist = dlist,
                                        .read = TRUE)
}

#### pf_backward_sampler_p_fst() with .read = FALSE:
if (FALSE) {
  out_pfbs <- pf_backward_sampler_p_fst(.history = out_pff$history,
                                        .dlist = dlist,
                                        .read = FALSE)
}

#### pf_backward_sampler_cpp()
# TO DO


#### End of code.
#########################
#########################
