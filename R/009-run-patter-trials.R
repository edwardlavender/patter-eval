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
obs[, depth_shallow_eps := 5]
obs[, depth_deep_eps := 5]

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
                         .trial = pf_opt_trial(.trial_resample_crit = Inf),
                         .control = control,
                         .record = record,
                         .verbose = TRUE)

toc()


#########################
#########################
#### Backward run comparisons

#### pf_backward_sampler_v(): 4m 21s (siam-linux20)
if (TRUE) {
  t0_start <- Sys.time()
  ssf()
  out_0 <-
    pf_backward_sampler_v(.history = out_pff$history,
                          .dpropose = pf_dpropose,
                          .obs = obs,
                          .dlist = dlist,
                          .dargs = dargs,
                          .record = record,
                          .verbose = FALSE)
  t0_end <- Sys.time()
  difftime(t0_end, t0_start)

  # Check diagnostics]
  library(prettyGraphics)
  (diag <- pf_diag_summary(out_0$history))
  p <- 1:250
  pretty_plot(diag$timestep[p], diag$nu[p],
              ylim = c(0, NA), type = "b", cex = 0.5)
  pretty_plot(diag$timestep[p], diag$ess[p],
              ylim = c(0, NA), type = "b", cex = 0.5,
              xlab = "Time (step)", ylab = "ESS")
  abline(v = obs$timestep[obs$detection == 1L], col = "red", lty = 3)

}

#### pf_backward_sampler_p():
if (TRUE) {
  dlist$algorithm$sim <- sim
  t1_start <- Sys.time()
  ssf()
  out_1 <- pf_backward_sampler_p(.history = out_pff$history,
                                 .dpropose = pf_dpropose_read,
                                 .obs = obs,
                                 .dlist = dlist,
                                 .dargs = list(),
                                 .record = record,
                                 .verbose = FALSE)
  t1_end <- Sys.time()
  difftime(t1_end, t1_start)

}

#### pf_backward_sampler_p_fst() with .read = TRUE:
if (TRUE) {
  t2_start <- Sys.time()
  ssf()
  out_2 <- pf_backward_sampler_p_fst(.history = out_pff$history,
                                     .dlist = dlist,
                                     .read = TRUE)
  t2_end <- Sys.time()
  difftime(t2_end, t2_start)
}

#### pf_backward_sampler_p_fst() with .read = FALSE:
if (TRUE) {
  t3_start <- Sys.time()
  ssf()
  out_3 <- pf_backward_sampler_p_fst(.history = out_pff$history,
                                     .dlist = dlist,
                                     .read = FALSE)
  t3_end <- Sys.time()
  difftime(t3_end, t3_start)
}

#### pf_backward_sampler_cpp()
if (TRUE) {
  # Compile C++ code
  Rcpp::sourceCpp(here_src("cpp", "pf_backward_sampler_cpp.cpp"))
  # Define list of coordinate matrices
  pxy <- lapply(out_pff$history, function(elm) {
    as.matrix(elm[, .(x_now, y_now)])
  })
  # Define path matrix
  # * One column per particle
  # * Two columns per time step
  t4_start <- Sys.time()
  ssf()
  out_4 <- pf_backward_sampler_cpp(particles = pxy,
                                   shape = sim$shape,
                                   scale = sim$scale,
                                   mobility = sim$mobility)
  t4_end <- Sys.time()
  difftime(t4_end, t4_start)
}


#### End of code.
#########################
#########################
