#########################
#########################
#### run-patter.R

#### Aims
# 1) Estimate UDs using patter (for performance & sensitivity simulations)

#### Prerequisites
# 1) Simulate data & prepare for algorithm implementation


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
library(microbenchmark)
library(parallel)
library(tictoc)
dv::src()

#### Load data
spatw      <- readRDS(here_input("spatw.rds"))
win        <- qs::qread(here_input("win.qs"))
sims_for_performance    <- readRDS(here_input("sims-performance.rds"))
sims_for_performance_ls <- split(sims_for_performance, sims_for_performance$id)


#########################
#########################
#### Julia connection

#### Define sims
type <- "performance"
sims <- sims_for_performance
sims[, n_particles := 5e4L]
# type <- "sensitivity"
# sims <- sims_for_sensitivity

#### Define multi-threading option
# We can multi-thread R or Julia code
# Multi-threading in R uses a socket cluster (forking crashes the R session)
multithread <- c("R", "Julia")
multithread <- multithread[1]

#### Connect to Julia
if (multithread == "R") {

  rsockets <- 10L

  setDTthreads(threads = 1)
  # julia_connect(.threads = 1)

} else {

  rsockets <- 1L

  setDTthreads(threads = Sys.getenv("OMP_NUM_THREADS"))
  julia_connect()

  set_seed()
  set_map(terra::unwrap(spatw))

}


#########################
#########################
#### Time trials


#########################
#### UD estimation

#### Results

# This code compares the speed of UD estimation via map_dens() with different options

# `.discretise  = TRUE` (with sigma = bw.diggle()):
# > There is a major speed benefit for `.discretise = TRUE`:
# > In `.map_coord.dt()`, `.discretise = TRUE` is approx. twice as fast (e.g., 5 s verus 10 s);
# > In map_dens(), total time is three times faster (12 s versus 71 s, with bw.diggle())

# > Using a coarser UD grid is an option for faster estimation
# > `sigma` (with `discretise = TRUE`):
# * bw.diggle(): 12 s (tends to undersmooth, uses cross validation)
# * bw.scott(): 1 s (tends to undersmooth, uses a rule of thumb)

# File size
# * We require 23 MB to store particles from one run of the forward filter:
# * 23.7 MB
# * Therefore, we require > 2000 GB (2 TB) to store particles from all filter runs!
# * 2844 GB (23.7 * 2 *  2 * 30000 / 1e3)
# * 23.7 * 2 (forward, smoothed particles) * 2 (ACPF, ACDCPF) * 30000 simulations

if (FALSE && multithread == "Julia") {

  #### Define simulation
  # (TO DO) Repeat the code with a few random choices to check result consistency
  sim       <- sims_for_performance[1, ]

  #### Read data
  acoustics <- read_acoustics(sim)
  archival  <- read_archival(sim)

  #### Algorithm components
  spat       <- terra::unwrap(spatw)
  t1         <- min(acoustics$timestamp)
  tT         <- max(acoustics$timestamp)
  timeline   <- seq(t1, tT, by = "2 mins")
  model_move <-
    move_xy(dbn_length =
              glue::glue("truncated(Gamma({sim$shape},
                                    {sim$scale}),
                              upper = {sim$mobility})"),
            dbn_angle = "Uniform(-pi, pi)")

  #### Filter args
  args <- list(.map = spat,
               .timeline = timeline,
               .state = "StateXY",
               .xinit = NULL,
               .xinit_pars = list(mobility = sim$mobility),
               .yobs = list(acoustics),
               .model_obs = "ModelObsAcousticLogisTrunc",
               .model_move = model_move,
               .n_particle = sim$n_particles,
               .n_resample = sim$n_particles,
               .verbose = FALSE)

  #### Forward filter
  # Run filter
  out_pff  <- do.call(pf_filter, args, quote = TRUE)
  # Check file size
  tmp      <- tempfile(fileext = ".qs")
  qs::qsave(out_pff$states, tmp)
  file.size(tmp) / 1e6

  #### Backward filter
  args$.direction <- "backward"
  out_pfb  <- do.call(pf_filter, args, quote = TRUE)

  #### Backward smoother
  out_smo <- pf_smoother_two_filter(.map = spat,
                                    .mobility = sim$mobility,
                                    .n_particle = 1000L)

  #### patter:::.map_coord.dt() speed (.discretise)
  # Use particles from forward filter
  tic()
  bm_pff_coord <-
    microbenchmark(
    continuous = patter:::.map_coord.dt(.map = spat, .coord = out_pff$states, .discretise = FALSE) |> nrow(),
    discrete = patter:::.map_coord.dt(.map = spat, .coord = out_pff$states, .discretise = TRUE) |> nrow(),
    times = 10L
  )
  toc()
  # Use particles from particle smoother (to check consistency of results)
  tic()
  bm_smo_coord <-
    microbenchmark(
      continuous = patter:::.map_coord.dt(.map = spat, .coord = out_smo$states, .discretise = FALSE) |> nrow(),
      discrete = patter:::.map_coord.dt(.map = spat, .coord = out_smo$states, .discretise = TRUE) |> nrow(),
      times = 10L
    )
  toc()

  #### UD estimation speed (.discretise)

  # `.discretise = FALSE`
  map_dens(.map = terra::unwrap(spatw),
           .owin = win,
           .coord = out_pff$states,
           .discretise = FALSE,
           .plot = FALSE,
           .verbose = TRUE,
           sigma = spatstat.explore::bw.diggle)$ud |> terra::plot()
  # `patter::map_dens()` called @ 2024-06-18 14:43:53...
  # ... 14:43:53: Processing `.map`...
  # ... 14:43:53: Building XYM...
  # ... 14:44:03: Defining `ppp` object...
  # ... 14:44:03: Estimating density surface...
  # ... 14:45:14: Scaling density surface...
  # ... 14:45:14: * Resampling density surface onto `.map`...
  # `patter::map_dens()` call ended @ 2024-06-18 14:45:14 (duration: ~1.35 min(s)).

  # `.discretise = TRUE`
  map_dens(.map = terra::unwrap(spatw),
           .owin = win,
           .coord = out_pff$states,
           .discretise = TRUE,
           .plot = FALSE,
           .verbose = TRUE,
           sigma = spatstat.explore::bw.diggle)$ud |> terra::plot()
  # `patter::map_dens()` called @ 2024-06-18 14:45:32...
  # ... 14:45:32: Processing `.map`...
  # ... 14:45:32: Building XYM...
  # ... 14:45:37: Defining `ppp` object...
  # ... 14:45:37: Estimating density surface...
  # ... 14:45:59: Scaling density surface...
  # ... 14:45:59: * Resampling density surface onto `.map`...
  # `patter::map_dens()` call ended @ 2024-06-18 14:45:59 (duration: ~27 sec(s)).

  #### UD estimation speed (sigma)
  map_dens(.map = terra::unwrap(spatw),
           .owin = win,
           .coord = out_pff$states,
           .discretise = TRUE,
           .plot = FALSE,
           .verbose = TRUE,
           sigma = spatstat.explore::bw.scott)$ud |> terra::plot()

}


#########################
#### Algorithm runs

# This code tests whether we should multi-thread R or Julia
# We run the algorithms:
# * Using R socket clusters, with Julia in single-threaded mode
# * Using R serial code, with Julia in multi-threaded mode
# * All tests run on Apple M2 MacBook with 10 rsockets or 10 Julia threads
# * Another option is to split algorithm runs from UD estimation
# * (and run the former in multi-threaded Julia and the latter in R,
# * but this is not currently tested)

if (FALSE) {

  # Define simulation test subset
  # > Now run the code in the following section.
  sims <- sims[1:100L, ]

  #### Guesses

  # Estimated multi-threaded R, using rsockets CPUs (bw.diggle):
  guess <- 158
  nsim  <- nrow(sims)
  (nsim * guess)/60/60/rsockets    # hours
  (nsim * guess)/60/60/24/rsockets # days

  # Estimated multi-threaded Julia, using rsockets threads:
  # > Even in the best case, we have to pay 60 * 2 s to estimate the UDs
  guess <- 120
  (nsim * guess)/60/60/1    # hours
  (nsim * guess)/60/60/24/1

  #### Empirical timings

  # 1 algorithm run, 1 thread (bw.diggle):
  # * Forward filter:    4 s
  # * Backward filter:   4 s
  # * Smoother:         90 s
  # * UD estimation:    60 s (2 * 30 for the forward & smoothed UDs)
  #                    158 s (* 2 = 316 s for ACPF & ACDCPF)

  # 1 algorithm run, 1 thread (bw.scott):
  # (not currently implemented)

  # Multi-threaded R (bw.diggle):
  # 100 simulations: 4854.023 s (1.34 hours)
  # Total ETA      : 16 days

  # Multi-threaded R (bw.scott):
  # (not currently implemented)

  # Multi-threaded Julia (bw.diggle):
  # > 2 hours

  # Multi-threaded Julia (bw.scott):
  # (not currently implemented)

}


#########################
#########################
#### Implementation

#### Initialise cluster chunks (log files)
# Define chunks to iterate over in parallel
chunks  <- patter:::cl_chunks(rsockets, nrow(sims))
nchunks <- length(chunks)
# Create a log file for each chunk
lapply(seq_len(nchunks), function(i) {
  file.create(here_output("log", "patter", type, paste0("log-", i, ".txt")))
}) |> invisible()

#### Initialise cluster
# We use a socket cluster, as the fork cluster crashes R
if (multithread == "Julia") {

  cl <- NULL

} else if (multithread == "R") {

  cluster <- "socket" # "fork"
  if (cluster == "fork") {

    cl <- rsockets
    stop("The fork cluster crashes R.")

  } else if (cluster == "socket") {

    # Spawn a separate Julia process on each core (~3 mins)
    tic()
    cl <- makeCluster(rsockets)
    ignore <- c("spatw", "sims_for_performance", "sims_for_performance_ls")
    export <- ls()[!(ls() %in% c())]
    clusterExport(cl, export)
    clusterEvalQ(cl, {

      # Load packages on each core
      library(collapse)
      library(dv)
      library(patter)
      library(data.table)
      library(dtplyr)
      library(dplyr, warn.conflicts = FALSE)
      library(tictoc)

      # Each core uses single threaded DT and Julia implementations
      setDTthreads(threads = 1)
      julia_connect(.threads = 1)

      # Set Julia objects on each core
      spatw <- readRDS(here_input("spatw.rds"))
      set_map(terra::unwrap(spatw))
      set_seed()
    })
    toc()

  }

}

#### Check cluster settings
# We are multi-threading:
nrow(sims)
multithread
# We are using n Julia threads
check_multithreading(multithread)
# We are using the following sigma bandwidth estimator:
formals(get_ud_patter)$sigma

#### Estimate UDs
gc()
tic()
sdt <-
  pbapply::pblapply(seq_len(nchunks), cl = cl, function(i) {

    #### Set up chunk
    # Logs
    # i <- 1
    log.txt <- here_output("log", "patter", type, paste0("log-", i, ".txt"))
    cat_log <- patter:::cat_init(.verbose = log.txt)
    # rstudioapi::navigateToFile(log.txt)
    cat_log(paste("CHUNK" , i, "\n"))
    t1_chunk <- Sys.time()
    cat_log(paste("Start:", as.character(t1_chunk), "\n"))
    # Grid
    spat <- terra::unwrap(spatw)
    # Simulations for chunk
    sims_for_chunk <- sims[chunks[[i]], ]

    #### Run simulations for chunk
    success <-
      lapply(split(sims_for_chunk, seq_len(nrow(sims_for_chunk))), function(sim) {
        # sim = sims_for_chunk[1, ]; spat = terra::unwrap(spatw)
        cat_log(paste0("\n", sim$row, ":\n"))
        t1 <- Sys.time()
        s <- workflow_patter(sim = sim,
                             spat = spat,
                             win = win)
        t2 <- Sys.time()
        difftime(t2, t1, units = "secs")
        cat_log(as.numeric(round(difftime(t2, t1, units = "secs"), 2)))
        s
      }) |> rbindlist()
    t2_chunk <- Sys.time()

    #### Record timings & return success
    cat_log("\n\n")
    cat_log(paste("End:", as.character(t2_chunk), "\n"))
    cat_log(paste("Duration:", round(difftime(t2_chunk, t1_chunk, units = "hours"), 2), "hours."))
    sink()
    success

  })
toc()

#### Record success
sdt <- rbindlist(sdt)
table(sdt$acpf)
table(sdt$acdcpf)
saveRDS(sdt, here_data("sims", "output", "success", paste0("patter-", type, ".rds")))


#########################
#########################
#### Quick visuals

if (interactive()) {

  qplot <- function(file) {
    terra_qplot(file)
    p <- read_path(sim)
    prettyGraphics::add_sp_path(p$x, p$y, length = 0, lwd = 0.1)
    m <- read_array(sim)
    points(m$receiver_x, m$receiver_y)
  }

  # Plot UDs
  pp <- par(mfrow = c(3, 2))
  qplot(ud_path)
  qplot(here_alg(sim, "coa", "120 mins", "ud.tif"))
  qplot(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-f.tif"))
  qplot(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-s.tif"))
  qplot(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-f.tif"))
  qplot(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-s.tif"))
  par(pp)

  # Plot skill (MB)
  ud_path <- here_alg(sim, "path", "ud.tif") |> terra::rast()
  ud_alg  <- c(here_alg(sim, "coa", "120 mins", "ud.tif") |> terra::rast(),
               here_alg(sim, "patter", "acpf", sim$alg_par, "ud-f.tif") |> terra::rast(),
               here_alg(sim, "patter", "acpf", sim$alg_par, "ud-s.tif") |> terra::rast(),
               here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-f.tif") |> terra::rast(),
               here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-s.tif") |> terra::rast())
  barplot(skill_by_alg(ud_alg, ud_path, .f = skill_me),
          names.arg = c("COA", "ACPF (F)", "ACPF (S)", "ACDCPF (F)", "ACDCPF (S)"))

}


#### End of code.
#########################
#########################
