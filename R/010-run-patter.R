#########################
#########################
#### run-patter.R

#### Aims
# 1) Estimate UDs using patter (for performance & sensitivity simulations)

#### Prerequisites
# 1) Simulate data & prepare for algorithm implementation

#### Log
# Performance simulations (n = 1181):
# * Machine         : SIA-LAVENDED
# * Total wall time : 26.5 hours (7.15 hrs for algorithms only) on 12 CPU
# * Convergence     : ACDPF (1 failure); ACDCPF (1 failure)
# * File transfer   : NA
# Sensitivity simulations
# * Batch size: 5000 simulations (5000 * 3 * 1.8 / 1e3 MB = 27 GB)
# * 1:1000        : PF-1; 7 hrs; copied;
# * 1001:2000     : PF-2; 32 hours; copied;
# * 2001:5000     : PF-1; in progress (working);
# * 5001:10000    : PF-3: in progress (RETRYING, TO CHECK);
# * 10001:15000   : PF-2: TO DO (RETRYING, TO CHECK);


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
library(JuliaCall)
library(microbenchmark)
library(parallel)
library(testthat)
library(tictoc)
dv::src()

#### Load data
Sys.time()
spatw      <- readRDS(here_input("spatw.rds"))
win        <- qs::qread(here_input("win.qs"))
sims_for_performance <- readRDS(here_input("sims-performance.rds"))
sims_for_sensitivity <- readRDS(here_input("sims-sensitivity.rds"))


#########################
#########################
#### Julia connection

#### Define sims
#
# Define type
# type      <- "performance"
type        <- "sensitivity"
performance <- ifelse(type == "performance", TRUE, FALSE)
#
# Performance simulations:
# sims     <- copy(sims_for_performance)
# batch    <- 1:nrow(sims)
#
# Sensitivity simulations:"
sims     <- copy(sims_for_sensitivity)
batch    <- 1:10L
#
# Subset sims by batch & define batch_id
sims     <- sims[batch, ]
batch_id <- min(batch)
#
# Cleanup:
rm(sims_for_performance, sims_for_sensitivity)
#
# (optional) Update number of particles
sims[, n_particles := NULL]
# sims[, n_particles := 5e4L]

# Visualise n particles
terra::plot(terra::unwrap(spatw))
points(terra::spatSample(terra::unwrap(spatw),
                         size = 30000L,
                         xy = TRUE, values = FALSE),
       pch = ".")

#### Define multi-threading option
# We can multi-thread R or Julia code
# Multi-threading in R uses a socket cluster (forking crashes the R session)
multithread <- c("R", "Julia")
multithread <- multithread[1]

#### Connect to Julia
if (multithread == "R") {

  rsockets <- 12L
  # rsockets <- 16L

  setDTthreads(threads = 1)
  # julia_connect(.threads = 1)

} else {

  rsockets <- 1L

  setDTthreads(threads = Sys.getenv("OMP_NUM_THREADS"))
  julia_connect()

  set_seed()
  set_map(terra::unwrap(spatw))

  JuliaCall::julia_command(ModelObsAcousticContainer)
  JuliaCall::julia_command(ModelObsAcousticContainer.logpdf_obs)

}


#########################
#########################
#### Test routines

# This code validates the workflow_patter() function used to run simulations (~20 mins).
# > All tests passed.

if (FALSE & multithread == "Julia") {

  #### Run tests for a random selection of simulations
  tic()
  for (i in c(1, 2, 140, 204, 595, 1001)) {

    #### Implement workflow_patter() with `test = TRUE`
    # i <- 1
    message(paste("Test", i))
    sim  <- sims[i, ]
    test <- workflow_patter(sim = sim,
                            spat = terra::unwrap(spatw),
                            win = win,
                            performance = performance,
                            test = TRUE)

    #### Verify alignment of input datasets
    # All time series should be defined between the first/last detection
    test$input$path
    test$input$acoustics
    test$input$archival
    start <- c(min(test$input$timeline),
               min(test$input$path$timestamp),
               min(test$input$acoustics$timestamp),
               min(test$input$archival$timestamp))
    end <- c(max(test$input$timeline),
             max(test$input$path$timestamp),
             max(test$input$acoustics$timestamp),
             max(test$input$archival$timestamp))
    expect_true(length(unique(start)) == 1L)
    expect_true(length(unique(end)) == 1L)
    expect_equal(test$input$timeline, test$input$path$timestamp)
    expect_equal(test$input$timeline, test$input$archival$timestamp)

    #### Verify incorporation of detections (~10 s): ok.
    cl_lapply(c("acpf", "acdcpf"), function(algorithm) {

      lapply(c("pff", "pfb", "smo"), function(fun) {

        # algorithm <- "acdcpf"; fun <- "pff"
        print(paste(algorithm, fun))

        # Define particles
        particles <- test$output[[algorithm]]$particles[[fun]]$states

        # Validate map_value
        particles[, map_value_terra := terra::extract(terra::unwrap(spatw), cbind(particles$x, particles$y))[, 1]]
        head(particles[, .(map_value, map_value_terra)])
        expect_true(max(abs(particles$map_value - particles$map_value_terra)) < 1e5)

        # Define detections
        detections <-
          sim |>
          read_acoustics() |>
          filter(obs == 1L) |>
          select(timestamp, receiver_x, receiver_y) |>
          as.data.table()

        # Define corresponding states
        states <-
          particles |>
          filter(timestamp %in% detections$timestamp) |>
          select(timestamp, x, y) |>
          as.data.table()

        # Verify that all states are within detection ranges (sim$gamma)
        # * Note there may be multiple detections at any one time step that we pair against multiple states
        any(duplicated(detections$timestamp))
        dist_from_receiver <-
          full_join(detections, states, relationship = "many-to-many") |>
          mutate(dist = terra::distance(cbind(receiver_x, receiver_y), cbind(x, y), lonlat = FALSE, pairwise = TRUE)) |>
          pull(dist)
        expect_true(all(dist_from_receiver <= sim$gamma))

      })

    })

    #### Verify archival data
    cl_lapply(c("pff", "pfb", "smo"), function(fun) {

      # Define particles
      particles <- test$output$acdcpf$particles[[fun]]$states

      # Define path
      path <- read_path(sim)

      # Verify that all state map_value is within ± 5 m of the observed depth
      v <- left_join(particles, path, by = "timestamp", suffix = c("_state", "_path")) |> as.data.table()
      v[, valid := depth >= map_value_state - 5 & depth <= map_value_state + 5]
      expect_true(all(v$valid))

    })

  }
  toc()

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

# Particle file size:
# * We require 23 MB to store particles from one run of the forward filter:
# * 23.7 MB
# * Therefore, we require > 2000 GB (2 TB) to store particles from all filter runs!
# * 2844 GB (23.7 * 2 *  2 * 30000 / 1e3)
# * 23.7 * 2 (forward, smoothed particles) * 2 (ACPF, ACDCPF) * 30000 simulations
# > We do not currently save particles

# UD file size
# * Each UD is ~1.8 MB and each simulation saves 4 tifs
# * We require ~ 200 GB for all outputs

if (FALSE && multithread == "Julia") {

  #### Define simulation
  # (TO DO) Repeat the code with a few random choices to check result consistency
  # sim <- sims[row == 205, ]
  # sim <- sims[row == 325, ]
  # sim <- sims[row == 47, ]  # 10,000 ok
  # sim <- sims[row == 145, ] # 10,000 ok
  # sim <- sims[row == 349, ] # 10,000 sometimes ok

  #### Read data
  acoustics  <- read_acoustics(sim); # tmp <- copy(acoustics)
  archival   <- read_archival(sim)

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

  #### Define observations
  # yobs      <- list(acoustics)
  # model_obs <- "ModelObsAcousticLogisTrunc"
  containers  <- assemble_acoustics_containers(acoustics,
                                               .direction = "forward",
                                               .mobility = sim$mobility)
  # yobs        <- list(acoustics, containers)
  # model_obs   <- c("ModelObsAcousticLogisTrunc", "ModelObsAcousticContainer")
  yobs        <- list(acoustics, containers, archival)
  model_obs   <- c("ModelObsAcousticLogisTrunc", "ModelObsAcousticContainer", "ModelObsDepthUniform")

  #### Filter args
  args <- list(.map = spat,
               .timeline = timeline,
               .state = "StateXY",
               .xinit = NULL,
               .xinit_pars = list(mobility = sim$mobility),
               .yobs = yobs,
               .model_obs = model_obs,
               .model_move = model_move,
               .n_particle = 10000L, #sim$n_particles,
               .n_resample = 10000, #sim$n_particles,
               .verbose = FALSE)

  #### Filter speeds (1 thread):
  # The number of particles is a balance between speed & convergence.
  # Other options are to tweak the resampling or use acoustic containers.
  # * 10,000 particles: 3 secs
  # * 25,000 particles: 9 secs
  # * 50,000 particles: 18 secs

  #### Forward filter
  # Run filter
  out_pff  <- do.call(pf_filter, args, quote = TRUE)
  # Check file size
  tmp      <- tempfile(fileext = ".qs")
  qs::qsave(out_pff$states, tmp)
  file.size(tmp) / 1e6
  unlink(tmp)

  #### Backward filter
  containers      <- assemble_acoustics_containers(acoustics,
                                                   .direction = "backward",
                                                   .mobility = sim$mobility)
  # args$.yobs      <- list(acoustics, containers)
  args$.yobs      <- list(acoustics, containers, archival)
  args$.direction <- "backward"
  out_pfb         <- do.call(pf_filter, args, quote = TRUE)

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
  sims <- sims[1:400L, ]

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
  # NB: These timings are approximate (routines updated subsequently)

  # 1 algorithm run, 1 thread (bw.diggle):
  # * Forward filter:    4 s
  # * Backward filter:   4 s
  # * Smoother:         90 s
  # * UD estimation:    60 s (2 * 30 for the forward & smoothed UDs)
  #                    158 s (* 2 = 316 s for ACPF & ACDCPF)

  # 1 algorithm run, 1 thread (bw.scott):
  # (not currently implemented)

  # Multi-threaded R (bw.diggle):
  # 100 simulations   : 1.34 hours (4854.023 s)
  # 1181 simulations  : 15.9 hours
  # Total ETA         : 16 days

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
  file.create(here_output("log", "patter", type, paste0(batch_id, "-log-", i, ".txt")))
}) |> invisible()

#### Initialise cluster
# We use a socket cluster, as the fork cluster crashes R
if (multithread == "Julia") {

  cl <- NULL
  check_multithreading(multithread)

} else if (multithread == "R") {

  cluster <- "socket" # "fork"
  if (cluster == "fork") {

    cl <- rsockets
    stop("The fork cluster crashes R.")

  # Spawn a separate Julia process on each core (~3 mins)
  } else if (cluster == "socket") {

    tic()

    # Initialise Julia outside of the cluster (see src/julia.R)
    julia_connect(.threads = 1L)

    # Make cluster
    cl <- makeCluster(rsockets)

    # clusterExport (data & parameters)
    export <- c(
        # items for clusterEvalQ
        "multithread", "ModelObsAcousticContainer", "ModelObsAcousticContainer.logpdf_obs",
        # Items for simulation
        "type", "batch_id", "sims", "chunks", "win", "performance")
    # export_ls <- as.list(export)
    # names(export_ls) <- export
    # sapply(export_ls, lobstr::obj_size) |> sort()
    clusterExport(cl, export)

    # clusterEvalQ (libraries, packages & Julia set up)
    clusterEvalQ(cl, {

      # Load packages on each core
      library(collapse)
      library(dv)
      library(patter)
      library(data.table)
      library(dtplyr)
      library(dplyr, warn.conflicts = FALSE)
      library(tictoc)
      dv::src()

      # Each core uses single threaded DT and Julia implementations
      setDTthreads(threads = 1)
      try_julia_connect()
      check_multithreading(multithread)

      # Set Julia objects on each core
      JuliaCall:::.julia$cmd("using RCall")
      JuliaCall::julia_command(ModelObsAcousticContainer)
      JuliaCall::julia_command(ModelObsAcousticContainer.logpdf_obs)
      set_map(terra::rast(here_data("sims", "input", "spat.tif")))
      set_seed()

      # Check memory used (101082480 bytes: ~100 MB)
      # stop(lobstr::mem_used())

      NULL
    })

    toc()

  }

}

#### Check cluster settings
# We are running:
glue::glue("{type} simulations (performance: {performance}) for batch {batch_id} ({min(batch)}:{max(batch)})")
# We are multi-threading:
multithread
# We are using nchunks:
nchunks
# For n simulations
nrow(sims)
# We are using the following sigma bandwidth estimator:
formals(get_ud_patter)$sigma
# The first files to be created should appear here:
here_alg(sims[1, ], "patter", "acpf", sims$alg_par[1])
# > Check all files are created after a few mins:
# > time.qs, ud-f.tif, ud-s.tif

#### Estimate UDs
gc()
Sys.time()
tic()
sdt <-
  pbapply::pblapply(seq_len(nchunks), cl = cl, function(i) {

    #### Set up chunk
    # Logs
    # i <- 1
    log.txt <- here_output("log", "patter", type, paste0(batch_id, "-log-", i, ".txt"))
    cat_log <- patter:::cat_init(.verbose = log.txt)
    # rstudioapi::navigateToFile(log.txt)
    cat_log(paste("CHUNK" , i, "\n"))
    t1_chunk <- Sys.time()
    cat_log(paste("Start:", as.character(t1_chunk), "\n"))
    # Grid
    spat <- terra::rast(here_input("spat.tif"))
    # Simulations for chunk
    sims_for_chunk <- sims[chunks[[i]], ]

    #### Run simulations for chunk
    success <-
      lapply(split(sims_for_chunk, seq_len(nrow(sims_for_chunk))), function(sim) {
        # sim = sims_for_chunk[1, ]; spat = terra::unwrap(spatw)
        cat_log(paste0("\n", sim$row, ":\n"))
        t1 <- Sys.time()
        s <- workflow_patter(sim = sim,
                             spat = spat, win = win,
                             performance = performance)
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
Sys.time()

#### Record success (200 trials):
#
# NB: These timings are approximate (routines updated subsequently)
#
# Without acoustic containers:
# * With 10,000 particles (~0.3 hrs):
# - ACPF  : 6 failures
# - ACDCPF: 28 failures
# * With 50,000 particles (~1.5 hrs):
# - ACPF  : success
# - ACDCPF: 8 failures
# * With 50,000 particles and infrequent resampling (.n_resample = 1000):
# - ACPF  : success
# - ACDCPF: 18 failures
#
# With acoustic containers:
# * With 10,000 particles (~0.5 hrs):
# - ACPF  : success
# - ACDCPF: 9 failures
# * With 50,000 particles (~1.75 hrs, computer busy):
# - ACPF  : success
# - ACDCPF: success
#
# Use different numbers of particles for the algorithms (with containers):
# * With 5,000 (ACPF) and 20,000 (ACDCPF) particles, 400 sims (~1 hr):
# - ACPF  : success
# - ACDCPF: 3 failures (row: 47, 145, 349)
# * With 10,000 (ACPF) and 30,000 (ACDCPF) particles, 400 sims (~1.75 hrs):
# - ACPF  : success
# - ACDCPF: 2 failures
#
# Conclusions
# * Increasing particles, acoustic containers & regular re sampling all help convergence
# * With small numbers of particles, acoustic containers help (both ACPF & ACDCPF)
# * With larger numbers of particles, ACPF converges but ACDCPF convergence benefits from containers
# * Re-sampling every time step appears to be beneficial for convergence
#
# Strategy
# * ACPF:
# > 5,000 particles was sufficient for convergence (with containers)
# > Use 5,000 particles + containers + resampling every time step
# * ACDCPF:
# > 10,000 particles insufficient
# > 20,000 particles almost always sufficient
# > Use 30,000 particles

#### Record success (simulations)
#
# Performance simulations (n = 1181)
# * Convergence: ACPF (1 failure: row 872); ACDCPF (1 failure: row 472);
# * Time: 1.11 days (12 socket clusters; MacBook)
#

sdt <- rbindlist(sdt)
table(sdt$acpf)
table(sdt$acdcpf)
saveRDS(sdt, here_data("sims", "output", "success", paste0("patter-", type, "-", batch_id, ".rds")))

#### Quick examination of convergence failures
# See also patter-debug.R
# sdt <- readRDS(here_data("sims", "output", "success", paste0("patter-", type, "-", batch_id, ".rds")))
if (FALSE) {

  sdt[acpf == FALSE, ]
  sdt[acdcpf == FALSE, ]
  sim <- sims[sims$row == 66, ]
  vis_path  <- get_path(sim, readRDS(here_input("paths.rds")), qs::qread(here_input("detections.qs")))
  vis_array <- get_array(sim, readRDS(here_input("arrays.rds")))
  terra::plot(spat)
  spatPoints(vis_path$x, vis_path$y)
  spatPoints(vis_array$receiver_x, vis_array$receiver_y, cex = 2, col = "red")

}


#########################
#########################
#### Quick visuals

if (FALSE && interactive()) {

  qplot <- function(file) {
    terra_qplot(file)
    p <- read_path(sim)
    prettyGraphics::add_sp_path(p$x, p$y, length = 0, lwd = 0.1)
    m <- read_array(sim)
    points(m$receiver_x, m$receiver_y)
  }

  # Define UD paths
  sim     <- sims[1000, ]
  ud_path <- here_alg(sim, "path", "ud.tif")
  ud_alg  <- c(here_alg(sim, "coa", "30 mins", "ud.tif"),
               here_alg(sim, "coa", "120 mins", "ud.tif"),
               here_alg(sim, "rsp", "default", "ud.tif"),
               here_alg(sim, "rsp", "custom", "ud.tif"),
               here_alg(sim, "patter", "acpf", sim$alg_par, "ud-f.tif"),
               here_alg(sim, "patter", "acpf", sim$alg_par, "ud-s.tif"),
               here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-f.tif") ,
               here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-s.tif"))

  # Plot UDs
  pp <- par(mfrow = c(5, 2))
  qplot(ud_path)
  lapply(ud_alg, qplot)
  par(pp)

  # Estimate file size
  (sapply(ud_alg, file.size) / 1e6 ) |> unname()

  # Plot skill (MB)
  ud_path <- ud_path |> terra::rast()
  ud_alg  <- sapply(ud_alg, function(f) terra::rast(f))
  barplot(skill_by_alg(ud_alg, ud_path, .f = skill_me),
          names.arg = c("COA (1)", "COA (2)",
                        "RSP (1)", "RSP (2)",
                        "ACPF (F)", "ACPF (S)",
                        "ACDCPF (F)", "ACDCPF (S)"))

}


#### End of code.
#########################
#########################
