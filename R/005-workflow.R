#########################
#########################
#### workflow.R

#### Aims:
# (1) Define a function to implement algorithms using simulated data
# (2) Record model skill

#### Prerequisites
# 1) This function is called via 006-run-sims.R

#### TO DO
# * Add RSP
# * Update pf_dens() function for improved speed
# * Update AC/ACPF routines
# * Add ACDCPF

workflow <- function(
    # input datasets
    sim,
    paths,
    detections,
    seed,
    # algorithm parameters
    sigma,
    # utility functions, passed as arguments
    unwrapr,
    mins,
    here_input,
    here_output,
    here_alg,
    pf_coords,
    skill_by_alg
    ) {


  #########################
  #########################
  #### Prepare inputs

  #### Record script duration
  t_onset <- Sys.time()

  #### Define log file
  log.txt <- here_alg(sim$id, "log.txt")
  file.create(log.txt)
  # rstudioapi::navigateToFile(log.txt)
  sink(log.txt)
  on.exit(sink(), add = TRUE)
  cat(paste0(Sys.time(), "\n"))
  print(sim)

  #### Define path (& realisation, if necessary)
  cat("* Data preparation...\n")
  path <-
    paths[[sim$path_type]] |>
    filter(path_id == sim$path_realisation) |>
    as.data.table()

  #### Define path-level algorithm inputs
  # Define observations
  # * detections is indexed by:
  #   - detection pars & path type pair
  #     - array design
  acoustics <-
    detections[[sim$system_type]][[sim$array_type]] |>
    filter(array_id == sim$array_realisation) |>
    filter(path_id == sim$path_type) |>
    as.data.table()
  if (nrow(acoustics) <= 5L) {
    cat("> Insufficient detections to proceed: ending simulation.\n")
    return(NULL)
  }
  # Revise path (focus on portion between first & last detections)
  path <-
    path |>
    filter(timestamp >= min(acoustics$timestamp) &
             timestamp <= max(acoustics$timestamp)) |>
    as.data.table()
  if (nrow(path) <= 5L) {
    cat("> Insufficient path points to proceed: ending simulation.\n")
    return(NULL)
  }

  #### Define array-level algorithm inputs
  overlaps   <-
    here_input("ac", sim$gamma, "overlaps.rds") |>
    readRDS()
  kernels    <-
    here_input("ac", sim$gamma, sim$alpha, sim$beta, "kernels.rds") |>
    readRDS() |>
    unwrapr()


  #########################
  #########################
  #### Truth

  #### Estimate UD for simulated path
  # We implement this once for each array/simulated path combination
  if (sim$performance) {
    cat("* Simulated patterns...\n")
    ud_path <-
      pf_dens(.xpf = grid,
              .coord = path[, .(x, y)],
              .plot = FALSE,
              .verbose = FALSE,
              sigma = sigma)
    if (is.null(ud_path)) {
      stop(paste0("`pf_dens()` failed for simulation ID ", sim$id, "."))
    }
    terra::writeRaster(ud_path,
                       here_alg(sim$id, "path", "ud.tif"),
                       overwrite = TRUE)
  }


  #########################
  #########################
  #### COA algorithm

  if (sim$performance) {

    #### Define .delta_t parameter
    cat("* COA algorithm...\n")
    # (optional) TO DO

    #### Implement COA algorithm for selected .delta_t parameters
    ud_coas <- lapply(c("30 mins", "120 mins"), function(.dt) {
      if (difftime(max(acoustics$timestamp), min(acoustics$timestamp), units = "mins") <= .dt) {
        return(NULL)
      }
      out_coa <- coa(acoustics,
                     .delta_t = .dt,
                     .plot_weights = FALSE)
      ud_coa  <- pf_dens(grid,
                         .coord = out_coa[, .(x = coa_x, y = coa_y)],
                         .plot = FALSE,
                         .verbose = FALSE,
                         sigma = sigma)
      if (!is.null(ud_coa)) {
        terra::writeRaster(ud_coa,
                           here_alg(sim$id, "coa", .dt, "ud.tif"),
                           overwrite = TRUE)
      }
      ud_coa
    })
  }


  #########################
  #########################
  #### RSPs

  if (sim$performance) {

    #### TO DO

  }


  #########################
  #########################
  #### flapper algorithms

  #### Define inputs
  # Define observations, account for grid resolution
  # This is fast enough that we don't pre-compute it
  cat("* flapper algorithms...\n")
  obs <- acs_setup_obs(acoustics,
                       .step = paste(sim$step, "mins"),
                       .mobility = sim$mobility)


  #########################
  #### ACPF

  #### Implement AC algorithm
  # User verbose = FALSE for speed
  # Use tryCatch because this may fail for mis-specified implementations
  cat("\t * AC algorithm...\n")
  out_ac <- tryCatch(
    acs(obs,
        .bathy = grid,
        .detection_overlaps = overlaps,
        .detection_kernels = kernels,
        .save_record = TRUE,
        .verbose = FALSE),
    error = function(e) e)
  if (inherits(out_ac, "error")) {
    cat(paste("\t AC algorithm failed with error:", out_ac))
    return(NULL)
  }

  #### Implement PF
  # Run forward simulation
  print("\t > Forward simulation...")
  t1_ac_pff <- Sys.time()
  set.seed(seed)
  out_pff <- pf_forward(obs,
                        .record = out_ac$record,
                        .kick = pf_kick,
                        .shape = sim$shape, .scale = sim$scale, .mobility = sim$mobility,
                        .rho = sim$rho, .sd = sim$sd,
                        .bathy = grid,
                        .n = sim$n_particles,
                        .save_history = TRUE,
                        .verbose = FALSE)
  if (!inherits(out_pff, "pf")) {
    print("\t \t > Convergence failure: ending simulation.")
    return(NULL)
  }
  t2_ac_pff <- Sys.time()

  # Run backward pass
  print("\t > Backward simulation...")
  t1_ac_pfb <- Sys.time()
  out_pfb <- pf_backward(out_pff$history,
                         .save_history = TRUE,
                         .verbose = FALSE)
  t2_ac_pfb <- Sys.time()
  # Implement smoothing
  t1_ac_ud <- Sys.time()
  print("\t > Smoothing...")
  out_pfp <- pf_coords(out_pfb$history)
  ud_acpf <- pf_dens(grid,
                     .coord = out_pfp,
                     .plot = FALSE,
                     .verbose = FALSE,
                     sigma = sigma)
  if (!is.null(ud_acpf)) {
    terra::writeRaster(ud_acpf,
                       here_alg(sim$id, "patter", "acpf", "ud.tif"),
                       overwrite = TRUE)
  }
  t2_ac_ud <- Sys.time()


  #########################
  #### ACDCPF

  #### Implement ACDCPF
  # TO DO


  #########################
  #########################
  #### Calculate skill

  #### Define template data.table
  cat("* Model skill...")
  skill <- data.table(
    id = sim$id,
    performance = sim$performance,
    alg = c("COA(30)", "COA(120)", "ACPF", "ACDCPF"),
    mb = NA_real_,
    me = NA_real_,
    rmse = NA_real_,
    R = NA_real_,
    d = NA_real_
  )

  #### Calculate skill metrics
  ud_acdcpf <- NULL
  if (sim$performance) {
    # We will calculate model skill scores for each algorithm
    uds <- list(ud_coas[[1]], ud_coas[[2]], ud_acpf, ud_acdcpf)
  } else {
    # We will calculate model skill scores for patter algorithms only
    uds <- list(NULL, NULL, ud_acpf, ud_acdcpf)
  }
  skill$mb   <- skill_by_alg(uds, ud_path, skill_mb)
  skill$me   <- skill_by_alg(uds, ud_path, skill_me)
  skill$rmse <- skill_by_alg(uds, ud_path, skill_rmse)
  skill$R    <- skill_by_alg(uds, ud_path, skill_R)
  skill$d    <- skill_by_alg(uds, ud_path, skill_d)
  # Save metrics
  saveRDS(skill, here_output("skill", paste0(sim$id, ".rds")))

  #### Record wall time
  cat("* Wall time...")
  t_end <- Sys.time()
  time <- data.table(id = sim$id,
                     ac_pff = mins(t2_ac_pff, t1_ac_pff),
                     ac_pfb = mins(t2_ac_pfb, t1_ac_pfb),
                     ac_ud = mins(t2_ac_ud, t1_ac_ud),
                     script_dur = mins(t_end, t_onset, units = "mins"))
  saveRDS(time, here_output("time", paste0(sim$id, ".rds")))


  #### Return NULL
  NULL
}


#### End of code.
#########################
#########################
