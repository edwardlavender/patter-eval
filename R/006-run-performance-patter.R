#########################
#########################
#### run-performance-patter.R

#### Aims
# 1) Estimate UDs using patter for performance simulations (correct parameters)

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
library(tictoc)
dv::src()

#### Load data
spatw      <- readRDS(here_input("spatw.rds"))
im         <- qs::qread(here_input("im.qs"))
win        <- qs::qread(here_input("win.qs"))
sims_for_performance    <- readRDS(here_input("sims-performance.rds"))
sims_for_performance_ls <- split(sims_for_performance, sims_for_performance$id)


#########################
#########################
#### Estimate UDs

#### Set up cluster
# Number of forks
cl <- 1L
# Define chunks to iterate over in parallel
chunks  <- patter:::cl_chunks(cl, nrow(sims_for_performance))
nchunks <- length(chunks)
# Create a log file for each chunk
lapply(seq_len(nchunks), function(i) {
  file.create(here_output("logs", paste0("log-", i, ".txt")))
}) |> invisible()

#### Time trials
# Estimated duration on {cl} CPUs, assuming simulations take {guess} seconds
guess <- 30 # 30 s
(nrow(sims_for_performance) * 30)/60/60/cl    # hours
(nrow(sims_for_performance) * 30)/60/60/24/cl # days

#### Implementation
gc()
success <-
  pbapply::pblapply(seq_len(nchunks), cl = NULL, function(i) {

    #### Set up chunk
    # Logs
    log.txt <- here_output("logs", paste0("log-", i, ".txt"))
    sink(log.txt)
    # rstudioapi::navigateToFile(log.txt)
    cat(paste("CHUNK" , i, "\n"))
    t1_chunk <- Sys.time()
    cat(paste("Start:", as.character(t1_chunk), "\n"))
    # Grid
    spat <- terra::unwrap(spatw)
    # Simulations for chunk
    sims_for_chunk <- sims[chunks[[i]], ]

    #### Run simulations for chunk
    success <-
      lapply(split(sims_for_chunk, seq_len(nrow(sims_for_chunk))), function(sim) {
        cat(paste0("\n", sim$row, ":\n"))
        t1 <- Sys.time()
        workflow_patter(sim = sim,
                        spat = spat,
                        im = im,
                        win = win)
        t2 <- Sys.time()
        cat(as.numeric(difftime(t2, t1, units = "secs")))
      }) |> unlist()
    t2_chunk <- Sys.time()

    #### Record timings & return success
    cat("\n\n")
    cat(paste("End:", as.character(t2_chunk), "\n"))
    cat(paste("Duration:", round(difftime(t2_chunk, t1_chunk, units = "hours"), 2), "hours."))
    sink()
    success

  }) |> unlist()

#### Validate convergence
# TO DO
table(success)


#########################
#########################
#### Quick visuals

qplot <- function(file) {
  terra_qplot(file)
  p <- read_path(sim)
  prettyGraphics::add_sp_path(p$x, p$y, length = 0, lwd = 0.1)
  m <- read_array(sim)
  points(m$receiver_easting, m$receiver_northing)
}

sim <- sims[4, ]
pp <- par(mfrow = c(2, 2))
qplot(here_alg(sim, "path", "ud.tif"))
qplot(here_alg(sim, "coa", "120 mins", "ud.tif"))
qplot(here_alg(sim, "patter", "acpf", sim$alg_par, "ud.tif"))
qplot(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud.tif"))
par(pp)


#### End of code.
#########################
#########################
