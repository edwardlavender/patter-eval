#########################
#########################
#### run-server.R

#### Aims
# 1) Server management
# * We implement the algorithms on the server
# * This script copies outputs onto the machine
# * Then cleans up the server for a subsequent run of the algorithms

#### Prerequisites
# 1) Mount patter-eval on the server (@ \\tsclient\patter-eval\)


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
library(data.table)
library(patter)
library(stringr)
library(testthat)
library(tictoc)
dv::src()

#### Load data
# sims <- readRDS(here_input("sims-performance.rds"))
sims <- readRDS(here_input("sims-sensitivity.rds"))


#########################
#########################
#### Management

#### (conditional) Clean up old files on server
# TO DO

#### Define sims
type     <- "sensitivity"
batch    <- 1:1000L
batch_id <- min(batch)

#### List patter outputs (~5 s)
tic()
# List logs
# logs    <- list.files(here_output("log", "patter", "performance"))
logs <- list.files(here_output("log", "patter", "sensitivity"))
head(logs)
# List outputs
outputs <- list.files(here_output("run"), recursive = TRUE)
outputs <- outputs[str_detect(outputs, "patter")]
head(outputs)
toc()

#### Compare computed outputs to expected outputs for a selected sims batch
# Run this code:
# * Before file transfer (on server);
# * After file transfer (on local machine);
# * Then it is safe to delete files on the server
# This is the set of directories that we expect to contain outputs (see here_alg()):
sims  <- sims[batch, ]
# algorithm <- "acpf"
algorithm   <- "acdcpf"
folders <- file.path(sims$combination, sims$array_type, sims$array_realisation,
                     sims$path_realisation, "patter", "acpf", sims$alg_par)
# These are the ones that have files attached (~3 s):
expected <- data.table(folder = folders, success = FALSE)
tic()
for (i in seq_row(expected)) {
  # print(i)
  # any(str_detect(outputs, expected$folder[1]))
  # Check for any outputs in the expected folder
  expected[i, success := any(str_detect(outputs, expected$folder[i]))]
}
toc()
table(expected$success)
# Check specifically for convergence.rds
expected <- file.path(sims$combination, sims$array_type, sims$array_realisation,
                      sims$path_realisation, "patter", "acpf", sims$alg_par,
                      "convergence.rds")
table(expected %in% outputs)

#### (lavended) Post-transfer checks
if (lavended()) {
  # Verify that expected outputs exist
  expect_true(all(expected %in% outputs))
  # Verify that we can read example SpatRasters
  tifs <- outputs[str_detect(outputs, ".tif")]
  int  <- sample.int(length(tifs), size = 10)
  tifs <- lapply(tifs[int], function(f) terra::rast(here_output("run", f)))
  terra::plot(tifs[[1]])
}

#### Transfer log.txt files
to_dir <- "\\\\tsclient\\patter-eval\\data\\sims\\output\\log\\patter\\sensitivity"
expect_true(dir.exists(to_dir))
# Transfer files
tic()
success <-
  sapply(logs, function(log.txt) {

    # Define file to copy
    # log.txt <- logs[1]
    print(log.txt)
    from <- file.path(here_output("log", "patter", "sensitivity"), log.txt)
    expect_true(file.exists(from))

    # Define file to create
    to     <- file.path(to_dir, log.txt)

    # Copy file
    file.copy(from, to, overwrite = TRUE)
  })
toc()
# Validate success
table(success)
expect_true(all(success))

#### Transfer success record (sdt) for batch
to_dir <- "\\\\tsclient\\patter-eval\\data\\sims\\output\\success"
expect_true(dir.exists(to_dir))
sdtname <- paste0("patter-", type, "-", batch_id, ".rds")
success <- file.copy(here_data("sims", "output", "success", sdtname),
                     file.path(to_dir, sdtname))
expect_true(success)

#### Transfer outputs
# Define to directory
to_dir <- "\\\\tsclient\\patter-eval\\data\\sims\\output\\run"
expect_true(dir.exists(to_dir))
# Set up cluster
cl <- NULL
# cl <- parallel::makeCluster(2L)
# parallel::clusterExport(cl, "here_output")
# Copy files
tic()
success <-
  cl_lapply(seq_len(length(outputs)),
            .cl = cl,
            .fun = function(i) {
              print(i)
              output <- outputs[i]
              from   <- file.path(here_output("run"), output)
              expect_true(file.exists(from))
              to <- file.path(to_dir, output)
              file.copy(from, to, overwrite = TRUE)
            })
toc()
# Validate success
table(success)
expect_true(all(success))

#### Clean up patter outputs (~5 s)
# This code is guarded by if(FALSE)
if (FALSE) {

  tic()
  # Logs
  # logs_path <- file.path(here_output("log", "patter", "performance"), logs)
  logs_path <- file.path(here_output("log", "patter", "sensitivity"), logs)
  expect_true(file.exists(logs_path[1]))
  unlink(logs_path)
  expect_false(file.exists(logs_path[1]))
  # Outputs
  outputs_path <- file.path(here_output("run"), outputs)
  expect_true(file.exists(outputs_path[1]))
  unlink(outputs_path)
  expect_false(file.exists(outputs_path[1]))
  toc()

}


#### End of code.
#########################
#########################
