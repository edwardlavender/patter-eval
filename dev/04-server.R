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
# This is the set of directories that we expect to contain outputs (see here_alg()):
# batch <- 1:1181L
batch <- 1:1000L
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

#### Transfer log.txt files
# TO DO

#### Transfer success record (sdt) for batch
# TO DO

#### Transfer outputs
# Set up cluster
cl <- parallel::makeCluster(2L)
parallel::clusterExport(cl, "here_output")
# Copy files
tic()
success <-
  cl_lapply(outputs,
            .cl = cl,
            .fun = function(output) {
              file.copy(file.path(here_output("run"), output),
                        file.path("\\tsclient\patter-eval\data\sims\output\run", output),
                        overwrite = TRUE)
            })
toc()
# Validate success
table(success)
expect_true(all(success))

#### Clean up patter outputs (~5 s)
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


#### End of code.
#########################
#########################
