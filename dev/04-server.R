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
library(patter)
library(testthat)
library(tictoc)
dv::src()


#########################
#########################
#### Management

#### (conditional) Clean up old files on server
# TO DO

#### List patter outputs (~5 s)
tic()
logs    <- list.files(here_output("log", "patter", "performance"),
                      full.names = TRUE)
outputs <- list.files(here_output("run"), recursive = TRUE)
outputs <- outputs[stringr::str_detect(outputs, "patter")]
toc()

#### (optional) Transfer log files
# TO DO

#### Transfer success record (sdt)
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
expect_true(file.exists(logs[1]))
unlink(logs)
expect_false(file.exists(logs[1]))
expect_true(file.exists(file.path(here_output("run"), outputs[1])))
unlink(file.path(here_output("run"), outputs))
expect_false(file.exists(file.path(here_output("run"), outputs[1])))
toc()


#### End of code.
#########################
#########################
