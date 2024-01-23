pf_backward_sampler <- function(.history,
                                .dens_step = dstep, ...,
                                .save_history = FALSE, .write_history = NULL,
                                .cl = NULL, .cl_varlist = NULL, .cl_chunks = TRUE,
                                .verbose = TRUE
) {

  #### Check user inputs
  t_onset <- Sys.time()
  check_inherits(.history, "list")
  if (!.save_history && is.null(.write_history)) {
    abort("`.save_history = FALSE` and `.write_history = NULL`. There is nothing to do.")
  }
  write_history_folder <- .pf_check_write_history(.write_history)

  #### Set up messages (as usual)
  cat_log <- cat_init(.verbose = .verbose)
  cat_log(call_start(.fun = "pf_backward_sampler", .start = t_onset))
  on.exit(cat_log(call_end(.fun = "pf_backward_sampler", .start = t_onset, .end = Sys.time())), add = TRUE)

  #### Set up loop
  cat_log("... Set up...")
  # Define whether or not to read history files
  if (inherits(.history[[1]], "data.frame")) {
    read_history <- FALSE
  } else {
    read_history <- TRUE
  }
  # Define constants
  n_step     <- length(.history)
  if (read_history) {
    .history[[n_step]] <- arrow::read_parquet(.history[[n_step]])
  }
  n_particle <- fnrow(.history[[n_step]])
  density    <- NULL

  #### Generate paths
  cat_log("... Generating path(s)...")
  paths <-
    cl_lapply(seq_len(n_particle),
              .cl = .cl,
              .varlist = .cl_varlist,
              .use_chunks = .cl_chunks,
              .fun = function(i) {

                #### Set up loop
                cat_log(paste("... ... On particle", i, "..."))
                cat_log("... ... ... Preparing to run sampler...")
                path <- dens <- list()
                path[[n_step]] <- .history[[n_step]][i, ]
                path[[n_step]][, density := 1]

                #### Run backwards sampler for a selected particle (i)
                cat_log("... ... ... Running sampler...")
                for (t in n_step:2) {
                  # Read history if necessary
                  cat_log(paste("... ... ... ... On time step", t, "..."))
                  if (read_history) {
                    # Drop history for t (to save memory) & read new history
                    .history[[t]] <- NA
                    .history[[t - 1]] <- arrow::read_parquet(.history[[t - 1]])
                  }
                  # Calculate step densities
                  dens <- .dens_step(.data_now = path[[t]],
                                     .data_past = .history[[t - 1]], ...)
                  # Sample a previous location
                  index <- sample.int(fnrow(.history[[t - 1]]), size = 1, prob = dens)
                  path[[t - 1]] <- .history[[t - 1]][index, ]
                  path[[t - 1]][, density := dens[index]]
                }

                #### Collate path
                cat_log("... ... ... Collating paths...")
                path <-
                  path |>
                  rbindlist() |>
                  mutate(path_id = i, .before = 1L) |>
                  as.data.table()

                #### Save path
                cat_log("... ... ... Recording path...")
                if (!is.null(.write_history)) {
                  .write_history$x    <- path
                  .write_history$sink <- file.path(write_history_folder, paste0(t, ".parquet"))
                  do.call(arrow::write_parquet, .write_history)
                }
                # Save path in memory
                if (.save_history) {
                  return(path)
                } else {
                  return(NULL)
                }

              })

  #### Return outputs
  cat_log("... Completing simulation...")
  if (is.null(paths[[1]])) {
    paths <- NULL
  } else {
    paths <-
      paths |>
      rbindlist() |>
      arrange(.data$path_id, .data$timestep) |>
      as.data.table()
  }
  t_end <- Sys.time()
  time <- list(start = t_onset,
               end = t_onset,
               duration = difftime(t_end, t_onset))
  out <- list(path = paths,
              time = time)

  #### Return outputs
  class(out) <- c(class(out), "pf_path")
  out
}
