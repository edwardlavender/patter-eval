#' @title Streamlined pf_backward_sampler_p() function

# TO DO
# * Update to handle varying numbers of particles & infrequent resampling, if required

pf_backward_sampler_p_fst <- function(.history,
                                      .dpropose = pf_dpropose,
                                      .obs = NULL,
                                      .dlist,
                                      .dargs = list(),
                                      .control = pf_opt_control(),
                                      .cl = NULL,
                                      .cl_varlist = NULL,
                                      .cl_chunk = cl_chunk(.cl),
                                      .record = pf_opt_record(),
                                      .verbose = getOption("patter.verbose")) {

  #### Check user inputs
  # TO DO

  #### Set up messages
  t_onset <- Sys.time()
  cat_log <- cat_init(.verbose = .verbose)
  cat_log(call_start(.fun = "pf_backward_sampler_p", .start = t_onset))
  on.exit(cat_log(call_end(.fun = "pf_backward_sampler_p", .start = t_onset,
                           .end = Sys.time())), add = TRUE)
  .history <- .pf_history_list(.history)
  read     <- .pf_history_read(.history)

  #### Set up loop
  # .history
  .history  <- .pf_history_list(.history)
  read      <- .pf_history_read(.history)
  inout     <- .pf_history_cols(.history = .history, .record = .record,
                                .input_cols = c("cell_now", "x_now", "y_now"))
  .record   <- inout$.record
  read_cols <- inout$read_cols
  write     <- .pf_history_write(.record)
  # Final (starting) particle samples for the backward sampler
  n_step <- length(.history)
  .history[[n_step]] <- .pf_history_elm(.history = .history, .elm = n_step,
                                        cols = read_cols)
  n_particle <- fnrow(.history[[n_step]])
  # Function arguments
  .dargs$.obs   <- .obs
  .dargs$.dlist <- .dlist
  .dargs$.drop  <- .control$drop
  # Global variables
  dens <- NULL

  #### Generate paths
  cat_log("... Generating path(s)...")
  paths <-
    cl_lapply(
      seq_len(n_particle),
      .cl = .cl,
      .varlist = .cl_varlist,
      .chunk = .cl_chunk,
      .fun = function(i) {

        #### Set up loop
        cat_log(paste("... ... On particle", i, "..."))
        cat_log("... ... ... Preparing to run sampler...")
        path <- density <- list()
        path[[n_step]] <-
          .history[[n_step]][i, ] |>
          select("timestep", "cell_now") |>
          as.data.table()

        #### Run backwards sampler for a selected particle (i)
        cat_log("... ... ... Running sampler...")
        for (t in n_step:2) {

          # Read history if necessary
          cat_log(paste("... ... ... ... On time step", t, "..."))
          tp <- t - 1L
          # if (read) {
          #   # Drop history for t (to save memory) & read new history
          #   .history[[t]] <- NA
          #   .history[[tp]] <- .pf_history_elm(.history = .history, .elm = tp,
          #                                     .read = TRUE, cols = read_cols)
          # }
          # Calculate step densities (cell_now -> cell_past)
          spat_dens <-
            here_input("density",
                       .dlist$algorithm$sim$combination,
                       paste0(path[[t]]$cell_now[1], ".qs")) |>
            qs::qread()
          dens <- spat_dens[match(.history[[tp]]$cell_now, spat_dens$cell)]$dens
          # Sample a previous location
          index <- sample.int(length(dens), size = 1L, prob = normalise(dens))
          # Record previous location
          path[[tp]] <- data.table(timestep = tp,
                                   cell_now = .history[[tp]]$cell_now[index])
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
        if (!is.null(.record$cols)) {
          path <- path |> select(all_of(.record$cols)) |> as.data.table()
        }
        .pf_write_particles(.particles = path, .sink = .record$sink,
                            .filename = t, .write = write)
        # Save path in memory
        if (.record$save) {
          return(path)
        } else {
          return(NULL)
        }
      })

  #### Return outputs
  if (is.null(paths[[1]])) {
    paths <- NULL
  } else {
    paths <-
      paths |>
      rbindlist() |>
      arrange(.data$path_id, .data$timestep) |>
      as.data.table()
  }
  .pf_backward_output(.start = t_onset,
                      .history = list(),
                      .path = paths,
                      .record = .record)
}
