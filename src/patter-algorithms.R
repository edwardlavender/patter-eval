#' @title Internal {patter} functions
# utils.add::load_internal_functions("patter")
check_dlist         <- patter:::check_dlist
cat_init            <- patter:::cat_init
call_start          <- patter:::call_start
call_end            <- patter:::call_end
call_time           <- patter:::call_time
.pf_history_list    <- patter:::.pf_history_list
.pf_history_read    <- patter:::.pf_history_read
.pf_history_cols    <- patter:::.pf_history_cols
.pf_history_elm     <- patter:::.pf_history_elm
.pf_history_write   <- patter:::.pf_history_write
.pf_write_particles <- patter:::.pf_write_particles

fnrow <- collapse::fnrow

#' @title Proposal densities

pf_dpropose_read <- function(.particles, .obs, .t, .dlist) {
  # Check inputs
  if (.t == max(.obs$timestep)) {
    check_dlist(.dlist = .dlist,
                .algorithm = "sim")
  }
  # Read precomputed densities from file
  spat_dens <-
    here_input("density",
               .dlist$algorithm$sim$combination,
               paste0(.particles$cell_now[1], ".qs")) |>
    qs::qread()
  # Match densities & return .particles with a 'dens' column
  # * This behaviour matches pf_dpropose()
  .particles[, dens := spat_dens[match(.particles$cell_past, spat_dens$cell)]]
  .particles[dens > 0, ]
}

#' @title Backward sampler (parallel implementation)

pf_backward_sampler_p <- function(.history,
                                  .dpropose = pf_dpropose,
                                  .obs = NULL,
                                  .dlist,
                                  .dargs = list(),
                                  .cl = NULL,
                                  .cl_varlist = NULL,
                                  .cl_chunk = cl_chunk(.cl),
                                  .record = pf_opt_record(),
                                  .verbose = getOption("patter.verbose")
) {

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
        path[[n_step]] <- .history[[n_step]][i, ]
        path[[n_step]][, dens := 1]

        #### Run backwards sampler for a selected particle (i)
        cat_log("... ... ... Running sampler...")
        for (t in n_step:2) {
          # Read history if necessary
          cat_log(paste("... ... ... ... On time step", t, "..."))
          tp <- t - 1L
          if (read) {
            # Drop history for t (to save memory) & read new history
            .history[[t]] <- NA
            .history[[tp]] <- .pf_history_elm(.history = .history, .elm = tp,
                                              cols = read_cols)
          }
          # Calculate step densities
          pnow <-
            dplyr::bind_cols(
              path[[t]] |>
                select("cell_now", "x_now", "y_now") |>
                as.data.table(),
              .history[[tp]] |>
                select(cell_past = "cell_now",
                       x_past = "x_now",
                       y_past = "y_now") |>
                as.data.table()
            ) |>
            as.data.table()
          .dargs$.particles <- pnow
          .dargs$.t         <- t
          prob              <- do.call(.dpropose, .dargs)$dens
          # Sample a previous location
          index <- sample.int(length(prob), size = 1L, prob = prob)
          path[[tp]] <- .history[[tp]][index, ]
          path[[t]][, dens := prob[index]]
        }

        #### Collate path
        cat_log("... ... ... Collating paths...")
        path[[1]][, dens := NA_real_]
        path <-
          path |>
          rbindlist() |>
          mutate(path_id = i, .before = 1L) |>
          as.data.table()

        #### Save path
        cat_log("... ... ... Recording path...")
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
