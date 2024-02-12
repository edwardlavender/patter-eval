#' @title Internal {patter} functions
# utils.add::load_internal_functions("patter")
check_dlist <- patter:::check_dlist
dist_2d     <- patter:::dist_2d
normalise   <- patter:::normalise
fnrow       <- collapse::fnrow

#' @title {patter} helpers

acs_setup_receiver_key <-
  function(.receiver) {
    paste(sort(unlist(.receiver)), collapse = "-")
  }

#' @title Proposal densities

pf_dpropose_read <- function(.particles, .obs, .t, .dlist, .drop) {
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
  .particles[, dens := spat_dens[match(.particles$cell_past, spat_dens$cell)]$dens]
  if (.drop) {
    .particles <- .particles[which(dens > 0), ]
  }
  .particles
}

#' @title (Deprecated) A streamlined pf_backward_sampler_p()
#' @description This is a streamlined version of `pf_backward_sampler_p()` that was implemented to test potential speed ups for sensitivity analyses. The speed of this function was too slow and it is not used.
#'
#' @param .history A named list of data.tables. Each data.table must contain a cell_now column.
#' @param .dlist A named list with a .dlist$algorithm$sim$combination element.
#' @param .read Logical; read densities from file or compute on the fly. Implemented to compare the speed of the two approaches.
#'
#' @details (optional) TO DO Note that this implementation does not currently varying numbers of particles or irregular resampling.

pf_backward_sampler_p_fst <- function(.history, .dlist, .read = TRUE) {

  n_step     <- length(.history)
  n_particle <- fnrow(.history[[n_step]])

  #### Build paths
  paths <-
    pbapply::pblapply(
      seq_len(n_particle),
      function(i) {

        #### Initiate loop over time steps
        path <- vector("list", length = n_step)
        if (.read) {
          path[[n_step]] <- .history[[n_step]][i, .(timestep, cell_now)]
        } else {
          path[[n_step]] <- .history[[n_step]][i, .(timestep, cell_now, x_now, y_now)]
        }

        #### Run backwards sampler for a selected particle (i)
        for (t in n_step:2) {

          # Define previous time step
          # print(glue::glue("Particle {i} (timestep {t}"))
          tp <- t - 1L

          # Read or calculate step densities
          if (.read) {
            spat_dens <-
              here_input("density",
                         .dlist$algorithm$sim$combination,
                         paste0(path[[t]]$cell_now[1], ".qs")) |>
              qs::qread()
            dens <- spat_dens[fastmatch::fmatch(.history[[tp]]$cell_now, spat_dens$cell)]$dens
            dens[is.na(dens)] <- 0
          } else {
            len <- dist_2d(.x0 = path[[t]]$x_now,
                           .y0 = path[[t]]$y_now,
                           .x1 = .history[[tp]]$x_now,
                           .y1 = .history[[tp]]$y_now)
            dens <- dtruncgamma(len,
                                .shape = .dlist$algorithm$sim$shape,
                                .scale = .dlist$algorithm$sim$scale,
                                .mobility = .dlist$algorithm$sim$mobility)
          }

          # Sample a previous location
          # TO DO
          # * Check normalise densities
          index <- sample.int(length(dens), size = 1L, prob = dens)
          # Record previous location
          if (.read) {
            path[[tp]] <- data.table(timestep = tp,
                                     cell_now = .history[[tp]]$cell_now[index])
          } else {
            path[[tp]] <- data.table(timestep = tp,
                                     cell_now = .history[[tp]]$cell_now[index],
                                     x_now = .history[[tp]]$x_now[index],
                                     y_now = .history[[tp]]$y_now[index])
          }

        }

        #### Collate & return path
        path |>
          rbindlist() |>
          mutate(path_id = i, .before = 1L) |>
          as.data.table()

      })

  #### Return list of paths
  paths

}
