#' @title {patter} internals

normalise           <- patter:::normalise
fnrow               <- collapse::fnrow

#' @title (Deprecated) Streamlined pf_backward_sampler_p() function
#' This is a streamlined version of pf_backward_sampler_p() that was implemented to test potential speed ups for sensitivity analyses.
#' Function speed was still too slow & this function is unused.
#'
#' @param .history A named list of data.tables. Each data.table must contain a cell_now column.
#' @param .dlist A named list with a .dlist$algorithm$sim$combination element.
#' @details TO DO. Update to handle varying numbers of particles & infrequent resampling, if required

pf_backward_sampler_p_fst <- function(.history, .dlist) {

  n_step     <- length(.history)
  n_particle <- fnrow(.history[[n_step]])

  #### Build paths
  paths <-
    pbapply::pblapply(
      seq_len(n_particle),
      function(i) {

        #### Initiate loop over time steps
        path <- vector("list", length = n_step)
        path[[n_step]] <- .history[[n_step]][i, .(timestep, cell_now)]

        #### Run backwards sampler for a selected particle (i)
        for (t in n_step:2) {

          # Define previous time step
          tp <- t - 1L

          # Read or calculate step densities
          spat_dens <-
            here_input("density",
                       .dlist$algorithm$sim$combination,
                       paste0(path[[t]]$cell_now[1], ".qs")) |>
            qs::qread()
          dens <- spat_dens[fastmatch::fmatch(.history[[tp]]$cell_now, spat_dens$cell)]$dens
          dens[is.na(dens)] <- 0
          # Sample a previous location
          index <- sample.int(length(dens), size = 1L, prob = normalise(dens))
          # Record previous location
          path[[tp]] <- data.table(timestep = tp,
                                   cell_now = .history[[tp]]$cell_now[index])
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
