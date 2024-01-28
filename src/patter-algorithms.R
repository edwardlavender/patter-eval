#' @title Internal {patter} functions
# utils.add::load_internal_functions("patter")
check_dlist         <- patter:::check_dlist

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
