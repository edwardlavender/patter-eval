#' @title Data generation/simulation helpers

gen_depth <- function(.xy) {
  250 - 100 * cos(sqrt(.xy$x^2 + .xy$y^2) / (500 * 2 * pi))
}

sim_depth <- function(.seabed) {
  z <- runif(length(.seabed), .seabed - 5, .seabed + 5)
  z[z < 0] <- 0
  z
}

#' @title Filter simulations to do
#' Given a data.table of simulations, this function identifies the set remaining without complete outputs that need to be run.

do_sims <- function(sims) {

  N       <- nrow(sims)
  sims    <- copy(sims)
  sims_ls <- split(sims, collapse::seq_row(sims))

  # Define 'done' simulations (TRUE, FALSE)
  done    <- pbapply::pblapply(sims_ls, function(sim) {

    # Define path to convergence.rds
    acpf_convergence.rds   <- here_alg(sim, "patter", "acpf", sim$alg_par, "convergence.rds")
    acdcpf_convergence.rds <- here_alg(sim, "patter", "acdcpf", sim$alg_par, "convergence.rds")

    # If either file does not exist, the simulation is not finished
    if (!file.exists(acpf_convergence.rds) | !file.exists(acdcpf_convergence.rds)) {
      return(FALSE)
    }

    # Read convergence
    acpf_convergence   <- readRDS(acpf_convergence.rds)
    acdcpf_convergence <- readRDS(acdcpf_convergence.rds)

    # If both convergence files report FALSE, the simulation is finished
    if (isFALSE(acpf_convergence) & isFALSE(acdcpf_convergence)) {
      return(TRUE)
    }

    # Define ud-s.tif files
    acpf.tif   <- here_alg(sim, "patter", "acpf", sim$alg_par, "ud-s.tif")
    acdcpf.tif <- here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-s.tif")

    # If both .tif files exist, the simulation is finished
    if (file.exists(acpf.tif) & file.exists(acdcpf.tif)) {
      return(TRUE)
    }

    # Otherwise, the simulation is not finished
    return(FALSE)

  }) |> unlist()

  # Return simulations to do (i.e., done = FALSE)
  pos <- which(done == FALSE)
  print(paste(length(pos), "/", N, "simulation(s) to do..."))
  sims[pos, ]

}
