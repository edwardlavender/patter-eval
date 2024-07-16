# Try julia_connect() iteratively

# Running Julia in parallel on a socket cluster is tricky:
# * In Julia, precompile all Julia dependencies (used by JuliaCall, Patter.jl, patter)
# -- This is implemented as precompilation in Julia happens in parallel which may interfere with the R process
# -- See GitHub thread below
# * Run example pf_filter() code normally (in case additional R setup is necessary)
# * Run RStudio/R as an administrator
# * Use julia_connect() before parallelisation
# -- See GitHub thread below
# * Use try_julia_connect() to connect to Julia in clusterEvalQ
# - (This includes a retry mechanism if the Julia process is locked)
# * Use JuliaCall:::.julia$cmd("using RCall") to set RCall
# -- See GitHub thread below
# * (optional) Try moving julia_connect() within the parallel code

# See also: https://github.com/Non-Contradiction/JuliaCall/issues/120

try_julia_connect <- function() {
  attempt <- 1L
  n_trial <- 10L
  # Iteratively attempt to connect to Julia
  # (Julia may be locked by another process)
  while (attempt <= n_trial) {
    connection <- tryCatch(patter::julia_connect(.threads = 1L),
                           error = function(e) e)
    if (inherits(connection, "error")) {
      if (attempt == n_trial) {
        stop(paste("`julia_connect()` failed after multiple attempts with the following error message: \n\n",
                   connection$message), call. = FALSE)
      }
      attempt <- attempt + 1L
      warning(paste("`julia_connect()` failed on ", attempt, " / ", n_trial, ". Retrying after 30 s..."),
              immediate. = TRUE)
      Sys.sleep(30)
    } else {
      attempt <- Inf
    }
  }
  invisible(connection)
}
