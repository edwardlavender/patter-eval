# Try julia_connect() iteratively

# Running Julia in parallel on a socket cluster is tricky:
# * In Julia, precompile all Julia dependencies (used by JuliaCall, Patter.jl, patter)
# - This is implemented as precompilation in Julia happens in parallel which may interfere with the R process
# * Run example pf_filter() code normally (in case additional R setup is necessary)
# * Run RStudio/R as an administrator
# * Use julia_connect() before parallelisation
# * Use try_julia_connect() to connect to Julia in clusterEvalQ (incl. a retry mechanism)
# * Use JuliaCall:::.julia$cmd("using RCall") to set RCall
# * (optional) Try moving julia_connect() within the parallel code

# See also: https://github.com/Non-Contradiction/JuliaCall/issues/120

try_julia_connect <- function() {
  attempt <- 1
  while (attempt < 10L) {
    connection <- tryCatch(patter::julia_connect(.threads = 1L),
                           error = function(e) e)
    if (inherits(connection, "error")) {
      attempt <- attempt + 1
      Sys.sleep(60)
    } else {
      attempt <- Inf
    }
  }
  invisible(connection)
}
