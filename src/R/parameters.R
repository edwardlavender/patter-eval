if (Sys.info()["nodename"] == "siam-linux20") {
  pkg_config <-
    '
  # Check & remove path set by R
  println(ENV["LD_LIBRARY_PATH"])
  delete!(ENV, "LD_LIBRARY_PATH")
  # Set default Julia path
  # * This is obtained by running Julia directly from the terminal
  # * using Libdl
  # * filter!(contains("curl"), dllist())
  ENV["LD_LIBRARY_PATH"] = "/opt/julias/julia-1.10/bin/../lib/julia/libcurl.so.4"
  # Validate settings
  println(ENV["LD_LIBRARY_PATH"])
  '
} else {
  pkg_config <- NULL
}

crs <- "+proj=utm +zone=1 +datum=WGS84"
sr  <- 10
