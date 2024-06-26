# Set global options
op <- options()
# Set terra colour palette
options(terra.pal = rev(terrain.colors(256)))
# On lavended, beep on error
if (Sys.info()["user"] == "lavended") {
  options(error = function(...) beepr::beep(7))
}

# Progress bar (chunk) options for cl_lapply()
if (Sys.info()["user"] == "lavended") {
  pbapply::pboptions("nout" = 10)
} else {
  pbapply::pboptions("nout" = 50)
}
