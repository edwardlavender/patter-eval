# Progress bar (chunk) options for cl_lapply()
pbapply::pboptions("nout" = 10)

# Set options
op <- options()
# Set terra colour palette
options(terra.pal = rev(terrain.colors(256)))
# On lavended, beep on error
if (Sys.info()["user"] == "lavended") {
  options(error = function(...) beepr::beep(7))
}


