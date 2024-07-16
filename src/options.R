# Set global options
op <- options()
# Set terra colour palette
options(terra.pal = rev(terrain.colors(256)))
# (optional) On lavended, beep on error
if (FALSE & Sys.info()["user"] == "lavended") {
  options(error = function(...) beepr::beep(7))
} else {
  options(error = NULL)
}

# Progress bar (chunk) options for cl_lapply()
if (Sys.info()["user"] == "lavended") {
  pbapply::pboptions("nout" = 10)
} else {
  pbapply::pboptions("nout" = 50)
}
