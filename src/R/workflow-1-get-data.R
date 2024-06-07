#' @title Data preparation functions

get_array <- function(sim, arrays) {
  arrays[[sim$array_type]] |>
    filter(array_id == sim$array_realisation) |>
    as.data.table()
}

get_acoustics <- function(sim, detections) {
  detections[[sim$combination]][[sim$array_type]] |>
    filter(array_id == sim$array_realisation) |>
    filter(path_id == sim$path_realisation) |>
    # Focus on detections only
    filter(obs == 1L) |>
    as.data.table()
}

get_path <- function(sim, paths, acoustics) {
  paths[[sim$path_type]] |>
    filter(path_id == sim$path_realisation) |>
    filter(timestamp >= min(acoustics$timestamp) &
             timestamp <= max(acoustics$timestamp)) |>
    as.data.table()
}
