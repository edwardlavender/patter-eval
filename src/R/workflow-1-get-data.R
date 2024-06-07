#' @title Data preparation functions

get_array <- function(sim, arrays) {
  arrays[[sim$array_type]] |>
    filter(array_id == sim$array_realisation) |>
    as.data.table()
}

get_acoustics <- function(sim, acoustics) {
  acoustics[[sim$combination]][[sim$array_type]] |>
    filter(array_id == sim$array_realisation) |>
    filter(path_id == sim$path_realisation) |>
    as.data.table()
}

get_detections <- function(sim, detections) {
  get_acoustics(sim, detections)
}


get_path <- function(sim, paths, detections) {
  paths[[sim$path_type]] |>
    filter(path_id == sim$path_realisation) |>
    filter(timestamp >= min(detections$timestamp) &
             timestamp <= max(detections$timestamp)) |>
    as.data.table()
}
