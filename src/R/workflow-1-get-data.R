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
  dets <- get_detections(sim, detections)
  paths[[sim$path_type]] |>
    filter(path_id == sim$path_realisation) |>
    filter(timestamp >= min(dets$timestamp) &
             timestamp <= max(dets$timestamp)) |>
    as.data.table()
}

get_archival <- function(sim, paths, detections) {
  path <- get_path(sim, paths, detections)
  path |>
    mutate(obs = depth,
           sensor_id = 1L,
           depth_shallow_eps = 5,
           depth_deep_eps = 5) |>
    select(timestamp, sensor_id, obs, depth_shallow_eps, depth_deep_eps) |>
    as.data.table()
}
