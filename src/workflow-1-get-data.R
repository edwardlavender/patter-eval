#' @title Data preparation functions

get_array <- function(sim, arrays) {
  arrays[[sim$array_type]] |>
    filter(array_id == sim$array_realisation) |>
    as.data.table()
}

get_detections <- function(sim, detections) {
  # This is used for RSP analyses (performance sims only)
  # So we keep the detection probability parameters used to simulate data
  detections[[sim$combination]][[sim$array_type]] |>
    filter(array_id == sim$array_realisation) |>
    filter(path_id == sim$path_realisation) |>
    filter(obs == 1L) |>
    as.data.table()
}

get_acoustics <- function(sim, acoustics) {
  # Define the full acoustic time series
  acoustics <-
    acoustics[[sim$combination]][[sim$array_type]] |>
    filter(array_id == sim$array_realisation) |>
    filter(path_id == sim$path_realisation) |>
    select("timestamp", "obs", "sensor_id", "receiver_x", "receiver_y") |>
    as.data.table()
  # Use detection probability parameters from sim
  # * This is required b/c the parameters used to simulate data
  # * ... differ from those used to model data for sensitivity simulations
  acoustics[, receiver_alpha := sim$alpha]
  acoustics[, receiver_beta := sim$beta]
  acoustics[, receiver_gamma := sim$gamma]
  # Identify the detections
  dets <-
    acoustics |>
    filter(obs == 1L) |>
    as.data.table()
  # Focus on the portion of acoustic data between the 1st & last detection
  acoustics |>
    filter(timestamp >= min(dets$timestamp) & timestamp <= max(dets$timestamp)) |>
    as.data.table()
}

get_path <- function(sim, paths, detections) {
  # Get detections
  dets <- get_detections(sim, detections)
  # Focus on the portion of the path between the 1st and last detections
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
