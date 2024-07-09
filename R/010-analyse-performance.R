#########################
#########################
#### analyse-performance.R

#### Aims:
# (1) Compares algorithm performance

#### Prerequisites
# 1) Implement simulations


#########################
#########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(collapse)
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(prettyGraphics)
library(testthat)
library(tictoc)
dv::src()

#### Load data
# Performance simulations
sims   <- readRDS(here_input("sims-performance.rds"))
skills <- readRDS(here_output("synthesis", "skill.rds"))


#########################
#########################
#### Computation time

#### Define a data.table of wall times (~2 s)
# Times are for single-threaded mode with parallelisation from R
# Times (for performance sims) recorded on SIA-LAVENDED
time <- cl_lapply(.x = split(sims, seq_row(sims)),
                  .fun = function(sim) {

                    lapply(c("acpf", "acdcpf"), function(alg) {

                      # Define convergence.rds & time.qs
                      file_convergence <-
                        here_alg(sim, "patter", alg, sim$alg_par, "convergence.rds")
                      file_time <-
                        here_alg(sim, "patter", alg, sim$alg_par, "time.qs")

                      # Validate file paths
                      # expect_true(file.exists(file_convergence))
                      # expect_true(file.exists(file_time))

                      # Read files
                      if (file.exists(file_time)) {
                        time <- qs::qread(file_time)
                      } else {
                        time <- data.table(id = sim$id,
                                           pff = NA_real_,
                                           pfbs = NA_real_,
                                           udf = NA_real_,
                                           uds = NA_real_,
                                           n = NA_integer_,
                                           success = FALSE)
                      }
                      time[, n := nrow(read_path(sim))]
                      time[, success := FALSE]
                      if (file.exists(file_convergence)) {
                        time[, success := readRDS(file_convergence)]
                      }
                      time[, algorithm := alg]
                      time

                    }) |> rbindlist()


                  }, .combine = rbindlist)

#### Summary statistics
head(time)
nrow(time)
# Convergence
time[success == FALSE, ]
# Forward filter time (per time step) in seconds
(time[, pff / n] * 60) |> utils.add::basic_stats(na.rm = TRUE, p = NULL)
# Smoothing time (per time step) in seconds
(time[, pfbs / n] * 60) |> utils.add::basic_stats(na.rm = TRUE, p = NULL)

#### Visualise wall time
# Define total wall time
ggtime_total <-
  time |>
  tidyr::pivot_longer(cols = c(pff, pfbs), names_to = "routine", values_to = "time") |>
  mutate(time = time * 60,
         Algorithm = factor(algorithm, levels = c("acpf", "acdcpf"), labels = c("ACPF", "ACDCPF")),
         routine = factor(routine, levels = c("pff", "pfbs"), labels = c("Forward filter", "Smoother")),
         statistic = "(A) Total wall time") |>
  as.data.table()
# Define wall time per time step
ggtime_step <- copy(ggtime_total)
ggtime_step[, statistic := "(B) Wall time per time step"]
ggtime_step[, time := time / n]
# Summarise wall time per time step
ggtime_step |>
  filter(routine == "Forward filter") |>
  group_by(algorithm) |>
  reframe(mean(time, na.rm = TRUE))
ggtime_step |>
  filter(routine == "Smoother") |>
  summarise(mean(time, na.rm = TRUE))
# Collect wall time statistics
ggtime <- rbind(ggtime_total, ggtime_step)
# Visualise wall time statistics
png(here_fig("performance", "wall-time.png"),
    height = 4, width = 12, units = "in", res = 600)
ggplot(ggtime) +
  geom_violin(aes(routine, time, fill = Algorithm)) +
  facet_wrap(~statistic, scales = "free") +
  xlab("Routine") + ylab("Time (s)") +
  scale_y_continuous(expand = c(0, 0)) +
  ggthemes::theme_clean() +
  theme(strip.text = element_text(face = "bold"))
dev.off()


#########################
#########################
#### Figure preparation

#### Define algorithms
# Define algorithms (subset)
algs <- data.frame(alg = c("Null", "COA(30)", "COA(120)", "RSP(1)", "RSP(2)", "ACPF(F)", "ACDCPF(F)", "ACPF(S)", "ACDCPF(S)"),
                   alg_name = c("Null", "COA", "COA", "RSP", "RSP", "ACPF", "ACDCPF", "ACPF", "ACDCPF"),
                   col_name = c("dimgrey", "red", "darkred", "orange", "darkorange", "slateblue1", "slateblue4", "lightgreen", "darkgreen")
                   )
algs$col <- scales::alpha(algs$col_name, 0.5)
algs$label <- seq_len(nrow(algs))
skills <-
  skills |>
  filter(alg %in% algs$alg) |>
  left_join(algs, by = "alg", relationship = "many-to-one") |>
  as.data.table()

#### Define selected receiver numbers
# On some plots, we only consider two different numbers of receivers
nr <- c(10, 100)

#### Define skill metrics & constant limits across all plots
metrics <- c("mb", "me", "rmse", "R", "d")
metrics_lim <-
  list(mb = range(skills$mb, na.rm = TRUE),
       me = range(skills$me, na.rm = TRUE),
       rmse = range(skills$rmse, na.rm = TRUE),
       R = c(0, 1),
       d = c(0, 1))

#### Define combination
# We create all plots for each hypothetical study system
# (i.e., combination of movement and detection probability parameters)
# Update combination == 1L or combination == 2L here to create plots for different study systems
# png_name() handles file names accordingly
comb <- 1L
unique(skills$combination)
skills_all <- copy(skills)
skills <-
  skills |>
  filter(combination == comb) |>
  as.data.table()
sims <-
  sims |>
  filter(combination == comb) |>
  as.data.table()
# Define function for combination-specific names
png_name <- function(name) {
  stopifnot(length(unique(skills$combination)) == 1L)
  paste0(name, "-",  skills$combination[1], ".png")
}


#########################
#########################
#### Maps of space use

#### Select simulations
# For specific combination, select two receiver levels, one array/path realisation
# This choice is informed by examination of the maps below
# Good examples: 1, 3, 10, 24, 29
sims_for_maps <-
  sims |>
  filter(n_receiver %in% nr) |>
  filter(array_realisation == 1L & path_realisation == 3L) |>
  arrange(n_receiver, arrangement) |>
  as.data.table()

#### Make maps (~5 mins)
if (FALSE) {

  # Create a map for each path_realisation_index
  tic()
  lapply(sort(unique(sims$path_realisation)), function(path_realisation_index) {

    #### Define simulations for path_realisation_index
    sims_for_maps <-
      sims |>
      filter(n_receiver %in% nr) |>
      filter(array_realisation == 1L & path_realisation == path_realisation_index) |>
      arrange(n_receiver, arrangement) |>
      as.data.table()

    #### Set up plot
    # The figure will be annotated outside of R
    # (optional) FLAG: adjust width & number of columns for algorithms
    here_maps <- function(...) here_fig("performance", "maps", ...)
    dir.create(here_maps(), recursive = TRUE)
    png(here_maps(png_name(paste0("map-", path_realisation_index))),
        height = 8, width = 12, units = "in", res = 600)
    pp <- par(mfrow = c(4, 9),
              oma = c(0, 0, 0, 0), mar = rep(1, 4))
    pbapply::pblapply(1:4, function(i) {

      #### Read UDs
      sim        <- sims_for_maps[i, ]
      if (all(is.na(sim))) {
        return(NULL)
      }
      # list.files(here_alg(sim), recursive = TRUE)
      ud_path    <- read_rast(here_alg(sim, "path", "ud.tif"))
      ud_coa_1   <- read_rast(here_alg(sim, "coa", "30 mins", "ud.tif"))
      ud_coa_2   <- read_rast(here_alg(sim, "coa", "120 mins", "ud.tif"))
      ud_rsp_1   <- read_rast(here_alg(sim, "rsp", "default", "ud.tif"))
      ud_rsp_2   <- read_rast(here_alg(sim, "rsp", "custom", "ud.tif"))
      ud_acpff   <- read_rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-f.tif"))
      ud_acdcpff <- read_rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-f.tif"))
      ud_acpfs   <- read_rast(here_alg(sim, "patter", "acpf", sim$alg_par, "ud-s.tif"))
      ud_acdcpfs <- read_rast(here_alg(sim, "patter", "acdcpf", sim$alg_par, "ud-s.tif"))
      # Select UDs
      # * (optional) FLAG: modify UDs included in this list
      uds <- list(ud_path,
                  ud_coa_1, ud_coa_2,
                  ud_rsp_1, ud_rsp_2,
                  ud_acpff, ud_acdcpff,
                  ud_acpfs, ud_acdcpfs
      )
      if (any(sapply(uds, is.null))) {
        return(NULL)
      }

      #### Calculate scaling parameter
      # We will scale UDs (within rows) to a maximum value of 1 for comparison
      # Scaling is implemented below to ensure correct HR calculation
      scale <- sapply(uds, \(ud) terra::global(ud, "max")$max) |> max()
      # uds   <- lapply(uds, \(ud) ud / scale)
      # sapply(uds, \(ud) terra::global(ud, "max"))

      #### Plot UDs for selected array
      m <- read_array(sim)
      lapply(uds, function(ud) {
        # Plot UD (scaled)
        rescale <- FALSE
        if (rescale) {
          sud <- ud / scale
          range <- c(0, 1)
        } else {
          sud <- ud
          range <- NULL
        }
        spatMap(sud, range = range, legend = FALSE, mar = NA)
        # (optional) Set speed to TRUE to check plot layout only
        speed <- FALSE
        if (speed) {
          spatAxes(ud)
          return(NULL)
        }
        # Add home range (based on unscaled UD)
        map_hr_home(ud, .add = TRUE,
                    border = "dimgrey", lwd = 0.75)
        # Plot UD (scaled)
        # terra::plot(ud)
        # Add receivers
        points(m$receiver_x, m$receiver_y,
               pch = 21,
               bg = scales::alpha("black", 0.75),
               col = scales::alpha("black", 0.75),
               cex = 0.5)
        # Add detection range(s)
        # * We use a scale bar instead
        if (FALSE) {
          cbind(m$receiver_x, m$receiver_y) |>
            terra::vect() |>
            terra::buffer(width = sim$gamma) |>
            terra::lines(lwd = 0.5, col = "dimgrey")
        }
        spatAxes(ud)
        terra::sbar(sim$gamma, lonlat = FALSE, col = "darkred")
        NULL
      }) |> invisible()

    }) |> invisible()
    par(pp)
    dev.off()

  })
  toc()

}


#########################
#########################
#### Barplots of error statistics
# (selected barplots for maps)

#### Define plot type (full, all algorithms; partial, a subset)
# * full = all algorithms (as in the supporting map)
# * partial = a subset of algorithms (as in the main text figure)
type  <- "full"
width <- 3
type  <- "partial"
width <- 2.5

#### Define dataset
skills_for_bars <-
  skills |>
  filter(id %in% sims_for_maps$id) |>
  arrange(n_receiver, arrangement, label) |>
  as.data.table()
if (type == "partial") {
  skills_for_bars <-
    skills_for_bars |>
    filter(alg %in% c("COA(30)", "RSP(1)", "ACPF(S)", "ACDCPF(S)")) |>
    as.data.table()
}

#### Build figure
png(here_fig("performance", png_name(paste0("map-barplot-", type))),
    height = 8, width = width, units = "in", res = 600)
pp <- par(mfrow = c(4, 1),
          oma = c(2, 4, 2, 2),
          mar = rep(2, 4))
pbapply::pblapply(1:4, function(i) {

  # Isolate skill metrics
  sim <- sims_for_maps[i, ]
  sk <-
    skills_for_bars |>
    filter(id == sim$id) |>
    as.data.table()

  # Create blank plot
  sk$x <- seq_len(nrow(sk))
  axis_ls <- pretty_plot(sk$x,
                         sk$me,
                         xlim = c(0.5, max(sk$x) + 0.5),
                         ylim = metrics_lim[["me"]],
                         pretty_axis = list(axis = list(list(side = 1,
                                                             at = sk$x,
                                                             labels = sk$label),
                                                        list(NULL))
                         ),
                         xlab = "", ylab = "",
                         type = "n")

  # Add bars
  for (j in seq_len(nrow(sk))) {
    rect(xleft = sk$x[j] - 0.5,
         ybottom = axis_ls[[2]]$lim[1],
         xright = sk$x[j] + 0.5,
         ytop = sk$me[j],
         col = sk$col[j])

  }

}) |> invisible()
dev.off()


#########################
#########################
#### Boxplots of error statistics

#### Select simulations
# Define arrays
combs <-
  CJ(arrangement = c("random", "regular"),
     n_receiver = nr) |>
  arrange(n_receiver, arrangement)
# Define simulations
sims_for_skill <-
  sims |>
  filter(n_receiver %in% nr) |>
  arrange(n_receiver, arrangement) |>
  as.data.table()

#### Visualise boxplots
png(here_fig("performance", png_name("boxplots")),
    height = 10, width = 15, units = "in", res = 600)
pp <- par(mfrow = c(4, length(metrics)),
          oma = c(1, 3, 1, 1), mar = c(2, 2, 2, 2))
lapply(1:4, function(i) {
  # Define skill scores across path realisations
  skill <- skills[id %in%
                  sims_for_skill[n_receiver == combs$n_receiver[i] & arrangement == combs$arrangement[i], ]$id
                  ]

  # Visualise distribution of error scores (by metric)
  lapply(metrics, function(metric) {
    # metric <- metrics[1]
    x <- skill$label # skill$alg
    y <- skill[[metric]]
    # p <- min(floor(log10(abs(y))), na.rm = TRUE)
    # y <- y / 10^p
    pretty_boxplot(x, y,
                   varwidth = TRUE,
                   col = skill$col,
                   ylim = metrics_lim[[metric]],
                   pretty_axis_args = list(side = 1:2, pretty = list(list(n = 100), list(n = 5))),
                   las = TRUE,
                   xlab = "", ylab = "")
    # Add helpful lines
    if (FALSE) {
      if (metric == "mb") {
        lines(c(0.5, 9.5), c(0, 0), lty = 2, lwd = 0.75, col = "darkgrey")
      }
      if (metric %in% c("R", "d")) {
        lines(c(0.5, 9.5), c(1, 1), lty = 2, lwd = 0.75, col = "darkgrey")
      }
    }

  }) |> invisible()
}) |> invisible()
dev.off()


#########################
#########################
#### Trends

# (optional) TO DO
# * use n_receivers or detection coverage here

#### (optional) Subset algorithms for improved clarity on plot
unique(skills$alg)
skills_for_trends <-
  skills_all |>
  filter(alg %in% c("Null", "COA(30)", "RSP(1)", "ACPF(S)", "ACDCPF(S)")) |>
  as.data.table()
combs <- unique(skills_for_trends$combination)
nc    <- length(unique(skills_for_trends$combination))
table(skills_for_trends$arrangement, skills_for_trends$n_receiver)

#### Build figure
png(here_fig("performance", "relationships.png"),
    height = 4 * nc, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(nc * 2, length(metrics)), oma = c(1, 3, 1, 1), mar = rep(2, 4))
lapply(combs, function(combination) {

  # Identify skill datasets
  # combination <- 1
  skill <- skills_for_trends[system_type == combination, ]

  # Loop over array arrangements
  lapply(split(skill, skill$arrangement), function(sk) {
    # Visualise skill ~ n_receivers
    lapply(metrics, function(metric) {
      # Create plot, distinguishing between algorithms
      # metric <- metrics[1]
      x <- sk$n_receiver + rnorm(nrow(sk), 0, 1)
      y <- sk[[metric]]
      pretty_plot(x, y,
                  ylim = metrics_lim[[metric]],
                  col = scales::alpha(sk$col, 0.25),
                  xlab = "", ylab = "",
                  cex = 0.5, lwd = 0.5)
      # Add trends for each algorithm
      lapply(split(sk, sk$alg, drop = TRUE), function(d) {
       # d <- split(sk, sk$alg)[[1]]
        add_smoother(d$n_receiver, d[[metric]], col = d$col, lwd = 2)
      }) |> invisible()
    }) |> invisible()

  }) |> invisible()
}) |> invisible()
dev.off()


#### Algorithm colour legend
png(here_fig("colour-scheme-algs.png"),
    height = 6, width = 3, units = "in", res = 600)
# (optional) Define subset of relevant algorithms
# algs_legend <- algs |> filter(alg %in% skills_for_trends$alg)
algs_legend <- algs
# Define legend labels as a list of expressions (NULL^1, COA^2) etc., as in maps
algs_legend_label <- list()
for (i in seq_len(nrow(algs_legend))) {
  algs_legend_label[[i]] <- bquote(.(algs_legend$alg_name[i])^.(algs_legend$label[i]))
}
# Make plot
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("top",
       pch = 21,
       lty = 1,
       legend = algs_legend_label,
       col = algs_legend$col,
       pt.bg = scales::alpha(algs_legend$col, 0.25),
       cex = 1.5,
       y.intersp = 1.5
       )
dev.off()


#### End of code.
#########################
#########################
