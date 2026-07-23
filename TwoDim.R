##################################
# 2d SIR model with one boundary #
##################################

## Source functions and load libraries
source("modelFunctions.R")
source("baseFunctions.R")
source("plotting.R")
require(lhs)
require(tidyr)
# Set seed for reproducibility
set.seed(1)

## Gillespie algorithm set-up
# The model is an SIRS model with fixed waning-immunity rate, making it
# functionally a 2-dimensional problem (but with interesting dynamics). The
# boundary is at beta = 0.
N <- list()
N$M <- c(S = 750, I = 250, R = 0)
Num <- sum(N$M)
N$Pre <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, byrow = TRUE)
N$Post <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
N$h <- function(x, t, th = rep(1, 2)) {
  Num = 1000
  return(
    c(th[1]*x[1]*x[2]/Num,
      th[2]*x[2],
      0.24*x[3])
  )
}

## Generate results: 20 points, each with 10 realisations
# Focus on the I compartment (output 2) at time t=15
reps <- 10
out_name <- "I"
out_index <- 2
t_point <- 15
# Define ranges and generate a LHD over the points
ranges <- list(beta = c(0, 1.5), gamma = c(0, 0.5))
training_points <- data.frame(t(apply(
  lhs::optimumLHS(10*length(ranges), length(ranges)),
  1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }
))) |> setNames(names(ranges))

## Get results
wave0_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(training_points)), function(i) {
  get_results(unlist(training_points[i,], use.names = FALSE), N, nreps = reps, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave0_output <- data.frame(wave0_results |> dplyr::group_by(across(all_of(names(ranges)))) |>
                             dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))

## Construct the emulators using `create_boundary_ems()`
# Explicitly create a variable for the base emulator and variance emulator for later
ems_wave1 <- create_boundary_ems(wave0_results, out_name, ranges, reps,
                                 SIR_functions, bb_data, N, c(0), c(1), out_index, t_point,
                                 thetas = NULL)
this_base_em <- ems_wave1$no_boundary$expectation$I$o_em
this_var_em <- ems_wave1$boundary_bulk$variance

## Create a 'big' grid to evaluate results on
bgn <- 50
big_grid <- expand.grid(
  beta = seq(ranges$beta[[1]], ranges$beta[[2]], length.out = bgn),
  gamma = seq(ranges$gamma[[1]], ranges$gamma[[2]], length.out = bgn)
)
# For the trained emulators, want to compare the results from the base (untrained)
# emulators, the emulators that only have access to boundary information, those 
# with only 'bulk' knowledge, and those which have both boundary and bulk info.
# Create data.frames with all of this information for both the expectation and
# variance emulators, reshaping for the purpose of plotting.
exp_df <- cbind.data.frame(
  big_grid,
  data.frame(
    Ebulk = ems_wave1$no_boundary$expectation[[out_name]]$get_exp(big_grid),
    Vbulk = ems_wave1$no_boundary$expectation[[out_name]]$get_cov(big_grid),
    Ebound = ems_wave1$boundary$expectation$get_exp(big_grid),
    Vbound = ems_wave1$boundary$expectation$get_cov(big_grid),
    Eboth = ems_wave1$boundary_bulk$expectation$get_exp(big_grid),
    Vboth = ems_wave1$boundary_bulk$expectation$get_cov(big_grid),
    Eno = ems_wave1$no_boundary$expectation[[out_name]]$o_em$get_exp(big_grid),
    Vno = ems_wave1$no_boundary$expectation[[out_name]]$o_em$get_cov(big_grid)
  )
)
exp_df_reshape <- tidyr::pivot_longer(exp_df, cols = !c(beta, gamma))
var_df <- cbind.data.frame(
  big_grid,
  data.frame(
    Ebulk = ems_wave1$no_boundary$variance[[out_name]]$get_exp(big_grid, check_neg = FALSE),
    Vbulk = ems_wave1$no_boundary$variance[[out_name]]$get_cov(big_grid),
    Ebound = ems_wave1$boundary$variance$get_exp(big_grid),
    Vbound = ems_wave1$boundary$variance$get_cov(big_grid),
    Eboth = ems_wave1$boundary_bulk$variance$get_exp(big_grid),
    Vboth = ems_wave1$boundary_bulk$variance$get_cov(big_grid),
    Eno = ems_wave1$no_boundary$variance[[out_name]]$o_em$get_exp(big_grid),
    Vno = ems_wave1$no_boundary$variance[[out_name]]$o_em$get_cov(big_grid)
  )
)
var_df_reshape <- tidyr::pivot_longer(var_df, cols = !c(1:2))
## Plot the results. Colour scale breaks have been chosen by-hand.
## Order of plotting is:
# Mean emulator: expectation then variance
# Variance emulator: expectation then variance
grid_plot(exp_df_reshape, "E", c("beta", "gamma"), "Mean", wave0_output, viridoption = "D",
          breaks = c(-250, 0, 100, 200, 300, 500, 750, 1000, 1200),
          labels = c("(-250, 0]", "[0, 100)", "[100, 200)", "[200, 300)",
                     "[300, 500)", "[500, 750)", "[750, 1000)", "[1000, 1200)")) +
  theme_minimal() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))
grid_plot(exp_df_reshape, "V", c("beta", "gamma"), "Mean", wave0_output, viridoption = "C",
          breaks = c(-100, 10, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000, 1e6),
          labels = c(TeX("$\\[0, 10)$"), TeX("$\\[10, 50)$"),
                     TeX("$\\[50, 100)$"), TeX("$\\[100, 200)$"),
                     TeX("$\\[200, 500)$"), TeX("$\\[500, 10^3)$"),
                     TeX(r"($\[10^3, 2\times 10^3)$)"), TeX(r"($\[2\times 10^3, 5\times 10^3)$)"),
                     TeX(r"($\[5\times 10^3, 10^4)$)"), TeX(r"($\[10^4, 5\times 10^4)$)"),
                     TeX(r"($\[5\times 10^4, 10^6)$)"))) +
  theme_minimal() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))
grid_plot(var_df_reshape, "E", c("beta", "gamma"), "Variance", wave0_output,
          viridoption = "D",
          breaks = c(-200, 0, 5, 10, 50, 100, 500, 1000, 2000, 5000)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))
grid_plot(var_df_reshape, "V", c("beta", "gamma"), "Variance", wave0_output,
          viridoption = "C",
          breaks = c(-1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8),
          labels = c(TeX("$\\[10^0, 10^1)$"),
                     TeX("$\\[10^1, 10^2)$"),
                     TeX("$\\[10^2, 10^3)$"),
                     TeX("$\\[10^3, 10^4)$"),
                     TeX("$\\[10^4, 10^5)$"),
                     TeX("$\\[10^5, 10^6)$"),
                     TeX("$\\[10^6, 10^7)$"),
                     TeX("$\\[10^7, 10^8)$"))) +
  theme_minimal() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))

# New Proposal
## Set up the parallelisation, and propose a new design of 20 points with another
# 200 reps to be shared around.
cl <- makeCluster(8); setDefaultCluster(cl = cl)
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(hmer))
clusterExport(cl, c("imspe", "part_inv", "R1", "r1"))
## We assume two repetitions will be placed per design point initially...
new_design <- point_design(training_points, this_base_em, this_var_em,
                           rep(10, 20), ranges, 40^2, 40, verbose = TRUE, return_scores = TRUE,
                           in_par = TRUE, boundary_col = 1, boundary_val = 0, nrepsadd = 2)
## ...then we allocate the full complement of reps.
new_design_with_reps <- rep_allocate(new_design$points, this_var_em, 400, 2)

## Plotting the proposal result
### A quick-and-dirty function to ensure that annotated repetition numbers don't
# exceed the bounds of the plot
get_loc <- function(y) {
  if (y > 0.48) return(y-0.012)
  else return(y+0.012)
}
## Plot the results: old design points are smaller and in grey. Background is the
# original emulator variance across the space (derived from exp_df above)
ggplot(data = subset(exp_df_reshape, name == "Vboth"), aes(x = beta, y = gamma)) +
  geom_raster(aes(fill = value), interpolate = TRUE) +
  scale_fill_gradientn(name = "Var", colours = viridis::viridis(17, option = "A"),
                       values = c(0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
  geom_point(data = new_design_with_reps, col = rep(c("grey40", "white"), each = 20), size = rep(c(0.8, 1), each = 20)) +
  geom_text(data = new_design_with_reps, aes(label = reps),
            y = sapply(new_design_with_reps$gamma, get_loc),
            col = rep(c("grey40", "white"), each = 20), size = rep(c(3, 4), each = 20)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0.01)) +
  scale_y_continuous(expand = c(0.01,0)) +
  labs(title = "New Design", x = TeX("$\\beta$"), y = TeX("$\\gamma$"))
## We can do the same with the predicted variance emulator expectation (aka stochasticity)
# Useful for looking at where the reps have been placed
ggplot(data = subset(var_df_reshape, name == "Eboth"), aes(x = beta, y = gamma)) +
  geom_raster(aes(fill = value), interpolate = TRUE) +
  scale_fill_gradientn(name = "Exp", colours = viridis::viridis(17, option = "B"),
                       values = c(0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
  geom_point(data = new_design_with_reps, col = rep(c("grey40", "white"), each = 20), size = rep(c(0.8, 1), each = 20)) +
  geom_text(data = new_design_with_reps, aes(label = reps),
            y = sapply(new_design_with_reps$gamma, get_loc),
            col = rep(c("grey40", "white"), each = 20), size = rep(c(3, 4), each = 20)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0.01)) +
  scale_y_continuous(expand = c(0.01,0)) +
  labs(title = "New Design", x = TeX("$\\beta$"), y = TeX("$\\gamma$"))

## Alternate plotting: show the updated emulator variance as each point is added.
# For this, we need variance data.frames for each intermediate design
design_dfs <- purrr::map(21:40, function(i) {
  subset_data <- new_design$points[1:(i-1),1:2]
  subset_data_mutate <- (subset_data |> dplyr::mutate(across(all_of(1), ~0)))
  v_em_vals <- this_var_em$get_exp(subset_data)
  v_em_vals[v_em_vals < 0] <- 1e-6
  rps <- new_design$points$reps
  start_mat <- R1(subset_data, subset_data, this_base_em, 0, 1) *
    this_base_em$get_cov(subset_data_mutate, full = TRUE) +
    diag(c(v_em_vals/rps[1:(i-1)]))
  start_inv <- tryCatch(chol2inv(chol(start_mat)), error = function(e) MASS::ginv(start_mat))
  x_mod <- new_design$points[i,1:2,drop=FALSE]
  outpt <- imspe(big_grid, this_base_em,
                  subset_data, x_mod, start_inv, return.raw = TRUE)
  return(outpt)
})
orig_df <- exp_df[,c("beta", "gamma", "Vboth")] |> setNames(c("beta", "gamma", "V"))
## Regularisation for (very small) negative values
orig_df[orig_df$V < 0, "V"] <- 1e-6
# A not necessarily nice combining of the initial data.frame and the new ones.
all_dfs <- list(orig_df)
for (i in 1:20) {
  this_df <- design_dfs[[i]]
  this_df[this_df$V < 0, "V"] <- 1e-6
  all_dfs[[i+1]] <- this_df
}
## Create the plots: 21 in total
gplots <- purrr::map(21:41, function(i) {
  df <- all_dfs[[i-20]]
  g <- ggplot(data = df, aes(x = beta, y = gamma)) +
    geom_raster(aes(fill = V), interpolate = TRUE) +
    scale_fill_gradientn(name = "Var", colours = viridis::viridis(17, option = "A"),
                         values = c(0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1),
                         limits = c(0, 8000))
  if (i == 21)
    g <- g + geom_point(data = new_design_with_reps[1:(i-1),],
                        col = rep(c('grey40', 'white'), times = c(20, i-21)),
                        size = rep(c(0.8, 1), times = c(20, i-21))) +
    geom_text(data = new_design_with_reps[1:(i-1),], aes(label = reps),
              y = sapply(new_design_with_reps[1:(i-1),"gamma"], get_loc),
              col = rep(c("grey40", "white"), times = c(20, i-21)),
              size = rep(c(3,4), times = c(20, i-21)))
  else
    g <- g + geom_point(data = new_design_with_reps[1:(i-1),],
                        col = rep(c('grey40', 'white', 'red'), times = c(20, i-22, 1)),
                        size = rep(c(0.8, 1), times = c(20, i-21))) +
    geom_text(data = new_design_with_reps[1:(i-1),], aes(label = reps),
              y = sapply(new_design_with_reps[1:(i-1),"gamma"], get_loc),
              col = rep(c("grey40", "white", 'red'), times = c(20, i-22, 1)),
              size = rep(c(3,4), times = c(20, i-21)))
  g <- g + theme_minimal() +
    scale_x_continuous(expand = c(0,0.01)) +
    scale_y_continuous(expand = c(0.01,0)) +
    labs(title = paste("Emulator Variance after", i-21, "new points proposed"),
         x = TeX("$\\beta$"), y = TeX("$\\gamma$"))
})
## Optional: save to a pdf
# pdf(file = "../../StochPropPlot.pdf")
# for (i in gplots) print(i)
# dev.off()

### Comparison of different methods of design
# We compare boundary emulators trained on different design methods:
# 1) A 'naive', augmented LHD, design with equal reps/point
# 2) The design above with equal reps/point
# 3) The full new design with reps chosen according to `rep_allocate()`.
## 1 - The basic design
# Setup
basic_lhs <- lhs::augmentLHS(
  t(apply(training_points, 1, function(x) {
    (x - purrr::map_dbl(ranges, ~.[[1]]))/purrr::map_dbl(ranges, diff)
  })), 20
)
basic_scaled <- data.frame(t(apply(basic_lhs, 1, function(x) {
  x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))
basic_added <- do.call("rbind.data.frame", purrr::map(21:40, function(i) {
  get_results(unlist(basic_scaled[i,], use.names = FALSE), N, nreps = 10, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
all_points_basic <- rbind.data.frame(wave0_results, basic_added)
basic_ems <- create_boundary_ems(all_points_basic, out_name, ranges, rep(10, 40),
                                 SIR_functions, bb_data, N, 0, c(1), out_index, t_point)
## 2 - The new design with uniform reps
# Setup
unif_rep_res <- do.call("rbind.data.frame", purrr::map(seq_len(nrow(new_design$points)), function(i) {
  get_results(unlist(new_design$points[i,1:2], use.names = FALSE), N, nreps = 10, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
unif_rep_ems <- create_boundary_ems(unif_rep_res, out_name, ranges, rep(10, 40),
                                     SIR_functions, bb_data, N, 0, c(1), out_index, t_point)
## 3 - The full design
# Setup
new_design_res <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(new_design$points)), function(i) {
  if (i <= 20 && new_design$points[i,3] == 10) return(NULL)
  get_results(unlist(new_design$points[i, 1:2], use.names = FALSE), N, nreps = ifelse(i <= 20, new_design$points[i,3]-10, new_design$points[i,3]),
              outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
all_res <- rbind.data.frame(wave0_results, new_design_res)
all_summary <- data.frame(all_res |> dplyr::group_by(across(all_of(names(ranges)))) |>
                              dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]]), reps = length(.data[[out_name]])))
new_design_ems <- create_boundary_ems(all_res, out_name, ranges, new_design$points$reps,
                                SIR_functions, bb_data, N, 0, c(1), out_index, t_point)
# Collating results: we include the predictive variance from the first wave of emulators too
em_names <- c("basic_ems", "unif_rep_ems", "new_design_ems")
all_var_df <- cbind.data.frame(cbind.data.frame(big_grid, exp_df$Vboth),
                               do.call('cbind.data.frame', purrr::map(em_names, function(nm) {
  these_ems <- get(nm)
  these_vars <- these_ems$boundary_bulk$expectation$get_cov(big_grid)
  these_vars[these_vars < 0] <- 1e-6
  return(these_vars)
}))) |> setNames(c("beta", "gamma", "Old", "Naive", "Uniform", "New"))
## Checking the mean predictive variance for each set
apply(all_var_df[,3:6], 2, mean)
## Comparative plot
comparison_plot(all_var_df, c("Old", "Uniform", "New", "Naive"), c("beta", "gamma"), "Variance", 
                breaks = c(0, 10, 50, 100, 200, 500, 1000, 10000, 100000),
                labels = c(
                  TeX(r"($\[0, 10)$)"), TeX(r"($\[10, 50)$)"),
                  TeX(r"($\[50, 100)$)"), TeX(r"($\[100, 200)$)"),
                  TeX(r"($\[200, 500)$)"), TeX(r"($\[500, 10^4)$)"),
                  TeX(r"($\[10^4, 10^5)$)"), TeX(r"($\[10^5, 10^6)$)")
                ),
                viridoption = "C") +
  theme_minimal() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))
