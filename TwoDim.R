###############################
# Emulation of SIR Model - 2d #
###############################
source("modelFunctions.R")
source("baseFunctions.R")
source("plotting.R")
library(lhs)
library(tidyr)
set.seed(2)

## Model set-up for Gillespie algorithm
Num <- 1000
N <- list()
N$M <- c(S = 750, I = 250, R = 0)

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

## Generating bulk results: 20 points, 10 reps each
reps = 10
## Focusing on number of infected at time t=15
out_name = "I"
out_index = 2
t_point = 15
## Ranges for the three parameters
ranges <- list(beta = c(0, 1.5), gamma = c(0, 0.5))
training_points <- data.frame(t(apply(
  lhs::optimumLHS(10*length(ranges), length(ranges)),
  1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }
))) |> setNames(names(ranges))
## Sample Plot
test_runs <- map(seq_len(nrow(training_points)), function(i) {
  get_results(as.numeric(training_points[i,]), N, 10, outs = c(out_name), times = 30, raw = TRUE)
})

run_sample <- sample(length(test_runs), 8)
## Infected
plot(1:10, 1:10, xlim = c(0, 30), ylim = c(0, 1000), type = 'n',
     main = "SIRS Model: Infected",
     xlab = "Time", ylab = "Number of Infected")
for (i in seq_along(run_sample)) {
  these_runs <- test_runs[[run_sample[i]]]
  for (j in seq_len(dim(these_runs)[3])) {
    lines(x = 0:30, y = these_runs[,,j][,2], col = i)
  }
}
abline(v = 15, lty = 2, col = "black")
## Recovered
plot(1:10, 1:10, xlim = c(0, 30), ylim = c(0, 1000), type = 'n',
     main = "SIRS Model: Recovered",
     xlab = "Time", ylab = "Number of Recovered")
for (i in seq_along(run_sample)) {
  these_runs <- test_runs[[run_sample[i]]]
  for (j in seq_len(dim(these_runs)[3])) {
    lines(x = 0:30, y = these_runs[,,j][,3], col = i)
  }
}
abline(v = 15, lty = 2, col = "black")
## Susceptible
plot(1:10, 1:10, xlim = c(0, 30), ylim = c(0, 1000), type = 'n',
     main = "SIRS Model: Susceptible",
     xlab = "Time", ylab = "Number of Susceptible")
for (i in seq_along(run_sample)) {
  these_runs <- test_runs[[run_sample[i]]]
  for (j in seq_len(dim(these_runs)[3])) {
    lines(x = 0:30, y = these_runs[,,j][,1], col = i)
  }
}
abline(v = 15, lty = 2, col = "black")

wave0_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(training_points)), function(i) {
  get_results(unlist(training_points[i,], use.names = FALSE), N, nreps = reps, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave0_output <- data.frame(wave0_results |> dplyr::group_by(across(all_of(names(ranges)))) |>
                             dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))

## Make the emulators
ems_wave1 <- create_boundary_ems(wave0_results, out_name, ranges, reps,
                                 SIR_functions, bb_data, N, c(0), c(1), out_index, t_point)
this_var_em <- ems_wave1$boundary_bulk$variance

## Create a 20x20 grid of points to evaluate mean emulator variance on
small_test_grid <- expand.grid(beta = seq(0, 1.5, length.out = 10), gamma = seq(0, 0.5, length.out = 10))

source("baseFunctions.R")
cl <- makeCluster(8); setDefaultCluster(cl = cl)
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(hmer))
clusterExport(cl, c("imspe", "part_inv", "small_test_grid", "R1", "r1"))
new_design_test <- point_design(training_points, ems_wave1$boundary$expectation, ems_wave1$boundary_bulk$variance,
                                rep(10, 20), ranges, small_test_grid, 40, verbose = TRUE, return_scores = TRUE,
                                in_par = TRUE)
plot(x = new_design_test$points$beta, y = new_design_test$points$gamma, pch = 16, col = rep(c("grey", "black"), each = 20))

design_with_reps <- rep_allocate(new_design_test$points, ems_wave1$boundary_bulk$variance, 400, 2)

## Create design for next wave of emulation
## Next 4 lines relevant if optimParallel is installed
cl <- makeCluster(8); setDefaultCluster(cl = cl)
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(hmer))
clusterExport(cl, c("new_point_score", "mean_em_var", "r1", "R1", "part_inv"))
## Next line relevant if future is installed
plan(multisession, workers = 8)
new_design <- design_subselect(training_points,
                          ems_wave1$no_boundary$expectation$I$o_em, this_var_em,
                          rep(10, 20), ranges, small_test_grid, rep_max = 400, pt_max = 40,
                          store_order = TRUE, verbose = TRUE, return_scores = TRUE, ntoadd = 2)

## Plotting results
plot(x = seq_along(new_design$rep_scores), y = new_design$rep_scores, type = 'l',
     main = "Predictive Variance", xlab = "Proposal", ylab = "Variance",
     ylim = c(0, max(new_design$rep_scores)))
lines(x = seq_along(new_design$pt_scores[!is.infinite(new_design$pt_scores)]), y = new_design$pt_scores[!is.infinite(new_design$pt_scores)], col = 'blue')
legend("topright", inset = 0.05, legend = c("New Point", "Extra Rep"),
       lty = 1, col = c("blue", "black"))

bgn <- 50
big_grid <- expand.grid(
  beta = seq(ranges$beta[[1]], ranges$beta[[2]], length.out = bgn),
  gamma = seq(ranges$gamma[[1]], ranges$gamma[[2]], length.out = bgn)
)
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

## Plotting the result of the proposal
get_loc <- function(y) {
  if (y > 0.48) return(y-0.012)
  else return(y+0.012)
}
ggplot(data = subset(exp_df_reshape, name == "Vboth"), aes(x = beta, y = gamma)) +
  geom_raster(aes(fill = value), interpolate = TRUE) +
  scale_fill_gradientn(name = "Var", colours = viridis::viridis(17, option = "A"),
                       values = c(0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
  geom_point(data = new_design$points, col = rep(c("grey40", "white"), each = 20), size = rep(c(0.8, 1), each = 20)) +
  geom_text(data = new_design$points, aes(label = reps),
            y = sapply(new_design$points$gamma, get_loc),
            col = rep(c("grey40", "white"), each = 20), size = rep(c(3, 4), each = 20)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0.01)) +
  scale_y_continuous(expand = c(0.01,0)) +
  labs(x = TeX("$\\beta$"), y = TeX("$\\gamma$"))
  

### Training new emulator and comparing to other proposal methods
# Three methods of point proposal are considered: one using the improved design with
# the corresponding suggested repetitions at each point (new_ems); one using a standard
# Latin hypercube design with uniformly many repetitions at each point (basic_ems), and
# one using the new design points but allocating equal numbers of repetitions to all
# points (unif_rep_ems).
old_vars <- ems_wave1$boundary_bulk$expectation$get_cov(big_grid)
old_vars[old_vars < 0] <- 1e-6

# Setup for the new design
wave1_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(new_design$points)), function(i) {
  if (i <= 20 && new_design$points[i,3] == 10) return(NULL)
  get_results(unlist(new_design$points[i, 1:2], use.names = FALSE), N, nreps = ifelse(i <= 20, new_design$points[i,3]-10, new_design$points[i,3]),
              outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave1_all <- rbind.data.frame(wave0_results, wave1_results)
wave1_output <- data.frame(wave1_all |> dplyr::group_by(across(all_of(names(ranges)))) |>
                             dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]]), reps = length(.data[[out_name]])))
new_ems <- create_boundary_ems(wave1_all, out_name, ranges, new_design$points$reps,
                               SIR_functions, bb_data, N, 0, c(1), out_index, t_point)

# Setup for the basic design
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

# Setup for the new design with uniform reps
unif_rep_res <- do.call("rbind.data.frame", purrr::map(seq_len(nrow(new_design$points)), function(i) {
  get_results(unlist(new_design$points[i,1:2], use.names = FALSE), N, nreps = 10, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
unif_rep_ems <- create_boundary_ems(unif_rep_res, out_name, ranges, rep(10, 40),
                                    SIR_functions, bb_data, N, 0, c(1), out_index, t_point)


## Comparing predictive variance across the space
new_vars <- new_ems$boundary_bulk$expectation$get_cov(big_grid)
new_vars[new_vars < 0] <- 1e-6
basic_vars <- basic_ems$boundary_bulk$expectation$get_cov(big_grid)
basic_vars[basic_vars < 0] <- 1e-6
unif_rep_vars <- unif_rep_ems$boundary_bulk$expectation$get_cov(big_grid)
unif_rep_vars[unif_rep_vars < 0] <- 1e-6

## Plotting the results: used here are the three emulator sets trained above with
# the different proposals, as well as the original wave 1 stoch KBE emulators
all_var_df <- cbind.data.frame(cbind.data.frame(big_grid, new_vars), cbind.data.frame(basic_vars, unif_rep_vars, old_vars)) |>
  setNames(c('beta', 'gamma', 'New', 'Naive', "Uniform", "Old"))

comparison_plot(all_var_df, c("Old", "Naive", "Uniform", "New"), c("beta", "gamma"), "Variance",
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
#### Paper Plots End Here ####


############################################
## Lasciate ogne speranza, voi ch'intrate ##
############################################
## Optional - animating the point proposal process
library(gganimate)
library(gifski)
library(scales)
# Function to convert proposed points into ordered data.frame of proposals over 'time'
create_transition_data_frame <- function(data, order, nadd = 1, original_points = NULL) {
  order <- rep(order, each = nadd)
  df_list <- purrr::map(seq_along(order), function(i) {
    relev_order <- order[1:i]
    order_count <- purrr::map_dbl(unique(relev_order), ~sum(relev_order == .))
    df <- data[unique(relev_order),]
    df$reps <- order_count
    df$time <- i
    df
  })
  return(do.call('rbind.data.frame', df_list))
}
trans_df <- create_transition_data_frame(new_design$points, new_design$order, nadd = 2)
for (i in seq_len(nrow(training_points))) {
  for (j in unique(trans_df$time)) {
    has_entry <- which(trans_df$time == j & trans_df$beta == training_points[i,"beta"] & trans_df$gamma == training_points[i,"gamma"])
    if (length(has_entry) == 0) {
      added_entry <- cbind.data.frame(training_points[i,1:2], data.frame(reps = 10, time = j)) |>
        setNames(c(names(training_points)[1:2], "reps", "time"))
      trans_df <- rbind.data.frame(trans_df, added_entry)
    }
    else {
      trans_df[has_entry,"reps"] <- trans_df[has_entry, "reps"] + 10
    }
  }
}
trans_df <- trans_df[order(trans_df$time),]
colour_mod <- c()
for (i in seq_len(nrow(trans_df))) {
  pt <- trans_df[i,1:2]
  in_train <- purrr::map_lgl(seq_len(nrow(training_points)), function(j) {
    all(pt == training_points[j, 1:2])
  })
  if (any(in_train)) colour_mod <- c(colour_mod, 1)
  else colour_mod <- c(colour_mod, 0)
}
trans_df$mod <- colour_mod
trans_df$mod <- factor(trans_df$mod, levels = c(0, 1))
exp_df_rename <- setNames(exp_df_reshape, c('b', 'g', 'n', 'v'))
anim_alt <- ggplot(data = trans_df, aes(x = beta, y = gamma)) +
  geom_raster(data = subset(exp_df_rename, n == "Vboth"),
              aes(x = b, y = g, fill = v), interpolate = TRUE) +
  scale_fill_viridis(name = "Var") +
  geom_point(data = trans_df, aes(x = beta, y = gamma, colour = mod, size = 2-(as.numeric(mod)-1)*0.5)) +
  geom_text(data = trans_df, aes(y = gamma + 0.01, label = reps, size = 5-2*(as.numeric(mod)-1), colour = mod)) +
  scale_colour_manual(values = c("1" = "grey80", "0" = "white")) +
  transition_time(time) +
  ease_aes("linear") +
  theme(legend.position = "none") +
  ggtitle("Repetitions Placed: {frame_time}")
animate(anim_alt, nframes = 200, end_pause = 25, height = 800, width = 800)
