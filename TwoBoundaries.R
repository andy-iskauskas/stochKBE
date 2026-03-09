###############################
# Emulation of SEIR Model     #
###############################

source("modelFunctions.R")
source("plotting.R")
library(lhs)
library(tidyr)
set.seed(1)

## SEIR Model specification
N_SEIR <- list(
  M = c(S = 900, E = 10, I = 90, R = 0),
  Pre = matrix(
    c(
      0,0,0,0, # 0 -> S
      1,0,0,0, # S -> 0
      1,0,0,0, # S -> E
      0,1,0,0, # E -> 0
      0,1,0,0, # E -> I
      0,0,1,0, # I -> 0 (natural)
      0,0,1,0, # I -> R
      0,0,0,1, # R -> 0
      0,0,0,1  # R -> S
    ),
    ncol = 4, byrow = TRUE
  ),
  Post = matrix(
    c(
      1,0,0,0,
      0,0,0,0,
      0,1,0,0,
      0,0,0,0,
      0,0,1,0,
      0,0,0,0,
      0,0,0,1,
      0,0,0,0,
      1,0,0,0
    ),
    ncol = 4, byrow = TRUE
  ),
  h = function(x, t, th = rep(1, 7)) {
    tot <- sum(x)
    return(
      c(
        th[1]*tot, # Birth
        th[2]*x[1], # S Death
        th[3]*x[1]*x[3]/tot, # Infection
        th[2]*x[2], # E Death
        th[4]*x[2], # Progression
        (th[2]+th[5])*x[3], # I Death
        th[6]*x[3], # Recovery
        th[2]*x[4], # R Death
        th[7]*x[4] # Waning Immunity
      )
    )
  }
)

## Generating bulk results: 120 points, 10 reps each
reps = 10
## Focusing on number of infected at time t=15
out_name = "I"
out_index = 3
t_point = 15
## Ranges for the parameters
ranges <- list(
  lambda = c(1e-5, 1e-4),
  mu = c(1e-5, 1e-4),
  beta = c(0, 1.5),
  epsilon = c(0, 0.21),
  alpha = c(0.01, 0.025),
  gamma = c(0.05, 0.08),
  omega = c(0.002, 0.004)
)
npts <- 140

training_points <- data.frame(t(apply(
  lhs::randomLHS(npts, length(ranges)),
  1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }
))) |> setNames(names(ranges))

## Sample Plot
test_runs <- map(seq_len(nrow(training_points)), function(i) {
  get_results(as.numeric(training_points[i,]), N_SEIR, 10, outs = c(out_name), times = 30, raw = TRUE)
})

run_sample <- sample(length(test_runs), 10)
## Susceptible
plot(1:10, 1:10, xlim = c(0, 30), ylim = c(0, 900), type = 'n',
     main = "SEIRS Model: Susceptible",
     xlab = "Time", ylab = "Number of Susceptible")
for (i in seq_along(run_sample)) {
  these_runs <- test_runs[[run_sample[i]]]
  for (j in seq_len(dim(these_runs)[3])) {
    lines(x = 0:30, y = these_runs[,,j][,1], col = i)
  }
}
abline(v = 15, lty = 2, col = "black")
## Exposed
plot(1:10, 1:10, xlim = c(0, 30), ylim = c(0, 1000), type = 'n',
     main = "SEIRS Model: Exposed",
     xlab = "Time", ylab = "Number of Exposed")
for (i in seq_along(run_sample)) {
  these_runs <- test_runs[[run_sample[i]]]
  for (j in seq_len(dim(these_runs)[3])) {
    lines(x = 0:30, y = these_runs[,,j][,2], col = i)
  }
}
abline(v = 15, lty = 2, col = "black")
## Infected
plot(1:10, 1:10, xlim = c(0, 30), ylim = c(0, 600), type = 'n',
     main = "SEIRS Model: Infected",
     xlab = "Time", ylab = "Number of Infected")
for (i in seq_along(run_sample)) {
  these_runs <- test_runs[[run_sample[i]]]
  for (j in seq_len(dim(these_runs)[3])) {
    lines(x = 0:30, y = these_runs[,,j][,3], col = i)
  }
}
abline(v = 15, lty = 2, col = "black")
## Recovered
plot(1:10, 1:10, xlim = c(0, 30), ylim = c(0, 1000), type = 'n',
     main = "SEIRS Model: Recovered",
     xlab = "Time", ylab = "Number of Recovered")
for (i in seq_along(run_sample)) {
  these_runs <- test_runs[[run_sample[i]]]
  for (j in seq_len(dim(these_runs)[3])) {
    lines(x = 0:30, y = these_runs[,,j][,4], col = i)
  }
}
abline(v = 15, lty = 2, col = "black")

wave0_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(training_points)), function(i) {
  get_results(unlist(training_points[i,], use.names = FALSE), N_SEIR, nreps = reps, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave0_output <- data.frame(wave0_results |> dplyr::group_by(across(all_of(names(ranges)))) |>
                             dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))

## Make the emulators
ems_wave1 <- create_boundary_ems(wave0_results, out_name, ranges, reps,
                                 SEIR_functions, bb_data, N_SEIR, 0, c(3,4), out_index, t_point,
                                 thetas = NULL)

## Create a (comparatively) large grid to evaluate on
big_grid <- expand.grid(
  beta = seq(ranges$beta[1], ranges$beta[2], length.out = 40),
  epsilon = seq(ranges$epsilon[1], ranges$epsilon[2], length.out = 40)
)
for (nm in names(ranges)) {
  if (nm != "beta" && nm != "epsilon")
    big_grid[,nm] <- training_points[1, nm]
}
big_grid <- big_grid[,names(ranges)]

## Expectation and Variance comparison
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
exp_df$Vbulk[exp_df$Vbulk < 0] <- 1e-6
exp_df$Vbound[exp_df$Vbound < 0] <- 1e-6
exp_df$Vboth[exp_df$Vboth < 0] <- 1e-6
exp_df$Vno[exp_df$Vno < 0] <- 1e-6
exp_df_reshape <- tidyr::pivot_longer(exp_df, cols = !c(1:7))
grid_plot(exp_df_reshape, "E", c("beta", "epsilon"), viridoption = "D",
          breaks = c(-100, 0, 10, 25, 50, 100, 200, 300, 400,  500, 1000)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))
grid_plot(exp_df_reshape, "V", c("beta", "epsilon"),
          breaks = c(0, 10, 50, 100, 200, 250, 260, 270,
                     280, 290, 300, 400, 500,
                     1000, 2000, 3000, 4000, 4500),
          viridoption = "C") +
  theme_minimal() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))


var_df <- cbind.data.frame(
  big_grid,
  data.frame(
    Ebulk = ems_wave1$no_boundary$variance[[out_name]]$get_exp(big_grid),
    Vbulk = ems_wave1$no_boundary$variance[[out_name]]$get_cov(big_grid),
    Ebound = ems_wave1$boundary$variance$get_exp(big_grid),
    Vbound = ems_wave1$boundary$variance$get_cov(big_grid),
    Eboth = ems_wave1$boundary_bulk$variance$get_exp(big_grid),
    Vboth = ems_wave1$boundary_bulk$variance$get_cov(big_grid),
    Eno = ems_wave1$no_boundary$variance[[out_name]]$o_em$get_exp(big_grid),
    Vno = ems_wave1$no_boundary$variance[[out_name]]$o_em$get_cov(big_grid, check_neg = FALSE)
  )
)
var_df_reshape <- tidyr::pivot_longer(var_df, cols = !c(1:7))
grid_plot(var_df_reshape, "E", c("beta", "epsilon"))
grid_plot(var_df_reshape, "V", c("beta", "epsilon"))

## Create a collection of points on which to evaluate emulator variance
test_lhs <- lhs::randomLHS(2000, length(ranges))
small_test_grid <- data.frame(t(apply(test_lhs, 1, function(x) {
  x * purrr:::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))

## Create design for next wave of emulation
## Uses parallelisation if available:
# futures: furrr::map (for rep_scores)
# optimParallel: optim (for point_scores)
this_var_em <- ems_wave1$boundary_bulk$variance
## This takes a while, even with parallelisation - included RData file, but
# uncomment to reproduce.
## Set up the cluster for optimParallel; load in required pieces
cl <- makeCluster(8); setDefaultCluster(cl = cl)
clusterEvalQ(cl, library("dplyr"))
clusterEvalQ(cl, library("hmer"))
clusterExport(cl, c("new_point_score", "mean_em_var", "r1", "R1",
                    "part_inv"))
## Create the plan for furrr
plan(multisession, workers = 8)
## Produce new design
new_design <- design_subselect(training_points,
                               ems_wave1$no_boundary$expectation$I$o_em, this_var_em,
                               rep(10, nrow(training_points)), ranges, small_test_grid,
                               rep_max = 2800, pt_max = 280,
                               store_order = TRUE, verbose = TRUE, return_scores = TRUE, ntoadd = 5,
                               boundary_col = c(3,4), boundary_val = 0
                               )
#load("TwoBoundPoints.RData")

## Plotting the results of the proposal
is_old <- rep(c(TRUE, FALSE), each = 140)
is_distant <- purrr::map_lgl(seq_len(nrow(new_design$points)), function(i) {
  sqrt(sum((new_design$points[i,-c(3,4,8)]-purrr::map_dbl(ranges[-c(3:4)], ~sum(.)*0.6))^2)) > 0.01
})
ggplot(data = subset(exp_df_reshape, name == "Vboth"), aes(x = beta, y = eps)) +
  geom_contour_filled(aes(z = value), alpha = 0.6) +
  scale_fill_viridis(discrete = TRUE, name = "Var",
                     labels = c("(0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]")) +
  geom_point(data = new_design$points, col = ifelse(is_old, "grey", "black"),
             pch = ifelse(is_old, 4, 16), size = ifelse(is_distant, 0.75, 1.5)) +
  geom_text(data = new_design$points, aes(y = eps + 0.0025,label = reps),
            col = ifelse(is_old, "grey", "black"), size = ifelse(is_distant, 2, 3))

## Training new emulators
wave1_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(new_design$points)), function(ind) {
  if (ind <= 140 && new_design$points[ind,8] == 10) return(NULL)
  get_results(unlist(new_design$points[ind, 1:7], use.names = FALSE), N_SEIR, nreps = ifelse(ind <= 140, new_design$points[ind,8]-10, new_design$points[ind,8]),
              outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave1_all <- rbind.data.frame(wave0_results, wave1_results)
new_ems <- create_boundary_ems(wave1_all, out_name, ranges, new_design$points$reps,
                               SEIR_functions, bb_data, N_SEIR, 0, c(3,4), out_index, t_point)

## Setup for the basic design
basic_lhs <- lhs::augmentLHS(
  t(apply(training_points, 1, function(x) {
    (x - purrr::map_dbl(ranges, ~.[[1]]))/purrr::map_dbl(ranges, diff)
  })), 140
)
basic_scaled <- data.frame(t(apply(basic_lhs, 1, function(x) {
  x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))
basic_added <- do.call("rbind.data.frame", purrr::map(141:280, function(i) {
  get_results(unlist(basic_scaled[i,], use.names = FALSE), N_SEIR, nreps = 10, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
all_points_basic <- rbind.data.frame(wave0_results, basic_added)
basic_ems <- create_boundary_ems(all_points_basic, out_name, ranges, rep(10, 280),
                                 SEIR_functions, bb_data, N_SEIR, 0, c(3,4), out_index, t_point)
# Setup for the new design with uniform reps
unif_rep_res <- do.call("rbind.data.frame", purrr::map(seq_len(nrow(new_design$points[141:280,])), function(ind) {
  get_results(unlist(new_design$points[ind+140,1:7], use.names = FALSE), N_SEIR, nreps = 10, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
unif_rep_all <- rbind.data.frame(wave0_results, unif_rep_res)
unif_rep_ems <- create_boundary_ems(unif_rep_all, out_name, ranges, rep(10, 280),
                                    SEIR_functions, bb_data, N_SEIR, 0, c(3,4), out_index, t_point)

## Comparing predictive variance across the space
new_vars <- new_ems$boundary_bulk$expectation$get_cov(big_grid)
new_vars[new_vars < 0] <- 1e-6
basic_vars <- basic_ems$boundary_bulk$expectation$get_cov(big_grid)
basic_vars[basic_vars < 0] <- 1e-6
unif_vars <- unif_rep_ems$boundary_bulk$expectation$get_cov(big_grid)
unif_vars[unif_vars < 0] <- 1e-6
old_vars <- exp_df$Vboth
old_vars[old_vars < 0] <- 1e-6

all_var_df <- cbind.data.frame(big_grid, new_vars, old_vars, unif_vars, basic_vars) |>
  setNames(c(names(ranges), "New", "Old", "Uniform", "Naive"))
comparison_plot(all_var_df, c("Old", "Naive", "Uniform", "New"), c("beta", "eps"), "Var")

### Checking against a large set of simulator runs
new_exps <- new_ems$boundary_bulk$expectation$get_exp(big_grid)
basic_exps <- basic_ems$boundary_bulk$expectation$get_exp(big_grid)
unif_exps <- unif_rep_ems$boundary_bulk$expectation$get_exp(big_grid)
old_exps <- exp_df$Eboth
## Including the high-repetition runs for comparison
actual_exps <- test_output[order(test_output$eps, test_output$beta),'exp']

all_exp_df <- cbind.data.frame(big_grid, new_exps, old_exps, unif_exps, basic_exps) |>
  setNames(c(names(ranges), "New", "Old", "Uniform", "Naive"))
comparison_plot(all_exp_df, c("Old", "Naive", "Uniform", "New"), c("beta", "eps"), "Expectation")

## Creating a diagnostic data.frame to plot
diff_df <- all_var_df
for (nm in c("New", "Old", "Naive", "Uniform")) {
  diff_df[,nm] <- abs(all_exp_df[,nm] - actual_exps)/sqrt(all_var_df[,nm])
}
for (i in seq_len(nrow(diff_df))) {
  if (diff_df[i, "beta"] == 0 || diff_df[i, "eps"] == 0) {
    diff_df[i, c("New", "Old", "Naive", "Uniform")] <- c(0,0,0,0)
  }
}
comparison_plot(diff_df, c("Old", "Naive", "Uniform", "New"), c("beta", "eps"), "Error")
