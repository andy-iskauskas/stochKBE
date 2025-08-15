###############################
# Emulation of SEIR Model     #
###############################

source("modelFunctions.R")
source("plotting.R")
library(lhs)
library(tidyr)
set.seed(42)

## SEIR Model specification
N_SEIR <- list(
  M = c(S = 940, E = 10, I = 50, R = 0),
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

## Generating bulk results: 20 points, 10 reps each
reps = 10
## Focusing on number of infected at time t=15
out_name = "I"
out_index = 3
t_point = 15
## Ranges for the three parameters
ranges <- list(
  lambda = c(1e-5, 1e-4),
  mu = c(1e-5, 1e-4),
  beta = c(0, 0.5),
  eps = c(0, 0.21),
  alpha = c(0.01, 0.025),
  gamma = c(0.05, 0.08),
  omega = c(0.002, 0.004)
)
set.seed(1)
training_points <- data.frame(t(apply(
  lhs::randomLHS(20*length(ranges), length(ranges)),
  1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }
))) |> setNames(names(ranges))
wave0_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(training_points)), function(i) {
  get_results(unlist(training_points[i,], use.names = FALSE), N_SEIR, nreps = reps, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave0_output <- data.frame(wave0_results |> dplyr::group_by(across(all_of(names(ranges)))) |>
                             dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))

## Make the emulators
ems_wave1 <- create_boundary_ems(wave0_results, out_name, ranges, reps,
                                 SEIR_functions, bb_data, N_SEIR, 0, c(3,4), out_index, t_point)

big_grid <- expand.grid(
  beta = seq(ranges$beta[1], ranges$beta[2], length.out = 20),
  eps = seq(ranges$eps[1], ranges$eps[2], length.out = 20)
)
for (nm in names(ranges)) {
  if (nm != "beta" && nm != "eps")
    big_grid[,nm] <- 0.6*sum(ranges[[nm]])
}
big_grid <- big_grid[,names(ranges)]

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
exp_df_reshape <- tidyr::pivot_longer(exp_df, cols = !c(1:7))
ggplot(data = exp_df_reshape[grepl("E", exp_df_reshape$name),], aes(x = beta, y = eps, z = value)) +
  geom_contour_filled() +
  facet_wrap(vars(name), nrow = 2)
ggplot(data = exp_df_reshape[grepl("V", exp_df_reshape$name),], aes(x = beta, y = eps, z = value)) +
  geom_contour_filled() +
  facet_wrap(vars(name), nrow = 2)

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
ggplot(data = var_df_reshape[grepl("E", var_df_reshape$name),], aes(x = beta, y = eps)) +
  geom_contour_filled(aes(z = value)) +
  facet_wrap(vars(name), nrow = 2)
ggplot(data = var_df_reshape[grepl("V", var_df_reshape$name),], aes(x = beta, y = eps)) +
  geom_contour_filled(aes(z = value)) +
  facet_wrap(vars(name), nrow = 2)


## Create a collection of points on which to evaluate emulator variance
test_lhs <- lhs::randomLHS(1000, length(ranges))
small_test_grid <- data.frame(t(apply(test_lhs, 1, function(x) {
  x * purrr:::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))

## Create design for next wave of emulation
## Uses parallelisation:
# futures: furrr::map (for rep_scores)
# optimParallel: optim (for point_scores)
this_var_em <- ems_wave1$boundary_bulk$variance
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
                               rep_max = 20*nrow(training_points), pt_max = 2*nrow(training_points),
                               store_order = TRUE, verbose = TRUE, return_scores = TRUE, ntoadd = 2,
                               boundary_col = c(3,4), boundary_val = 0
                               )

## Plotting the results
is_old <- rep(c(TRUE, FALSE), each = 140)
is_distant <- purrr::map_lgl(seq_len(nrow(new_design$points)), function(i) {
  sqrt(sum((new_design$points[i,-c(3,4,8)]-purrr::map_dbl(ranges[-c(3:4)], ~sum(.)*0.6))^2)) > 0.01
})
is_distant
ggplot(data = subset(exp_df_reshape, name == "Vboth"), aes(x = beta, y = eps)) +
  geom_contour_filled(aes(z = value)) +
  scale_fill_viridis(discrete = TRUE, name = "Var",
                     labels = c("(0, 5]", "(5, 10]", "(10, 15]", "(15, 20]", "(20, 25]", "(25, 30]")) +
  geom_point(data = new_design$points, col = ifelse(is_old, "grey", "black"),
             pch = ifelse(is_old, 4, 16), size = ifelse(is_distant, 0.75, 1.5)) +
  geom_text(data = new_design$points, aes(y = eps + 0.0025,label = reps),
            col = ifelse(is_old, "grey", "black"), size = ifelse(is_distant, 2, 3))

## Training a new emulator on the proposal
wave1_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(new_design$points)), function(i) {
  get_results(unlist(new_design$points[i, 1:7], use.names = FALSE), N_SEIR, nreps = new_design$points[i,8],
              outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave1_output <- data.frame(wave1_results |> dplyr::group_by(across(all_of(names(ranges)))) |>
                             dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))
new_ems <- create_boundary_ems(wave1_results, out_name, ranges, reps,
                               SEIR_functions, bb_data, N_SEIR, 0, c(3,4), out_index, t_point)

# Setup for the basic design
training_points <- data.frame(t(apply(
  lhs::randomLHS(20*length(ranges), length(ranges)),
  1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }
))) |> setNames(names(ranges))

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
unif_rep_res <- do.call("rbind.data.frame", purrr::map(seq_len(nrow(new_design$points)), function(i) {
  get_results(unlist(new_design$points[i,1:7], use.names = FALSE), N_SEIR, nreps = 10, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
unif_rep_ems <- create_boundary_ems(unif_rep_res, out_name, ranges, rep(10, 280),
                                    SEIR_functions, bb_data, N_SEIR, 0, c(3,4), out_index, t_point)

## Comparing predictive variance across the space
new_vars <- new_ems$boundary_bulk$expectation$get_cov(small_test_grid)
basic_vars <- basic_ems$boundary_bulk$expectation$get_cov(small_test_grid)
unif_rep_vars <- unif_rep_ems$boundary_bulk$expectation$get_cov(small_test_grid)
# Mean predictive variance
mean(new_vars)
mean(basic_vars)
mean(unif_rep_vars)
mean(basic_vars)/mean(new_vars)

relev_new_vars <- new_ems$boundary_bulk$expectation$get_cov(big_grid)
relev_basic_vars <- basic_ems$boundary_bulk$expectation$get_cov(big_grid)
relev_unif_vars <- unif_rep_ems$boundary_bulk$expectation$get_cov(big_grid)

all_var_df <- cbind.data.frame(cbind.data.frame(big_grid, relev_new_vars), cbind.data.frame(exp_df$Vboth, relev_unif_vars, relev_basic_vars)) |>
  setNames(c(names(ranges), 'New', 'Old', "Uniform", "Naive"))
reshape_var_df <- tidyr::pivot_longer(all_var_df, cols = !names(ranges))
reshape_var_df$name <- factor(reshape_var_df$name, levels = c("Old", "Naive", "Uniform", "New"))

ggplot(data = reshape_var_df, aes(x = beta, y = eps, z = value)) +
  geom_contour_filled() +
  scale_fill_viridis(discrete = TRUE, name = "Variance") +
  facet_wrap(vars(name), nrow = 2)
