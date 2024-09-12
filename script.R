##########################
# Emulation of SIR Model #
##########################

source("modelFunctions.R")
source("plotting.R")
library(lhs)
library(tidyr)
set.seed(42)

## Model set-up for Gillespie algorithm
Num <- 1000
N <- list()
N$M <- c(750, 250, 0)

N$Pre <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, byrow = TRUE)
N$Post <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
N$h <- function(x, t, th = rep(1, 3)) {
  Num = 1000
  return(
    c(th[1]*x[1]*x[2]/Num,
      th[2]*x[2],
      th[3]*x[3])
  )
}

## Generating bulk results: 60 points, 10 reps each
reps = 10
## Focusing on number of infected at time t=15
out_name = "I"
out_index = 2
t_point = 15
## Ranges for the three parameters
ranges <- list(beta = c(0, 1.5), gamma = c(0, 0.5), omega = c(0, 0.5))
set.seed(42)
training_points <- data.frame(t(apply(
  lhs::optimumLHS(20*length(ranges), length(ranges)),
  1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }
))) |> setNames(names(ranges))
wave0_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(training_points)), function(i) {
  get_results(unlist(training_points[i,], use.names = FALSE), nreps = reps, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave0_output <- data.frame(wave0_results |> dplyr::group_by(beta, gamma, omega) |>
                            dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))

## Make the emulators
ems_wave1 <- create_boundary_ems(wave0_results, out_name, ranges, reps,
                                 b_exp, b_cov, bb_exp, bb_cov, bb_imp)

## Checking that the boundary uncertainty makes sense
test_points <- data.frame(t(apply(
  lhs::maximinLHS(50, length(ranges)),
  1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }
))) |> setNames(names(ranges))
r_beta <- test_points$beta
## Checking the uncertainty of the stochasticity surface emulator
just_bound <- ems_wave1$boundary$variance$get_cov(test_points)[order(r_beta)]
just_bulk <- ems_wave1$no_boundary$variance[[out_name]]$get_cov(test_points)[order(r_beta)]
bound_bulk <- ems_wave1$boundary_bulk$variance$get_cov(test_points)[order(r_beta)]
beta_for_plot <- r_beta[order(r_beta)]
plot(x = beta_for_plot, y = just_bound,
     ylim = range(c(just_bound, just_bulk, bound_bulk)),
     main = "Variance Emulator Uncertainty", ylab = "Uncertainty",
     xlab = "beta", type = 'l')
lines(x = beta_for_plot, y = just_bulk, col = 'blue')
lines(x = beta_for_plot, y = bound_bulk, col = 'red')
legend("topright", inset = 0.05, lty = 1, col = c('black', 'blue', 'red'),
       legend = c("Boundary", "Bulk", "Boundary and Bulk"))
## Doing the same with the mean surface
just_bound_exp <- ems_wave1$boundary$expectation$get_cov(test_points)[order(r_beta)]
just_bulk_exp <- ems_wave1$no_boundary$expectation[[out_name]]$get_cov(test_points)[order(r_beta)]
bound_bulk_exp <- ems_wave1$boundary_bulk$expectation$get_cov(test_points)[order(r_beta)]
plot(x = beta_for_plot, y = just_bound_exp,
     ylim = range(c(just_bound_exp, just_bulk_exp, bound_bulk_exp)),
     main = "Mean Emulator Uncertainty", ylab = "Uncertainty",
     xlab = "beta", type = 'l')
lines(x = beta_for_plot, y = just_bulk_exp, col = 'blue')
lines(x = beta_for_plot, y = bound_bulk_exp, col = 'red')
legend("topright", inset = 0.05, lty = 1, col = c('black', 'blue', 'red'),
       legend = c("Boundary", "Bulk", "Boundary and Bulk"))

## Create a large(-ish) grid for testing large-scale results
bgn <- 20
big_grid <- expand.grid(
  beta = seq(ranges$beta[[1]], ranges$beta[[2]], length.out = bgn),
  gamma = seq(ranges$gamma[[1]], ranges$gamma[[2]], length.out = bgn),
  omega = seq(ranges$omega[[1]], ranges$omega[[2]], length.out = bgn)
)
## Choosing a 'slice' on which to plot results
omega_point <- unique(big_grid$omega)[10]

## Might take a while to run (~5 minutes): evaluating emulator predictions across
# the grid for both variance emulators and mean emulators, in four different cases:
# 1) Using only bulk information (bulk)
# 2) Using only boundary information (bound)
# 3) Using bulk and boundary information (both)
# 4) Using only prior emulators with no Bayes Linear update (no)
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
    Vno = ems_wave1$no_boundary$variance[[out_name]]$o_em$get_cov(big_grid)
  )
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
## Reshaping the data.frames for ggplot-ing
var_df_reshape <- tidyr::pivot_longer(var_df, cols = !c(beta, gamma, omega))
exp_df_reshape <- tidyr::pivot_longer(exp_df, cols = !c(beta, gamma, omega))
## Actually plotting!
grid_plot(var_df_reshape, "V", omega_point, "Variance", wave0_output, 0.05)
grid_plot(var_df_reshape, "E", omega_point, "Variance", wave0_output, 0.05)
grid_plot(exp_df_reshape, "V", omega_point, "Mean", wave0_output, 0.05)
grid_plot(exp_df_reshape, "E", omega_point, "Mean", wave0_output, 0.05)

## Creating a 'target' for implausibility calculation
fake_point <- data.frame(beta = 0.1, gamma = 0.15, omega = omega_point)
out_value <- get_results(unlist(fake_point[1,], use.names = F), nreps = 500,
                         outs = c(out_name), times = t_point)
targ <- list(
  val = mean(out_value[[paste0(out_name, t_point)]]),
  sigma = sqrt(sd(out_value[[paste0(out_name, t_point)]])^2 +
                 (0.05*mean(out_value[[paste0(out_name, t_point)]])/3)^2)
)
## Create a finer grid on the omega-slice for implausibility calculations
grid_for_imp <- expand.grid(
  beta = seq(0, 1.5, by = 0.05),
  gamma = seq(0, 0.5, by = 0.01),
  omega = omega_point
)
imp_df <- cbind.data.frame(
  grid_for_imp,
  data.frame(
    Ibulk = ems_wave1$no_boundary$expectation[[out_name]]$implausibility(grid_for_imp, targ),
    Ibound = ems_wave1$boundary$expectation$implausibility(grid_for_imp, targ),
    Iboth = ems_wave1$boundary_bulk$expectation$implausibility(grid_for_imp, targ),
    Ino = abs(ems_wave1$no_boundary$expectation[[out_name]]$o_em$get_exp(grid_for_imp) - targ$val)/
      sqrt(ems_wave1$no_boundary$expectation[[out_name]]$o_em$get_cov(grid_for_imp) + targ$sigma^2)
  )
)
imp_df_reshape <- tidyr::pivot_longer(imp_df, cols = !c(beta, gamma, omega))
grid_plot(imp_df_reshape, "I", omega_point)

## Future Design
g_inverse <- inverse(integ, xmin = ranges$beta[[1]], xmax = ranges$beta[[2]],
                     theta = ems_wave1$no_boundary$expectation[[out_name]]$o_em$corr$hyper_p$theta)
g_norm <- g_inverse(ranges$beta[[2]])
## Generate a large LHD and perform warping relative to boundary beta = 0
test_design <- lhs::randomLHS(800*length(ranges), length(ranges))
test_lhs <- data.frame(t(apply(test_design, 1, function(x) {
  x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))
warp_lhs <- test_lhs
warp_lhs$beta <- sapply(warp_lhs$beta, function(x) {
  gin <- g_inverse(x)
  return(gin * 1.5/g_norm)
})
warp_imps <- ems_wave1$boundary_bulk$expectation$implausibility(warp_lhs, targ)
lhs_restrict <- warp_lhs[warp_imps <= 3,]

## Current method: take a maximin sample from the non-implausible warped points and choose reps
max_samp <- maximin_sample(lhs_restrict, 30, nms = names(ranges))
max_set <- design_subselect(max_samp, wave0_output[,names(ranges)],
                            ems_wave1$boundary_bulk$expectation,
                            ems_wave1$boundary_bulk$variance,
                            30, 300, ems_wave1$no_boundary$expectation[[out_name]]$o_em,
                            ranges)
imp_plot <- grid_plot(subset(imp_df_reshape, name == "Iboth"), "I", omega_point, "", max_set[,names(ranges)], 0.1) +
  geom_text(data = subset(max_set, abs(omega - omega_point) < 0.1),
                     aes(x = beta, y = gamma + 0.01, label = reps), size = 3) +
  geom_text(data = subset(max_set, abs(omega - omega_point) > 0.1),
            aes(x = beta, y = gamma + 0.01, label = reps), colour = 'grey40', size = 2)

## Doing it naively: rejection sampling, no boundary, no rep selection
imps_naive <- ems_wave1$boundary_bulk$expectation$implausibility(test_lhs, targ)
restrict_naive <- test_lhs[imps_naive <= 3,]
naive_prop <- maximin_sample(restrict_naive, 30)

### Running the new proposal points
wave1_new_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(max_set)), function(i) {
  get_results(max_set[i,names(ranges)], max_set[i, 'reps'], outs = c(out_name), times = t_point)
})) |> setNames(c(names(ranges), out_name))
wave1_new_output <- data.frame(wave1_new_results |> dplyr::group_by(beta, gamma, omega) |>
                                 dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))
wave1_naive_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(naive_prop)), function(i) {
  get_results(naive_prop[i,names(ranges)], 10, outs = c(out_name), times = t_point)
})) |> setNames(c(names(ranges), out_name))
wave1_naive_output <- data.frame(wave1_naive_results |> dplyr::group_by(beta, gamma, omega) |>
                                 dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))

new_ems <- create_boundary_ems(wave1_new_results, out_name, ranges, max_set$reps,
                               b_exp, b_cov, bb_exp, bb_cov, bb_imp)
naive_ems <- create_boundary_ems(wave1_naive_results, out_name, ranges, 10,
                                 b_exp, b_cov, bb_exp, bb_cov, bb_imp)

## Check implausibilities
imp_df_new <- cbind.data.frame(
  grid_for_imp,
  data.frame(
    Ibulk = new_ems$no_boundary$expectation[[out_name]]$implausibility(grid_for_imp, targ),
    Ibound = new_ems$boundary$expectation$implausibility(grid_for_imp, targ),
    Iboth = new_ems$boundary_bulk$expectation$implausibility(grid_for_imp, targ),
    Ino = abs(new_ems$no_boundary$expectation[[out_name]]$o_em$get_exp(grid_for_imp) - targ$val)/
      sqrt(new_ems$no_boundary$expectation[[out_name]]$o_em$get_cov(grid_for_imp) + targ$sigma^2)
  )
)
imp_df_new[is.na(imp_df_new)] <- 20
## Implausibilities at second wave are max between first and second at each point
imp_df_new$Ibulk <- purrr::map2_dbl(imp_df$Ibulk, imp_df_new$Ibulk, max)
imp_df_new$Ibound <- purrr::map2_dbl(imp_df$Ibound, imp_df_new$Ibound, max)
imp_df_new$Iboth <- purrr::map2_dbl(imp_df$Iboth, imp_df_new$Iboth, max)
imp_df_new$Ino <- purrr::map2_dbl(imp_df$Ino, imp_df_new$Ino, max)
imp_df_new_reshape <- tidyr::pivot_longer(imp_df_new, cols = !c(beta, gamma, omega))

grid_plot(subset(imp_df_new_reshape, name == "Iboth"), "I", omega_point, "New", fake_point, 0.1)

imp_df_naive <- cbind.data.frame(
  grid_for_imp,
  data.frame(
    Ibulk = naive_ems$no_boundary$expectation[[out_name]]$implausibility(grid_for_imp, targ),
    Ibound = naive_ems$boundary$expectation$implausibility(grid_for_imp, targ),
    Iboth = naive_ems$boundary_bulk$expectation$implausibility(grid_for_imp, targ),
    Ino = abs(naive_ems$no_boundary$expectation[[out_name]]$o_em$get_exp(grid_for_imp) - targ$val)/
      sqrt(naive_ems$no_boundary$expectation[[out_name]]$o_em$get_cov(grid_for_imp) + targ$sigma^2)
  )
)
imp_df_naive[is.na(imp_df_naive)] <- 20
imp_df_naive$Ibulk <- purrr::map2_dbl(imp_df$Ibulk, imp_df_naive$Ibulk, max)
imp_df_naive$Ibound <- purrr::map2_dbl(imp_df$Ibound, imp_df_naive$Ibound, max)
imp_df_naive$Iboth <- purrr::map2_dbl(imp_df$Iboth, imp_df_naive$Iboth, max)
imp_df_naive$Ino <- purrr::map2_dbl(imp_df$Ino, imp_df_naive$Ino, max)
imp_df_naive_reshape <- tidyr::pivot_longer(imp_df_naive, cols = !c(beta, gamma, omega))

grid_plot(subset(imp_df_naive_reshape, name == "Iboth"), "I", omega_point, "Naive", fake_point, 0.1)
