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
wave0_output <- data.frame(wave0_results |> dplyr::group_by(across(all_of(names(ranges)))) |>
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
test_design <- lhs::randomLHS(50*length(ranges), length(ranges))
test_lhs <- data.frame(t(apply(test_design, 1, function(x) {
  x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))
warp_lhs <- test_lhs
warp_lhs$beta <- sapply(warp_lhs$beta, function(x) {
  gin <- g_inverse(x)
  return(gin * 1.5/g_norm)
})
warp_imps <- ems_wave1$boundary_bulk$expectation$implausibility(warp_lhs, targ)
#lhs_restrict <- warp_lhs[warp_imps <= 3,]
# Use all points, ignoring implausibility
lhs_restrict <- warp_lhs

## Current method: take a maximin sample from the non-implausible warped points and choose reps
max_samp <- maximin_sample(lhs_restrict, 30, nms = names(ranges))
# max_set <- design_subselect(max_samp, wave0_output[,names(ranges)],
#                             ems_wave1$boundary_bulk$expectation,
#                             ems_wave1$boundary_bulk$variance,
#                             30, 300, ems_wave1$no_boundary$expectation[[out_name]]$o_em,
#                             ranges)

max_set <- design_subselect(lhs_restrict, wave0_output[,names(ranges)],
                            ems_wave1$boundary_bulk$expectation,
                            ems_wave1$boundary_bulk$variance,
                            30, 300, ems_wave1$no_boundary$expectation[[out_name]]$o_em,
                            ranges, store_order = TRUE, verbose = TRUE)


## Checking 2d proposal
test_design_2d <- lhs::randomLHS(100*length(ranges), 2)
test_lhs_2d <- data.frame(t(apply(test_design_2d, 1, function(x) {
  x * purrr::map_dbl(ranges, diff)[1:2] + purrr::map_dbl(ranges, ~.[[1]])[1:2]
}))) |> setNames(c('beta', 'gamma'))
warp_lhs_2d <- test_lhs_2d
warp_lhs_2d$beta <- sapply(warp_lhs_2d$beta, function(x) {
  gin <- g_inverse(x)
  return(gin * 1.5/g_norm)
})
warp_lhs_2d$omega <- omega_point
warp_imps <- ems_wave1$boundary_bulk$expectation$implausibility(warp_lhs_2d, targ)
source("modelFunctions.R")
max_set_2d <- design_subselect(warp_lhs_2d, wave0_output[,names(ranges)],
                            ems_wave1$boundary_bulk$expectation,
                            ems_wave1$boundary_bulk$variance,
                            30, 300, ems_wave1$no_boundary$expectation[[out_name]]$o_em,
                            ranges, ntoadd = 1, store_order = TRUE,
                            verbose = TRUE, return_scores = TRUE, fixed_vals = list(omega = omega_point))

plot(x = 1:60, y = max_set_2d$rep_scores, type = 'l', ylim = c(0, max(c(max_set_2d$rep_scores, max_set_2d$pt_scores))),
     main = "Scores for new point vs extra rep", xlab = "Proposal Index", ylab = "Score")
lines(x = 1:60, y = max_set_2d$pt_scores, col = "blue")
legend("topright", legend = c("Extra Rep", "New Point"), col = c("black", "blue"), lty = 1)

## Animating point proposal
## Not strictly necessary, but quite cool
## Requires gganimate, gifski, and scales
library(gganimate)
library(gifski)
library(scales)
create_transition_data_frame <- function(data, order, nadd = 1) {
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
trans_df <- create_transition_data_frame(max_set$points, max_set$order)
imp_df_rename <- setNames(imp_df_reshape, c("b", "g", "o", "n", "v"))
imp_breaks_mod <- imp_breaks
imp_breaks_mod[length(imp_breaks_mod)] <- 1000
anim <- ggplot(data = trans_df, aes(x = beta, y = gamma)) +
  geom_raster(data = subset(imp_df_rename, n == "Iboth" & o == omega_point), aes(x = b, y = g, fill = v), interpolate = TRUE) +
  scale_fill_gradientn(colours = redgreen, limits = c(0, 1000), breaks = imp_breaks_mod, values = rescale(imp_breaks_mod, c(0,1), c(0, 1000)), name = "I", labels = imp_names) +
  geom_contour(data = subset(imp_df_rename, n == "Iboth" & o == omega_point), aes(x = b, y = g, z = v), breaks = c(0, 3, Inf), colour = "black") +
  geom_point(data = subset(trans_df, abs(omega - omega_point) < 0.1), size = 2) +
  geom_point(data = subset(trans_df, abs(omega - omega_point) >= 0.1), colour = "grey40", size = 2) +
  geom_text(data = subset(trans_df, abs(omega - omega_point) < 0.1), aes(y = gamma + 0.01, label = reps), size = 5) +
  geom_text(data = subset(trans_df, abs(omega - omega_point) >= 0.1), aes(y = gamma + 0.01, label = reps), size = 4, colour = "grey40") +
  transition_time(time) +
  ease_aes("linear") +
  theme(legend.position = "none") +
  ggtitle("Repetitions Placed: {frame_time}")
animate(anim, nframes = 300, end_pause = 25, height = 800, width = 800)
anim_save("NewProposalAnimation.gif")

exp_df_rename <- setNames(exp_df_reshape, c('b', 'g', 'o', 'n', 'v'))
anim_alt <- ggplot(data = trans_df, aes(x = beta, y = gamma)) +
  geom_raster(data = subset(exp_df_rename, n == "Vboth" & o == omega_point),
              aes(x = b, y = g, fill = v), interpolate = TRUE) +
  # geom_contour(data = subset(imp_df_rename, n == "Iboth" & o == omega_point),
  #              aes(x = b, y = g, z = v), breaks = c(0, 3, Inf), colour = "white") +
  scale_fill_viridis(name = "Var") +
  geom_contour(data = subset(imp_df_rename, n == "Iboth" & o == omega_point), aes(x = b, y = g, z = v), breaks = c(0, 3, Inf), colour = "black") +
  geom_point(data = subset(trans_df, abs(omega - omega_point) < 0.1), size = 2, colour = "white") +
  geom_point(data = subset(trans_df, abs(omega - omega_point) >= 0.1), colour = "grey80", size = 2) +
  geom_text(data = subset(trans_df, abs(omega - omega_point) < 0.1), aes(y = gamma + 0.01, label = reps), size = 5, colour = "white") +
  geom_text(data = subset(trans_df, abs(omega - omega_point) >= 0.1), aes(y = gamma + 0.01, label = reps), size = 4, colour = "grey80") +
  transition_time(time) +
  ease_aes("linear") +
  theme(legend.position = "none") +
  ggtitle("Repetitions Placed: {frame_time}")
animate(anim_alt, nframes = 300, end_pause = 25, height = 800, width = 800)
anim_save("AltProposalAnimation.gif")

trans_df_2d <- create_transition_data_frame(max_set_2d$points, max_set_2d$order, nadd = 1)
exp_df_rename <- setNames(exp_df_reshape, c("b", "g", "o", "n", "v"))
anim_2d <- ggplot(data = trans_df_2d, aes(x = beta, y = gamma)) +
  geom_raster(data = subset(exp_df_rename, n == "Eboth" & o == omega_point),
              aes(x = b, y = g, fill = v), interpolate = TRUE) +
  # geom_contour(data = subset(imp_df_rename, n == "Iboth" & o == omega_point), aes(x = b, y = g, z = v), breaks = c(0, 3, Inf), colour = "white") +
  scale_fill_viridis(name = "Exp") +
  geom_point(size = 2, colour = "grey90") +
  geom_text(aes(y = gamma + 0.01, label = reps), size = 5, colour = "grey90") +
  transition_time(time) +
  ease_aes("linear") +
  #theme(legend.position = "none") +
  ggtitle("Repetitions Placed: {frame_time}")
## Testing animation
#animate(anim_2d, nframes = 300, end_pause = 25)
## For exporting
animate(anim_2d, nframes = 300, end_pause = 25, height = 800, width = 800)
anim_save("NewProposal2dExp.gif")

imp_plot <- grid_plot(subset(imp_df_reshape, name == "Iboth"), "I", omega_point, "", max_set$points[,names(ranges)], 0.1) +
  geom_text(data = subset(max_set$points, abs(omega - omega_point) < 0.1),
                     aes(x = beta, y = gamma + 0.01, label = reps), size = 3) +
  geom_text(data = subset(max_set$points, abs(omega - omega_point) > 0.1),
            aes(x = beta, y = gamma + 0.01, label = reps), colour = 'grey40', size = 2)
var_plot <- grid_plot(subset(exp_df_reshape, name == "Vboth"), "V", omega_point, "", max_set$points[,names(ranges)], 0.1,
                      breaks = c(0, 10, 20, 50, 500, 1000, 5000, Inf)) +
  geom_text(data = subset(max_set$points, abs(omega - omega_point) < 0.1),
            aes(x = beta, y = gamma + 0.01, label = reps), size = 3) +
  geom_text(data = subset(max_set$points, abs(omega - omega_point) > 0.1),
            aes(x = beta, y = gamma + 0.01, label = reps), colour = "grey40", size = 2)
suggested_imps <- ems_wave1$boundary_bulk$expectation$implausibility(max_set[,names(ranges)], targ)
suggested_vars <- ems_wave1$boundary_bulk$expectation$get_cov(max_set[,names(ranges)])

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


var_vals <- new_ems$boundary_bulk$expectation$get_cov(imp_df[,names(ranges)])
new_var_df <- cbind.data.frame(imp_df[,names(ranges)], var_vals) |> setNames(c(names(ranges), "V"))
new_var_df[new_var_df$V < 0, "V"] <- 1e-6
ggplot(data = new_var_df, aes(x = beta, y = gamma)) +
  geom_contour_filled(aes(z = sqrt(V)), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf))
