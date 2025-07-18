###############################
# Emulation of SIR Model - 2d #
###############################

source("modelFunctions.R")
source("plotting.R")
library(lhs)
library(tidyr)
set.seed(42)

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
set.seed(1)
training_points <- data.frame(t(apply(
  lhs::optimumLHS(10*length(ranges), length(ranges)),
  1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }
))) |> setNames(names(ranges))
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
small_test_grid <- expand.grid(beta = seq(0, 1.5, length.out = 20), gamma = seq(0, 0.5, length.out = 20))
## Create design for next wave of emulation
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

bgn <- 20
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
grid_plot(exp_df_reshape, "V", NULL, "Mean", wave0_output, 0.05)

final_vars <- cbind.data.frame(small_test_grid, mean_em_var(small_test_grid, ems_wave1$no_boundary$expectation$I$o_em,
                                                            this_var_em, new_design$points[,1:2], new_design$points$reps)) |>
  setNames(c("beta", "gamma", "Var"))
compare_vars <- cbind.data.frame(final_vars, exp_df$Vboth) |>
  setNames(c("beta", "gamma", "after", "before"))
compare_var_reshape <- tidyr::pivot_longer(compare_vars, cols = !c(beta, gamma))
ggplot(data = compare_var_reshape, aes(x = beta, y = gamma, fill = value)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis(name = "Var", breaks = c(0, 100, 200, 300, 400, 500, 1000, 2000, 5000)) +
  facet_grid(rows = vars(name))

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
anim_save("PropAnim2dAlt.gif", animation = anim_alt)
save.image(file = "May13.RData")


### Training new emulator and comparing to other proposal methods
# Three methods of point proposal are considered: one using the improved design with
# the corresponding suggested repetitions at each point (new_ems); one using a standard
# Latin hypercube design with uniformly many repetitions at each point (basic_ems), and
# one using the new design points but allocating equal numbers of repetitions to all
# points (unif_rep_ems).

# Setup for the new design
wave1_results <- do.call('rbind.data.frame', purrr::map(seq_len(nrow(new_design$points)), function(i) {
  get_results(unlist(new_design$points[i,c('beta', 'gamma')], use.names = FALSE), nreps = new_design$points[i,"reps"], outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
wave1_output <- data.frame(wave1_results |> dplyr::group_by(across(all_of(names(ranges)))) |>
                             dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))
new_ems <- create_boundary_ems(wave1_results, out_name, ranges, reps,
                               b_exp, b_cov, bb_exp, bb_cov, bb_imp)

# Setup for the basic design
basic_points <- lhs::maximinLHS(20, 2)
basic_scaled <- data.frame(t(apply(basic_points, 1, function(x) {
  x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))
all_basic <- rbind.data.frame(training_points, basic_scaled)
basic_added <- do.call("rbind.data.frame", purrr::map(seq_len(nrow(basic_scaled)), function(i) {
  get_results(unlist(basic_scaled[i,], use.names = FALSE), nreps = 10, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
all_points_basic <- rbind.data.frame(wave0_results, basic_added)
basic_ems <- create_boundary_ems(all_points_basic, out_name, ranges, rep(10, 40),
                                 b_exp, b_cov, bb_exp, bb_cov, bb_imp)

# Setup for the new design with uniform reps
unif_rep_res <- do.call("rbind.data.frame", purrr::map(seq_len(nrow(new_design$points)), function(i) {
  get_results(unlist(new_design$points[i,c("beta", "gamma")], use.names = FALSE), nreps = 10, outs = c(out_name), times = 15)
})) |> setNames(c(names(ranges), out_name))
unif_rep_ems <- create_boundary_ems(unif_rep_res, out_name, ranges, rep(10, 40),
                                    b_exp, b_cov, bb_exp, bb_cov, bb_imp)

## Comparing predictive variance across the space
new_vars <- new_ems$boundary_bulk$expectation$get_cov(small_test_grid)
new_vars[new_vars < 0] <- 1e-6
basic_vars <- basic_ems$boundary_bulk$expectation$get_cov(small_test_grid)
basic_vars[basic_vars < 0] <- 1e-6
unif_rep_vars <- unif_rep_ems$boundary_bulk$expectation$get_cov(small_test_grid)
# Mean predictive variance
mean(new_vars)
mean(basic_vars)
mean(unif_rep_vars)
mean(basic_vars)/mean(new_vars)

## Plotting the results: used here are the three emulator sets trained above with
# the different proposals, as well as the original wave 1 stoch KBE emulators
all_var_df <- cbind.data.frame(cbind.data.frame(small_test_grid, new_vars), cbind.data.frame(basic_vars, unif_rep_vars, exp_df$Vboth)) |>
  setNames(c('beta', 'gamma', 'New', 'Old', "Uniform", "Naive"))
reshape_var_df <- tidyr::pivot_longer(all_var_df, cols = !c(beta, gamma))
reshape_var_df$name <- factor(reshape_var_df$name, levels = c("Old", "Naive", "Uniform", "New"))

ggplot(data = reshape_var_df, aes(x = beta, y = gamma, z = value)) +
  geom_contour_filled(breaks = c(0, 50, 100, 150, 200, 250, 300, 400, 500, 600, 1000, 5000)) +
  scale_fill_viridis(discrete = TRUE, name = "Variance") +
  facet_wrap(vars(name), nrow = 2)
