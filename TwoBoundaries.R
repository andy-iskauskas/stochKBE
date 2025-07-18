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
                                 SEIR_functions, bb_data, N_SEIR, c(0,0), c(3,4), out_index, t_point)

big_grid <- expand.grid(
  beta = seq(ranges$beta[1], ranges$beta[2], length.out = 20),
  eps = seq(ranges$eps[1], ranges$eps[2], length.out = 20)
)
for (nm in names(ranges)) {
  if (nm != "beta" && nm != "eps")
    big_grid[,nm] <- mean(ranges[[nm]])
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
    Vno = ems_wave1$no_boundary$variance[[out_name]]$o_em$get_cov(big_grid)
  )
)
var_df_reshape <- tidyr::pivot_longer(var_df, cols = !c(1:7))
ggplot(data = var_df_reshape[grepl("E", var_df_reshape$name),], aes(x = beta, y = eps, z = value)) +
  geom_contour_filled() +
  facet_wrap(vars(name), nrow = 2)