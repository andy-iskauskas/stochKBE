#######################
#  GENERIC FUNCTIONS  #
#######################

## Currently needs the github version of hmer
# devtools::install_github("andy-iskauskas/hmer")

library(purrr)
library(dplyr)
library(hmer)
library(MASS)
## Gillespie algorithm for obtaining model realisations
# N is a list containing initial compartment numbers,
# pre- and post-transition matrices, and hazard function
# Returns the final compartment numbers at for times 1,...,T
gillespied=function (N, T=400, dt=1, ...)
{
  tt=0
  n=T%/%dt
  x=N$M
  S=t(N$Post-N$Pre)
  u=nrow(S)
  v=ncol(S)
  xmat=matrix(0,ncol=u,nrow=n)
  i=1
  target=0
  repeat {
    h=N$h(x, tt, ...)
    h0=sum(h)
    if (h0<1e-10)
      tt=1e99
    else
      tt=tt+rexp(1,h0)
    while (tt>=target) {
      xmat[i,]=x
      i=i+1
      target=target+dt
      if (i>n)
        return(xmat)
    }
    j=sample(v,1,prob=h)
    x=x+S[,j]
  }
}

## Obtain sanitised results from Gillespie algorithm
# Takes params (a collection of input parameters), nreps (number of
# repetitions at the parameter set), outs (the output compartments
# returned), and times (the times to record). raw gives the option of
# returning the full, unsanitised, data.
get_results <- function(params, obj, nreps = 100, outs, times, raw = FALSE) {
  params <- unlist(params, use.names = FALSE)
  tseq <- 0:max(times)
  arra <- array(0, dim = c(max(tseq)+1, ncol(obj$Pre), nreps))
  for(i in 1:nreps) arra[,,i] <- gillespied(obj ,T=max(times) + 1 + 0.001,dt=1,th=params)
  if(raw) return(arra)
  collected <- list()
  for (i in 1:nreps) {
    relev <- c(arra[times+1, which(names(obj$M) %in% outs), i])
    names <- unlist(purrr::map(outs, ~paste0(., times, sep = "")))
    relev <- setNames(relev, names)
    collected[[i]] <- relev
  }
  input_dat <- setNames(data.frame(matrix(rep(params, nreps), ncol = length(params), byrow = TRUE)), names(params))
  return(cbind(input_dat, do.call('rbind', collected)))
}

## Create boundary emulators
# Takes data_raw (the full collection of data from Gillespie, eg), out_name (the
# name of the output to emulate), ranges (parameter ranges), and reps (either a
# single value if all parameter sets used the same number of realisations, or a
# vector of numerics of length equal to the number of parameter sets). Returns
# a list of emulators: prior, bulk only, boundary only, and bulk-boundary.
create_boundary_ems <- function(data_raw, out_name, ranges, reps,
                                analytics, bb_data, model,
                                vals = c(0), indices = c(1), out_index, t) {
  data <- data.frame(data_raw |> dplyr::group_by(across(all_of(names(ranges)))) |>
                       dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))
  if (length(reps) == 1) reps <- rep(reps, nrow(data))
  no_bound_ems <- hmer::emulator_from_data(data_raw, out_name, ranges,
                                           emulator_type = "variance", 
                                           # specified_priors = list(
                                           #   expectation = list(delta = c(0.01)),
                                           #   variance = list(delta = c(0.01))
                                           # )
                                           )
  prior_var_em <- no_bound_ems$variance[[out_name]]$o_em
  prior_exp_em <- no_bound_ems$expectation[[out_name]]$o_em
  boundary_em <- hmer::Proto_emulator$new(
    ranges, out_name, analytics$b_exp, analytics$b_cov,
    em = prior_var_em, analytic = function(y) analytics$analytic_sd(y, t, out_index, init_vals = model$M),
    vals = vals, indices = indices
  )
  gillesp_var <- analytics$b_cov(data[,names(ranges)], em = prior_var_em, full = TRUE, vals = vals, indices = indices) +
    purrr::map_dbl(seq_len(nrow(data)), ~prior_var_em$s_diag(data[.,], reps[.]))
  gillesp_var_inv <- tryCatch(chol2inv(chol(gillesp_var)), error = function(e) MASS::ginv(gillesp_var))
  gillesp_exp_diff <- data$var - analytics$b_exp(data[,names(ranges)], prior_var_em, function(y) analytics$analytic_sd(y, t, out_index, init_vals = model$M),
                                                 vals = vals, indices = indices)
  
  boundary_bulk_em <- hmer::Proto_emulator$new(
    ranges, out_name, bb_data$bb_exp, bb_data$bb_cov, pre_em = prior_var_em,
    b_em = boundary_em, dat = data[,names(ranges)],
    binv = gillesp_var_inv, bmod = gillesp_exp_diff,
    analytic = function(y) analytics$analytic_sd(y, t, out_index, init_vals = model$M),
    vals = vals, indices = indices, bcov = analytics$b_cov
  )
  boundary_em_mean <- hmer::Proto_emulator$new(
    ranges, out_name, analytics$b_exp, analytics$b_cov, em = prior_exp_em,
    analytic = function(y) analytics$analytic_mean(y, t, out_index, init_vals = model$M),
    vals = vals, indices = indices
  )
  
  gillesp_e_var <- analytics$b_cov(data[,names(ranges)], em = prior_exp_em, full = TRUE, vals = vals, indices = indices) +
    purrr::map_dbl(seq_len(nrow(data)), ~boundary_em$get_exp(data[.,])/reps[.])
  gillesp_e_var_inv <- tryCatch(chol2inv(chol(gillesp_e_var)), error = function(e) MASS::ginv(gillesp_e_var))
  gillesp_e_exp_diff <- data$exp - analytics$b_exp(data[,names(ranges)], prior_exp_em, function(y) analytics$analytic_mean(y, t, out_index, init_vals = model$M),
                                                   vals = vals, indices = indices)
  
  boundary_bulk_em_mean <- hmer::Proto_emulator$new(
    ranges, out_name, bb_data$bb_exp, bb_data$bb_cov, pre_em = prior_exp_em,
    implausibility_func = bb_data$bb_imp, v_em = boundary_bulk_em,
    b_em = boundary_em_mean, dat = data[,names(ranges)],
    binv = gillesp_e_var_inv, bmod = gillesp_e_exp_diff,
    analytic = function(y) analytics$analytic_mean(y, t, out_index, init_vals = data$M),
    vals = vals, indices = indices, bcov = analytics$b_cov
  )
  
  return(
    list(no_boundary = no_bound_ems,
         boundary = list(variance = boundary_em, expectation = boundary_em_mean),
         boundary_bulk = list(variance = boundary_bulk_em, expectation = boundary_bulk_em_mean)
    )
  )
}

## Warping functions
## g^-1(x)
# This assumes an exponential squared, separable correlation structure
integ <- function(y, theta, xmin) {
  y - xmin + sqrt(pi/2) * theta * (pnorm(2*xmin/theta) - pnorm(2*y/theta))
}
## Inverter for g^-1(x)
inverse = function (f, lower = -5, upper = 10, xmin, xmax, theta) {
  function (y) uniroot((function (x) integ(x, theta, xmin) - y), lower = lower, upper = upper)[[1]]
}

## Scoring functions for variance reduction
# Calculate the mean emulator variance across the space, using a representative
# collection of points.
mean_em_var <- function(pt, pre_em, var_em, data, reps, boundary_col = 1, boundary_val = 0) {
  em_prior_var <- pre_em$u_sigma^2
  datamutate <- (data |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
  ptmutate <- (pt |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
  datarminus <- pre_em$get_cov(datamutate, full = TRUE)/em_prior_var
  term1 <- em_prior_var * (1-r1(pt, pre_em)^2)
  term2a <- R1(pt, data, pre_em) * pre_em$get_cov(ptmutate, datamutate, full = TRUE)
  term2c <- R1(data, pt, pre_em) * pre_em$get_cov(datamutate, ptmutate, full = TRUE)
  var_em_vals <- var_em$get_exp(data)
  var_em_vals[var_em_vals < 0] <- 1e-6
  term2binv <- R1(data, data, pre_em) * pre_em$get_cov(datamutate, datamutate, full = TRUE) +
    diag(c(1/reps * var_em_vals))
  term2b <- tryCatch(chol2inv(chol(term2binv)),
                     error = function(e) MASS::ginv(term2binv))
  term2 <- term2a %*% term2b %*% term2c
  complete <- term1 - diag(term2)
  complete[complete < 0] <- 1e-6
  return(complete)
}
# Calculate the score due to including a new design point
new_point_score <- function(point, points, pre_em, var_em, reps, grid, ntoadd = 1) {
  new_data <- rbind.data.frame(points, point)
  if (is.null(point))
    new_reps <- reps
  else
    new_reps <- c(reps, ntoadd)
  vars <- mean(mean_em_var(grid, pre_em, var_em, new_data, new_reps))
  return(vars)
}
# Calculate the score due to adding a repetition at an existing design point
new_rep_score <- function(index, points, pre_em, var_em, reps, grid, ntoadd = 1) {
  new_reps <- reps
  new_reps[index] <- new_reps[index]+ntoadd
  vars <- mean(mean_em_var(grid, pre_em, var_em, points, new_reps))
  return(vars)
}

## Chooses an 'optimal' design
# Given a set of data points, progressively adds points to a candidate
# set or adds a rep to an existing point in the candidate set. Calculations of
# expected improvement are based on evaluating the mean emulator variance
# over a grid of points, testgrid: the more dense the grid, the more accurate the
# estimate, but the more computationally intensive it is. Returns the collection
# of points, along with a column denoting how many reps are to be run at each point.
design_subselect <- function(data, pre_em, var_em, reps, ranges, testgrid,
                              rep_max, pt_max, ntoadd = 1,
                              store_order = FALSE, verbose = FALSE,
                              return_scores = FALSE) {
  if (store_order) order_vector <- c()
  if (return_scores) {
    r_scores <- c()
    p_scores <- c()
  }
  find_next_point <- function(data, pre_em, var_em, reps, ranges, ntoadd) {
    opt_func <- function(x) {
      x_mod <- data.frame(matrix(x, nrow = 1)) |> setNames(names(ranges))
      x_mod <- x_mod[,names(ranges), drop = FALSE]
      new_point_score(x_mod, data, pre_em, var_em, reps, testgrid, ntoadd)
    }
    optimised <- optim(purrr::map_dbl(ranges, ~runif(1, .[[1]], .[[2]])), opt_func,
                       lower = purrr::map_dbl(ranges, ~.[[1]]+0.01*diff(.)),
                       upper = purrr::map_dbl(ranges, ~.[[2]]-0.01*diff(.)),
                       method = "L-BFGS-B", control = list(trace = FALSE))
    this_val <- optimised$value
    this_pt <- data.frame(matrix(optimised$par, nrow = 1)) |> setNames(names(ranges))
    dists <- apply(data, 1, function(x) {
      sum((x-this_pt)^2)
    })
    if (this_val < 0) this_val <- Inf
    if (any(dists < 1e-6)) this_val <- NaN
    return(list(val = this_val, point = this_pt))
  }
  find_next_rep <- function(data, pre_em, var_em, reps, ntoadd) {
    rep_vals <- purrr::map_dbl(seq_len(nrow(data)), function(i) {
      new_rep_score(i, data, pre_em, var_em, reps, testgrid, ntoadd)
    })
    return(list(val = min(rep_vals), index = which.min(rep_vals)))
  }
  failsafe <- 0
  while(sum(reps) < rep_max && failsafe < 1000) {
    rep_suggest <- find_next_rep(data, pre_em, var_em, reps, ntoadd)
    if (nrow(data) < pt_max)
      pt_suggest <- find_next_point(data, pre_em, var_em, reps, ranges, ntoadd)
    else
      pt_suggest <- NULL
    if (verbose) {
      print_str <- paste0("Proposal ", failsafe, ":")
      if (!is.null(pt_suggest)) {
        print_str <- paste0(print_str, " Point score ", signif(pt_suggest$val, 4), ";")
      }
      print_str <- paste0(print_str, " Rep score ", signif(rep_suggest$val, 4))
      if (!is.null(pt_suggest) && !is.nan(pt_suggest$val)) {
        if (pt_suggest$val < rep_suggest$val)
          print_str <- paste0(print_str, " - New point chosen;")
        else
          print_str <- paste0(print_str, " - Extra rep chosen;")
      }
      else
        print_str <- paste0(print_str, " - Extra rep chosen;")
      print(print_str)
    }
    if (return_scores) {
      r_scores <- c(r_scores, rep_suggest$val)
      if (is.null(pt_suggest)) p_scores <- c(p_scores, Inf)
      else p_scores <- c(p_scores, pt_suggest$val)
    }
    if (is.null(pt_suggest) || is.nan(pt_suggest$val) || rep_suggest$val < pt_suggest$val) {
      reps[rep_suggest$index] <- reps[rep_suggest$index] + ntoadd
      if (store_order) order_vector <- c(order_vector, rep_suggest$index)
    }
    else {
      data <- rbind.data.frame(data, pt_suggest$point)
      reps <- c(reps, ntoadd)
      if (store_order) order_vector <- c(order_vector, nrow(data))
    }
    failsafe <- failsafe + 1
  }
  pts_with_reps <- cbind.data.frame(data, reps) |> setNames(c(names(data), "reps"))
  if (store_order) {
    if (return_scores)
      return(list(points = pts_with_reps, order = order_vector, rep_scores = r_scores, pt_scores = p_scores))
    return(list(points = pts_with_reps, order = order_vector))
  }
  return(pts_with_reps)
}
