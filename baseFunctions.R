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
get_results <- function(params, nreps = 100, outs, times, raw = FALSE) {
  params <- unlist(params, use.names = FALSE)
  tseq <- 0:max(times)
  arra <- array(0, dim = c(max(tseq)+1, 3, nreps))
  for(i in 1:nreps) arra[,,i] <- gillespied(N,T=max(times) + 1 + 0.001,dt=1,th=params)
  if(raw) return(arra)
  collected <- list()
  for (i in 1:nreps) {
    relev <- c(arra[times+1, which(c("S", "I", "R") %in% outs), i])
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
                                bound_exp, bound_cov, bound_bulk_exp, bound_bulk_cov,
                                bound_bulk_imp) {
  data <- data.frame(data_raw |> dplyr::group_by(beta, gamma, omega) |>
                       dplyr::summarise(exp = mean(.data[[out_name]]), var = var(.data[[out_name]])))
  if (length(reps) == 1) reps <- rep(reps, nrow(data))
  no_bound_ems <- hmer::emulator_from_data(data_raw, out_name, ranges,
                                           emulator_type = "variance", 
                                           specified_priors = list(
                                             expectation = list(delta = c(0)),
                                             variance = list(delta = c(0))
                                           ))
  prior_var_em <- no_bound_ems$variance[[out_name]]$o_em
  prior_exp_em <- no_bound_ems$expectation[[out_name]]$o_em
  boundary_em <- hmer::Proto_emulator$new(
    ranges, out_name, bound_exp, bound_cov,
    em = prior_var_em, analytic = function(y) analytic_sd(y, t_point, out_index),
    val = 0
  )
  gillesp_var <- bound_cov(data[,names(ranges)], em = prior_var_em, full = TRUE) +
    purrr::map_dbl(seq_len(nrow(data)), ~prior_var_em$s_diag(data[.,], reps[.]))
  gillesp_var_inv <- tryCatch(chol2inv(chol(gillesp_var)), error = function(e) MASS::ginv(gillesp_var))
  gillesp_exp_diff <- data$var - b_exp(data[,names(ranges)], prior_var_em, function(y) analytic_sd(y, 15, 3))
  
  boundary_bulk_em <- hmer::Proto_emulator$new(
    ranges, out_name, bound_bulk_exp, bound_bulk_cov, pre_em = prior_var_em,
    b_em = boundary_em, dat = data[,names(ranges)],
    binv = gillesp_var_inv, bmod = gillesp_exp_diff,
    analytic = function(y) analytic_sd(y, t_point, out_index),
    val = 0
  )
  boundary_em_mean <- hmer::Proto_emulator$new(
    ranges, out_name, bound_exp, bound_cov, em = prior_exp_em,
    analytic = function(y) analytic_mean(y, t_point, out_index), val = 0
  )
  
  gillesp_e_var <- b_cov(data[,names(ranges)], em = prior_exp_em, full = TRUE) +
    purrr::map_dbl(seq_len(nrow(data)), ~boundary_em$get_exp(data[.,])/reps[.])
  gillesp_e_var_inv <- tryCatch(chol2inv(chol(gillesp_e_var)), error = function(e) MASS::ginv(gillesp_e_var))
  gillesp_e_exp_diff <- data$exp - b_exp(data[,names(ranges)], prior_exp_em, function(y) analytic_mean(y, 15, 3))
  
  boundary_bulk_em_mean <- hmer::Proto_emulator$new(
    ranges, out_name, bound_bulk_exp, bound_bulk_cov, pre_em = prior_exp_em,
    implausibility_func = bound_bulk_imp, v_em = boundary_bulk_em,
    b_em = boundary_em_mean, dat = data[,names(ranges)],
    binv = gillesp_e_var_inv, bmod = gillesp_e_exp_diff,
    analytic = function(y) analytic_mean(y, t_point, out_index), val = 0
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
## Calculates w(X,X'); as per Binois (2019)
w_func <- function(x, xp, theta, ranges) {
  (sqrt(2*pi/4) * theta)^length(x) * prod(
    purrr::map_dbl(seq_along(x), function(i) {
      exp(-(x[i]-xp[i])^2/(2*theta^2)) * (2*pnorm(sqrt(2) * (2*ranges[[i]][2]-x[i]-xp[i])/(sqrt(2)*theta)) -
                                            2*pnorm(sqrt(2) * (2*ranges[[i]][1]-x[i]-xp[i])/sqrt(2)*theta))
    })
  )
}
## Finds the variance reduction under the introduction of a new point, x'
new_point_score <- function(point, points, em_exp, em_var, ranges, kinv, wmat) {
  r_val <- max(0, em_var$get_exp(point)) # r(x)
  sig <- em_exp$get_cov(point) + r_val #sigma2(x)
  k_val <- em_exp$get_cov(point, points, full = TRUE) #k(x)
  em_theta <- em_exp$add_args$pre_em$corr$hyper_p$theta
  w_vec <- purrr::map_dbl(seq_len(nrow(points)), function(i) w_func(unlist(point), unlist(points[i,]), em_theta, ranges)) # w(x)
  (k_val %*% kinv %*% wmat %*% kinv %*% t(k_val) + w_func(unlist(point), unlist(point), em_theta, ranges) - 2*w_vec %*% kinv %*% t(k_val))/sig
}
## Finds the variance reduction under the addition of a replicate to an existing point
new_rep_score <- function(points, reps, index, em_var, kinv, wmat) {
  r_val <- max(em_var$get_exp(points[index,]), 0)
  kout <- outer(kinv[index,,drop=TRUE], kinv[,index,drop=TRUE], "*")
  trace_val <- sum(diag(kout %*% wmat))
  quot <- reps[index]*(reps[index]+1)*r_val - kinv[index,index]
  return(trace_val/quot)
}

## Chooses an 'optimal' design
# Given a set of (non-implausible) points, progressively adds points to a candidate
# set or adds a rep to an existing point in the candidate set. Returns the collection
# of points, along with a column denoting how many reps are to be run at each point.
design_subselect <- function(points, data, em_exp, em_var, pt_max, rep_max, prior_em, ranges) {
  pt_var_exp <- em_exp$get_cov(points)
  pt_var_var <- em_var$get_exp(points)
  pt_var_var[pt_var_var < 0] <- 1e-6
  find_next_point <- function(candidates, data, em_exp, em_var, ranges, kinv, wmat) {
    point_scores <- purrr::map_dbl(seq_len(nrow(candidates)), function(i) {
      new_point_score(candidates[i,], data, em_exp, em_var, ranges, kinv, wmat)
    })
    return(list(val = max(point_scores), index = which.max(point_scores)))
  }
  find_next_rep <- function(points, reps, em_var, kinv, wmat) {
    rep_scores <- purrr::map_dbl(seq_len(nrow(points)), function(i) {
      new_rep_score(points, reps, i, em_var, kinv, wmat)
    })
    return(list(val = max(rep_scores), index = which.max(rep_scores)))
  }
  if (nrow(points) == pt_max) {
    current_pts <- seq_len(nrow(points))
    current_reps <- rep(2, nrow(points))
  }
  else {
    current_pts <- c()
    current_reps <- c()
  }
  fixed_indices <- seq_len(nrow(points))
  failsafe <- 0
  while (sum(current_reps) < rep_max && failsafe < 500) {
    current_points <- points[current_pts,]
    current_data <- rbind.data.frame(data, current_points)
    if (length(current_reps) != 0) {
      rep_k_mat <- em_exp$get_cov(current_points, full = TRUE)
      rep_k_inv <- MASS::ginv(rep_k_mat)
      rep_w_mat <- matrix(do.call('rbind', purrr::map(seq_len(nrow(current_points)), function(i) {
        purrr::map_dbl(seq_len(nrow(current_points)), function(j) {
          w_func(unlist(current_points[i,]), unlist(current_points[j,]),
                 prior_em$corr$hyper_p$theta, ranges)
        })
      })), nrow = nrow(current_points), byrow = TRUE)
      rep_suggest <- find_next_rep(current_points, current_reps, em_var, rep_k_inv, rep_w_mat)
    }
    else
      rep_suggest <- NULL
    if (length(current_pts) < pt_max) {
      pt_k_mat <- em_exp$get_cov(current_data, full = TRUE)
      pt_k_inv <- MASS::ginv(pt_k_mat)
      pt_w_mat <- matrix(do.call('rbind', purrr::map(seq_len(nrow(current_data)), function(i) {
        purrr::map_dbl(seq_len(nrow(current_data)), function(j) {
          w_func(unlist(current_data[i,]), unlist(current_data[j,]),
                 prior_em$corr$hyper_p$theta, ranges)
        }
        )})), nrow = nrow(current_data), byrow = TRUE)
      if (length(current_pts) == 0)
        pt_suggest <- find_next_point(points, data, em_exp, em_var, ranges, pt_k_inv, pt_w_mat)
      else
        pt_suggest <- find_next_point(points[-current_pts,], current_data, em_exp, em_var, ranges, pt_k_inv, pt_w_mat)
    }
    else
      pt_suggest <- NULL
    if (is.null(pt_suggest) || (!is.null(rep_suggest) && rep_suggest$val > pt_suggest$val)) {
      current_reps[rep_suggest$index] <- current_reps[rep_suggest$index] + 1
    }
    else {
      current_pts <- c(current_pts, fixed_indices[pt_suggest$index])
      fixed_indices <- fixed_indices[-pt_suggest$index]
      current_reps <- c(current_reps, 2)
    }
    failsafe <- failsafe + 1
  }
  selected_pts <- points[current_pts,]
  pts_with_reps <- cbind.data.frame(selected_pts, current_reps) |> setNames(c(names(selected_pts), "reps"))
  return(pts_with_reps)
}
