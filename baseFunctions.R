#######################
#  GENERIC FUNCTIONS  #
#######################

## Currently needs the github version of hmer
# devtools::install_github("andy-iskauskas/hmer")

library(purrr)
library(dplyr)
library(hmer)
library(MASS)

## Optional: future and optimParallel for parallelisation of proposal
has_par_optim <- suppressWarnings(require(optimParallel))
has_future <- suppressWarnings(require(future))

### HELPER FUNCTIONS
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

## Partition Inverse for augmented matrix
# Given an nxn covariance matrix with pre-calculated inverse ainv, and an 
# additional point to be included, calculates the inverse of the augmented
# (n+1)x(n+1) covariance matrix, where b is the covariance between the new point
# and all previous points, and c is the variance of the new point.
# This is a more efficient way to calculate the inverse for all but the smallest
# covariance matrices (any more than ~70 points).
part_inv <- function(ainv, b, c) {
  invfact <- c - mahalanobis(t(b), center = FALSE, cov = ainv, inverted = TRUE)
  multi <- c(ainv %*% b)
  elem1 <- ainv + outer(multi, multi, "*")/invfact
  elemsym <- -multi/invfact
  out_mat <- matrix(0, nrow = nrow(ainv)+1, ncol = ncol(ainv)+1)
  out_mat[1:nrow(ainv), 1:ncol(ainv)] <- elem1
  out_mat[nrow(ainv)+1, 1:ncol(ainv)] <- t(elemsym)
  out_mat[1:nrow(ainv), ncol(ainv)+1] <- elemsym
  out_mat[nrow(ainv)+1, ncol(ainv)+1] <- 1/invfact
  return(out_mat)
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
    names <- unlist(map(outs, ~paste0(., times, sep = "")))
    relev <- setNames(relev, names)
    collected[[i]] <- relev
  }
  input_dat <- setNames(data.frame(matrix(rep(params, nreps), ncol = length(params), byrow = TRUE)), names(params))
  return(cbind(input_dat, do.call('rbind', collected)))
}

### PUBLIC FUNCTIONS
#' Create Boundary Emulators
#' 
#' Takes model run data and boundary information, and generates boundary emulators
#' 
#' For a stochastic model which admits (analytic) known boundaries, one can use
#' this information to produce more accurate emulators across the entire space
#' without requiring additional training points. This function takes existing
#' training points and information about the boundary and creates the corresponding
#' (Bayes linear) emulators for a stochastic known-boundary (KBE) system.
#' 
#' @param data_raw The model runs, where each row corresponds to one realisation
#' of the model at a given parameter combination
#' @param out_name The name of the output to emulate
#' @param ranges The parameter ranges, as a list of (lower upper) pairs
#' @param reps The number of repetitions afforded to each parameter combination
#' @param analytics A collection of analytic functions specific to the model, arising
#' as a result of the known boundary (see examples in `modelFunctions.R`)
#' @param bb_data A collection of functions for boundary-bulk prediction (provided
#' in `modelFunctions.R` for this particular case)
#' @param model The specifics used for the Gillespie algorithm to perform model runs
#' @param vals The boundary locations for each parameter
#' @param indices The corresponding parameter indices to which `vals` apply
#' @param out_index The index of the output in the list of model outputs
#' @param t The model time of evaluation for the emulated output
#' @param thetas The correlation lengths for the emulators; if NULL, they're
#' determined via bounded MAP estimation.
#' 
#' @returns A appropriate list of emulators for no-boundary, boundary, boundary + bulk.
create_boundary_ems <- function(data_raw, out_name, ranges, reps,
                                analytics, bb_data, model,
                                vals = c(0), indices = c(1), out_index, t,
                                thetas = c(1,1)) {
  get_summary <- function(data, input_names, out_name) {
    data_uids <- apply(data[,input_names], 1, rlang::hash)
    unique_uids <- unique(data_uids)
    out_arr <- array(0, dim = c(length(unique_uids), length(input_names)+2))
    for (i in seq_along(unique_uids)) {
      which_dat <- data[data_uids == unique_uids[i],]
      which_mean <- mean(which_dat[,out_name])
      which_var <- var(which_dat[,out_name])
      out_arr[i,] <- unlist(c(which_dat[1,input_names], which_mean, which_var), use.names = FALSE)
    }
    return(data.frame(out_arr) |> setNames(c(input_names, "exp", "var")))
  }
  data <- get_summary(data_raw, names(ranges), out_name)
  if (length(reps) == 1) reps <- rep(reps, nrow(data))
  if (!is.null(thetas))
    no_bound_ems <- hmer::emulator_from_data(data_raw, out_name, ranges,
                                             emulator_type = "variance",
                                             order = 1, beta.var = FALSE,
                                             specified_priors = list(
                                               variance = list(hyper_p = list(theta = thetas[1]), delta = 0),
                                               expectation = list(hyper_p = list(theta = thetas[2]), delta = 0)
                                             )
    )
  else
    no_bound_ems <- hmer::emulator_from_data(data_raw, out_name, ranges,
                                             emulator_type = "variance",
                                             order = 1, beta.var = FALSE,
                                             specified_priors = list(delta = 0)
    )
  prior_var_em <- no_bound_ems$variance[[out_name]]$o_em
  prior_exp_em <- no_bound_ems$expectation[[out_name]]$o_em
  boundary_em <- hmer::Proto_emulator$new(
    ranges, out_name, analytics$b_exp, analytics$b_cov,
    em = prior_var_em, analytic = function(y) analytics$analytic_sd(y, t, out_index, init_vals = model$M),
    vals = vals, indices = indices
  )
  gillesp_var <- analytics$b_cov(data[,names(ranges)], em = prior_var_em, full = TRUE, vals = vals, indices = indices) +
    map_dbl(seq_len(nrow(data)), ~prior_var_em$s_diag(data[.,], reps[.]))
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
    map_dbl(seq_len(nrow(data)), ~boundary_em$get_exp(data[.,])/reps[.])
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

#' Rep Allocation
#' 
#' Given a design, determines the 'optimal' means of allocating reps.
#' 
#' For a given collection of points and a budget of repetitions to be performed,
#' the optimal strategy might not to be to equally spread the repetitions across
#' the points. For instance, if a model is more stochastic in particular places,
#' it may be beneficial to allocate more repetitions to those parts of space at
#' the expense of the less variable parts. This function allocates the budget of
#' repetitions based on the principle of balancing the stochastic variability,
#' E[f_v(x)]/n, over all design points (f_v(x) here is the emulator for the
#' stochasticity).
#' 
#' Note that for an augmented design, the repetitions allocated can be afforded
#' to new design points or as an addition to old design points. The points provided
#' to this function therefore comprise both the new *and* old design points.
#' 
#' @param points The design points, including original design points and initial rep number
#' @param var_em The emulator for the stochasticity
#' @param rep_max The maximum number of repetitions to be allocated
#' @param ntoadd The number of repetitions to be allocated at each iteration
#' @param new_index The row-index of the first new design point.
#' 
#' @returns A data.frame of design points with a column corresponding to the reps.
rep_allocate <- function(points, var_em, rep_max, ntoadd, new_index = floor(nrow(points)/2)+1) {
  v_em_vals <- var_em$get_exp(points)
  rep_vals <- points$reps
  rep_vals[seq(new_index, length(rep_vals))] <- max(ntoadd,2)
  total_reps <- sum(rep_vals)
  while(total_reps < rep_max) {
    old_v_vals <- v_em_vals/rep_vals
    new_v_vals <- v_em_vals/(rep_vals + ntoadd)
    rep_add_index <- which.max(old_v_vals - new_v_vals)
    rep_vals[rep_add_index] <- rep_vals[rep_add_index] + ntoadd
    total_reps = total_reps + ntoadd
  }
  points$reps <- rep_vals
  return(points)
}

#' IMSPE Function
#'
#' Calculates 'integrated mean-squared prediction error', aka mean emulator
#' variance, generated due to the inclusion of a new training point. The mean
#' emulator variance computed is an approximation of the truth, based on evaluating
#' the emulator variance over a large grid of points. The results account for any
#' boundary information.
#' 
#' It can also be used to determine emulator variance across the space (for plotting,
#' for example) via the `return.raw` argument.
#' 
#' @param pt The grid on which to evaluate points
#' @param pre_em The base emulator, unexposed to training points or the boundary
#' @param data The training points currently used
#' @param new_point The new point to be included in the training data
#' @boundary_col The column (or columns) of the data frame for which a boundary exists
#' @boundary_val The corresponding locations of the boundaries; one per boundary_col element
#' @param return.raw If TRUE, returns the emulator variance at each (grid) point
#' 
#' @returns Either the mean emulator variance, or a data.frame of points and emulator variances.
imspe <- function(pt, pre_em, v_em, data,
                  new_point, invmat, rep,
                  boundary_col = 1, boundary_val = 0,
                  return.raw = FALSE) {
  n_data <- rbind.data.frame(data, new_point)
  data_mutate <- (n_data |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
  pt_mutate <- (pt |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
  term1 <- pre_em$u_sigma^2 * (1-r1(pt, pre_em, boundary_val, boundary_col)^2)
  term2a <- R1(pt, n_data, pre_em, boundary_val, boundary_col) * 
    pre_em$get_cov(pt_mutate, data_mutate, full = TRUE)
  b <- R1(data, new_point, pre_em, boundary_val, boundary_col) *
    pre_em$get_cov(data_mutate[-nrow(data_mutate),], data_mutate[nrow(data_mutate),,drop=FALSE])
  c <- as.numeric(R1(new_point, new_point, pre_em, boundary_val, boundary_col) * 
                    pre_em$get_cov(data_mutate[nrow(data_mutate),,drop=FALSE])) +
    as.numeric(v_em$get_exp(new_point)/rep)
  term2b <- part_inv(invmat, b, c)
  diag_res <- mahalanobis(term2a, center = FALSE, cov = term2b, inverted = TRUE)
  complete <- term1 - diag_res
  if (return.raw) {
    return(cbind.data.frame(pt, complete) |> setNames(c(names(pt), "V")))
  }
  complete <- complete[complete > 0]
  return(mean(complete))
}

#' Design choice for stochastic known-boundary models
#' 
#' Chooses a new collection of point locations for an emulated stochastic system.
#' 
#' Given an existing design and emulators trained on those design points, with
#' understanding of known boundaries in the system, this function aims to pick a
#' good design of points for future emulators to be trained on. The design process
#' assumes that 'infinite' repetitions are available when considering design sites;
#' the choice of where to put a (finite) repetition allocation is dealt with in
#' `rep_allocate()`.
#' 
#' @param data The training points provided to the emulators
#' @param b_em The 'basic' emulator: no boundary knowledge and not trained on `data`
#' @param v_em The (trained) variance emulator for the system
#' @param reps The  initial reps provided for each training point
#' @param ranges The parameter ranges, as a list of pairs (lower, upper)
#' @param testgrid_pts The number of points to evaluate on (for `imspe`)
#' @param pt_max The maximum number of new design points to select
#' @param boundary_col The parameters for which a boundary exists
#' @param boundary_val The corresponding locations of each boundary; one per parameter
#' @param nrepsadd The number of repetitions to add (initially) to each design point
#' @param verbose If TRUE, prints out details of the imspe at each selection
#' @param return_scores If TRUE, provides the imspe scores as well as the points
#' @param in_par If TRUE, tries to parallelise the computation
#' @param nrandomrestart Determines the number of different seeding locations for optim
#' 
#' @returns Either a data.frame of proposed points, or (if `return_scores` is TRUE)
#' a list consisting of this data.frame and a vector of imspe scores.
point_design <- function(data, b_em, v_em, reps, ranges, testgrid_pts,
                         pt_max = 2*nrow(data), rep_max = 2*sum(reps),
                         boundary_col, boundary_val, nrepsadd = reps[1],
                         verbose = FALSE, return_scores = FALSE, in_par = FALSE,
                         nrandomrestart = 10) {
  if (return_scores)
    p_scores <- c()
  tlhs <- lhs::randomLHS(testgrid_pts, length(ranges))
  testgrid <- data.frame(t(apply(tlhs, 1, function(x) {
    x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
  }))) |> setNames(names(ranges))
  grid_var <- v_em$get_exp(testgrid)/(floor(rep_max/pt_max))
  grid_var[grid_var < 0] <- 1e-6
  ave_variance <- mean(grid_var)
  # Possibly, at this point, allocate reps to existing points (if required)
  # to get them to similar variance
  if (verbose) print("Pre-allocation of reps to existing design points:")
  t_rep <- sum(reps)
  reps <- pmax(reps, round(v_em$get_exp(data)/ave_variance))
  if (verbose) print(paste("An additional", sum(reps)-t_rep, "repetitions added to existing design points."))
  find_next_point <- function(data, b_em, v_em, reps, ranges) {
    datamutate <- (data |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
    v_em_vals <- v_em$get_exp(data)
    v_em_vals[v_em_vals < 0] <- 1e-6
    start_mat <- R1(data, data, b_em, boundary_val, boundary_col) * b_em$get_cov(datamutate, full = TRUE) +
      diag(c(v_em_vals/reps))
    start_inv <- tryCatch(chol2inv(chol(start_mat)), error = function(e) MASS::ginv(start_mat))
    opt_func <- function(x) {
      x_mod <- data.frame(matrix(x, nrow = 1)) |> setNames(names(ranges))
      x_mod <- x_mod[,names(ranges), drop = FALSE]
      which_n <- max(2, round(v_em$get_exp(x_mod)/ave_variance))
      imspe(testgrid, b_em, v_em, data, x_mod, start_inv, which_n, boundary_col, boundary_val)
    }
    if (has_par_optim && in_par) {
      possible_points <- purrr::map(seq_len(nrandomrestart), function(i) {
        optimised <- optimParallel(map_dbl(ranges, ~runif(1, .[[1]], .[[2]])), opt_func,
                                   lower = map_dbl(ranges, ~.[[1]]),
                                   upper = map_dbl(ranges, ~.[[2]]),
                                   control = list(trace = FALSE))
        list(val = optimised$value, pt = data.frame(matrix(optimised$par, nrow = 1)) |> setNames(names(ranges)))
      })
      vals <- purrr::map_dbl(possible_points, "val")
      this_val <- min(vals)
      this_pt <- possible_points[[which.min(vals)]]$pt
      this_rep <- max(2, round(v_em$get_exp(this_pt)/ave_variance))
    }
    else {
      possible_points <- purrr::map(seq_len(nrandomrestart), function(i) {
        optimised <- optim(map_dbl(ranges, ~runif(1, .[[1]], .[[2]])), opt_func,
                           lower = map_dbl(ranges, ~.[[1]]),
                           upper = map_dbl(ranges, ~.[[2]]),
                           method = "L-BFGS-B", control = list(trace = FALSE))
        list(val = optimised$value, pt = data.frame(matrix(optimised$par, nrow = 1)) |> setNames(names(ranges)))
      })
      vals <- purrr::map_dbl(possible_points, "val")
      this_val <- min(vals)
      this_pt <- possible_points[[which.min(vals)]]$pt
      this_rep <- max(2, round(v_em$get_exp(this_pt)/ave_variance))
    }
    dists <- apply(data, 1, function(x) {
      sum((x-this_pt)^2)
    })
    if (this_val < 0) this_val <- Inf
    if (any(dists < 1e-6)) this_val <- NaN
    return(list(val = this_val, point = this_pt, rep = this_rep))
  }
  counter <- 0
  while(nrow(data) < pt_max) {
  # while(nrow(data) < pt_max && sum(reps) < rep_max) {
    pt_suggest <- find_next_point(data, b_em, v_em, reps, ranges)
    if (is.nan(pt_suggest$val)) next
    if (verbose) {
      print_str <- paste0("Proposal ", counter+1, ": Point score ", signif(pt_suggest$val, 4))
      print(print_str)
    }
    if (return_scores) {
      p_scores <- c(p_scores, pt_suggest$val)
    }
    data <- rbind.data.frame(data, pt_suggest$point)
    reps <- c(reps, pt_suggest$rep)
    counter <- counter + 1
  }
  # if (sum(reps) >= rep_max && nrow(data) < pt_max && verbose) {
  #   print("Repetition budget expended before new point budget.")
  # }
  pts_with_reps <- cbind.data.frame(data, reps) |> setNames(c(names(data), "reps"))
  if (return_scores)
    return(list(points = pts_with_reps, pt_scores = p_scores))
  return(pts_with_reps)
}

########################
#### NOW DEPRECATED ####
########################
## Scoring functions for variance reduction
# Calculate the mean emulator variance across the space, using a representative
# collection of points.
# mean_em_var <- function(pt, pre_em, var_em, data, reps, boundary_col = 1, boundary_val = 0, 
#                         invmat, new_point = NULL, new_rep = NULL, new_method = (nrow(data) > 70),
#                         return.raw = FALSE) {
#   if (!is.null(new_point)) {
#     n_data <- rbind.data.frame(data, new_point)
#     datamutate <- (n_data |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
#   }
#   else {
#     datamutate <- (data |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
#   }
#   ptmutate <- (pt |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
#   term1 <- pre_em$u_sigma^2 * (1-r1(pt, pre_em, boundary_val, boundary_col)^2)
#   if (!is.null(new_point) && new_method) {
#     npmutate <- (new_point |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
#     term2a <- R1(pt, n_data, pre_em, boundary_val, boundary_col) * pre_em$get_cov(ptmutate, datamutate, full = TRUE)
#     var_em_val_add <- var_em$get_exp(new_point)
#     if (var_em_val_add < 0) var_em_val_add <- 1e-6
#     b <- R1(data, new_point, pre_em, boundary_val, boundary_col) * pre_em$get_cov(datamutate[-nrow(datamutate),], npmutate)
#     c <- as.numeric(R1(new_point, new_point, pre_em, boundary_val, boundary_col) * pre_em$get_cov(npmutate) +
#       var_em_val_add/reps[length(reps)])
#     term2b <- part_inv(invmat, b, c)
#   }
#   else if (!is.null(new_rep) && new_method) {
#     rep_ind <- which(new_rep != 0)
#     modif <- new_rep[rep_ind]/(reps[rep_ind]*(reps+new_rep)[rep_ind])
#     term2a <- R1(pt, data, pre_em, boundary_val, boundary_col) * pre_em$get_cov(ptmutate, datamutate, full = TRUE)
#     var_em_val <- var_em$get_exp(data[rep_ind,,drop=FALSE])
#     if (var_em_val < 0) var_em_val <- 1e-6
#     uval <- as.numeric(sqrt(modif * var_em_val))
#     denom <- 1 - uval^2 * invmat[rep_ind,rep_ind]
#     num <- uval^2 * outer(c(invmat[rep_ind,]), c(invmat[rep_ind,]), "*")
#     term2b <- invmat + num/denom
#   }
#   else {
#     if (!is.null(new_rep)) reps <- reps + new_rep
#     if (!is.null(new_point)) data <- rbind.data.frame(data, new_point)
#     term2a <- R1(pt, data, pre_em, boundary_val, boundary_col) * pre_em$get_cov(ptmutate, datamutate, full = TRUE)
#     var_em_vals <- var_em$get_exp(data)
#     var_em_vals[var_em_vals < 0] <- 1e-6
#     term2binv <- R1(data, data, pre_em, boundary_val, boundary_col) * pre_em$get_cov(datamutate, full = TRUE) +
#       diag(c(1/reps * var_em_vals))
#     term2b <- tryCatch(chol2inv(chol(term2binv)),
#                        error = function(e) MASS::ginv(term2binv))
#   }
#   diag_res <- mahalanobis(term2a, center = FALSE, cov = term2b, inverted = TRUE)
#   complete <- term1 - diag_res
#   if (return.raw) {
#     return(cbind.data.frame(pt, complete) |> setNames(c(names(pt), "V")))
#   }
#   complete <- complete[complete > 0]
#   #complete[complete < 0] <- 1e-6
#   return(mean(complete))
# }
# Calculate the score due to including a new design point
# new_point_score <- function(point, points, pre_em, var_em, reps, grid, invmat, ntoadd = 1,
#                             boundary_col, boundary_val) {
#   #new_data <- rbind.data.frame(points, point)
#   if (is.null(point)) {
#     warning("Point is NULL.")
#     new_reps <- reps
#   }
#   else
#     new_reps <- c(reps, ntoadd)
#   vars <- mean_em_var(grid, pre_em, var_em, points, new_reps, boundary_col, boundary_val, invmat, point)
#   return(vars)
# }
# Calculate the score due to adding a repetition at an existing design point
# new_rep_score <- function(index, points, pre_em, var_em, reps, grid, invmat, ntoadd = 1,
#                           boundary_col, boundary_val) {
#   new_reps <- rep(0, length(reps))
#   new_reps[index] <- ntoadd
#   vars <- mean_em_var(grid, pre_em, var_em, points, reps, boundary_col,
#                       boundary_val, invmat, new_rep = new_reps)
#   return(vars)
# }
## Chooses an 'optimal' design
# Given a set of data points, progressively adds points to a candidate
# set or adds a rep to an existing point in the candidate set. Calculations of
# expected improvement are based on evaluating the mean emulator variance
# over a grid of points, testgrid: the more dense the grid, the more accurate the
# estimate, but the more computationally intensive it is. Returns the collection
# of points, along with a column denoting how many reps are to be run at each point.
# design_subselect <- function(data, pre_em, var_em, reps, ranges, testgrid,
#                               rep_max, pt_max, ntoadd = 1,
#                               store_order = FALSE, verbose = FALSE,
#                               return_scores = FALSE, boundary_col = 1, boundary_val = 0,
#                              rep_favour_factor = 1) {
#   if (store_order) order_vector <- c()
#   if (return_scores) {
#     r_scores <- c()
#     p_scores <- c()
#   }
#   find_next_point <- function(data, pre_em, var_em, reps, ranges, ntoadd, boundary_col, boundary_val) {
#     datamutate <- (data |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
#     v_em_vals <- var_em$get_exp(data)
#     v_em_vals[v_em_vals < 0] <- 1e-6
#     start_inv <- R1(data, data, pre_em, boundary_val, boundary_col) * pre_em$get_cov(datamutate, full = TRUE) +
#       diag(c(v_em_vals/reps))
#     start <- tryCatch(
#       chol2inv(chol(start_inv)),
#       error = function(e) MASS::ginv(start_inv)
#     )
#     opt_func <- function(x) {
#       x_mod <- data.frame(matrix(x, nrow = 1)) |> setNames(names(ranges))
#       x_mod <- x_mod[,names(ranges), drop = FALSE]
#       new_point_score(x_mod, data, pre_em, var_em, reps, testgrid, start, ntoadd,
#                       boundary_col, boundary_val)
#     }
#     if (has_par_optim)
#       optimised <- optimParallel(map_dbl(ranges, ~runif(1, .[[1]], .[[2]])), opt_func,
#                          lower = map_dbl(ranges, ~.[[1]]+0.01*diff(.)),
#                          upper = map_dbl(ranges, ~.[[2]]-0.01*diff(.)),
#                          control = list(trace = FALSE))
#     else
#       optimised <- optim(map_dbl(ranges, ~runif(1, .[[1]], .[[2]])), opt_func,
#                          lower = map_dbl(ranges, ~.[[1]]+0.01*diff(.)),
#                          upper = map_dbl(ranges, ~.[[2]]-0.01*diff(.)),
#                          method = "L-BFGS-B", control = list(trace = FALSE))
#     this_val <- optimised$value
#     this_pt <- data.frame(matrix(optimised$par, nrow = 1)) |> setNames(names(ranges))
#     dists <- apply(data, 1, function(x) {
#       sum((x-this_pt)^2)
#     })
#     if (this_val < 0) this_val <- Inf
#     if (any(dists < 1e-6)) this_val <- NaN
#     return(list(val = this_val, point = this_pt))
#   }
#   find_next_rep <- function(data, pre_em, var_em, reps, ntoadd, boundary_col, boundary_val) {
#     datamutate <- (data |> dplyr::mutate(across(all_of(boundary_col), ~boundary_val)))
#     v_em_vals <- var_em$get_exp(data)
#     v_em_vals[v_em_vals < 0] <- 1e-6
#     start_inv <- R1(data, data, pre_em, boundary_val, boundary_col) * pre_em$get_cov(datamutate, full = TRUE) +
#       diag(c(1/v_em_vals * reps))
#     start <- tryCatch(
#       chol2inv(chol(start_inv)),
#       error = function(e) MASS::ginv(start_inv)
#     )
#     if (has_future)
#       rep_vals <- furrr::future_map_dbl(seq_len(nrow(data)), function(i) {
#         new_rep_score(i, data, pre_em, var_em, reps, testgrid, start, ntoadd, boundary_col, boundary_val)
#       })
#     else
#       rep_vals <- map_dbl(seq_len(nrow(data)), function(i) {
#         new_rep_score(i, data, pre_em, var_em, reps, testgrid, start, ntoadd, boundary_col, boundary_val)
#       })
#     return(list(val = min(rep_vals), index = which.min(rep_vals)))
#   }
#   failsafe <- 0
#   while(sum(reps) < rep_max && failsafe < 10000) {
#     rep_suggest <- find_next_rep(data, pre_em, var_em, reps, ntoadd, boundary_col, boundary_val)
#     if (nrow(data) < pt_max) {
#       pt_suggest <- find_next_point(data, pre_em, var_em, reps, ranges, ntoadd, boundary_col, boundary_val)
#       pt_suggest$val <- rep_favour_factor * pt_suggest$val
#     }
#     else
#       pt_suggest <- NULL
#     if (verbose) {
#       print_str <- paste0("Proposal ", failsafe+1, ":")
#       if (!is.null(pt_suggest)) {
#         print_str <- paste0(print_str, " Point score ", signif(pt_suggest$val, 4), ";")
#       }
#       print_str <- paste0(print_str, " Rep score ", signif(rep_suggest$val, 4))
#       if (!is.null(pt_suggest) && !is.nan(pt_suggest$val)) {
#         if (pt_suggest$val < rep_suggest$val)
#           print_str <- paste0(print_str, " - New point chosen;")
#         else
#           print_str <- paste0(print_str, " - Extra rep chosen;")
#       }
#       else
#         print_str <- paste0(print_str, " - Extra rep chosen;")
#       print(print_str)
#     }
#     if (return_scores) {
#       r_scores <- c(r_scores, rep_suggest$val)
#       if (is.null(pt_suggest)) p_scores <- c(p_scores, Inf)
#       else p_scores <- c(p_scores, pt_suggest$val)
#     }
#     if (is.null(pt_suggest) || is.nan(pt_suggest$val) || rep_suggest$val < pt_suggest$val) {
#       reps[rep_suggest$index] <- reps[rep_suggest$index] + ntoadd
#       if (store_order) order_vector <- c(order_vector, rep_suggest$index)
#     }
#     else {
#       data <- rbind.data.frame(data, pt_suggest$point)
#       reps <- c(reps, ntoadd)
#       if (store_order) order_vector <- c(order_vector, nrow(data))
#     }
#     failsafe <- failsafe + 1
#   }
#   pts_with_reps <- cbind.data.frame(data, reps) |> setNames(c(names(data), "reps"))
#   if (store_order) {
#     if (return_scores)
#       return(list(points = pts_with_reps, order = order_vector, rep_scores = r_scores, pt_scores = p_scores))
#     return(list(points = pts_with_reps, order = order_vector))
#   }
#   return(pts_with_reps)
# }
