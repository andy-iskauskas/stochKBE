############################
# Model-Specific Functions #
############################

source("baseFunctions.R")

## Analytic Mean
# Gives the analytic mean for the given SIR model,
# according to the exact solutions.
analytic_mean <- function(params, t, index, init_vals = N$M) {
  params <- unlist(params, use.names = FALSE)
  gamma <- params[2]
  omega <- params[3]
  S0 <- init_vals[1]
  I0 <- init_vals[2]
  R0 <- init_vals[3]
  exp_g <- exp(-gamma * t)
  exp_o <- exp(-omega * t)
  frac <- ifelse(gamma == omega, 0, 1/(omega - gamma))
  out_vec <- c(
    I0 * frac * (gamma * exp_o - omega * exp_g) +
      R0 * (1 - exp_o) + S0 + I0,
    I0 * exp_g,
    R0 * exp_o + gamma * frac * I0 * (exp_g - exp_o)
  )
  if (missing(index)) return(out_vec)
  return(out_vec[index])
} 

## Analytic SD
# Gives the analytic standard deviation for the SIR model,
# according to the exact solutions.
analytic_sd <- function(params, t, index1, index2, init_vals = N$M) {
  params <- unlist(params, use.names = FALSE)
  gamma <- params[2]
  omega <- params[3]
  S0 <- init_vals[1]
  I0 <- init_vals[2]
  R0 <- init_vals[3]
  exp_g <- exp(-gamma * t)
  exp_o <- exp(-omega * t)
  frac <- ifelse(gamma == omega, 0, 1/(omega - gamma))
  lambda <- gamma * frac
  vS <- -I0*frac^2 * (gamma * exp_o - omega * exp_g) *
    (omega * (1 - exp_g) - gamma * (1 - exp_o)) + R0 * exp_o * (1 - exp_o)
  vI <- I0 * exp_g * (1 - exp_g)
  vR <- -lambda^2 * I0 * (exp_g - exp_o)^2 + lambda * I0 * (exp_g - exp_o) +
    R0 * exp_o * (1 - exp_o)
  cSI <- -I0 * frac * exp_g * (omega * (1 - exp_g) - gamma * (1 - exp_o))
  cIR <- -lambda * I0 * exp_g * (exp_g - exp_o)
  cRS <- -lambda * I0 * (exp_g - exp_o) * (omega * (1 - exp_g) - gamma * (1 - exp_o)) -
    R0 * exp_o * (1 - exp_o)
  out_mat <- matrix(c(
    vS, cSI, cRS, cSI, vI, cIR, cRS, cIR, vR
  ), nrow = 3)
  if (missing(index1) && missing(index2)) return(out_mat)
  if (missing(index2)) return(out_mat[index1, index1])
  if (missing(index1)) return(out_mat[index2, index2])
  return(out_mat[index1, index2])
}

## Shorthand function for r_1(x)
r1 <- function(pts, em, val = 0) {
  if (!is.data.frame(pts)) pts <- data.frame(pts)
  pts <- pts |> dplyr::mutate(across(seq_len(length(pts)-1)+1, ~0))
  ref_pt <- data.frame(matrix(c(val, rep(0, length(pts)-1)), nrow = 1)) |>
    setNames(names(pts))
  return(
    em$get_cov(pts, ref_pt, full = TRUE)/em$u_sigma^2
  )
}

## Shorthand function for R_1(x,x')=r_1(x-x')-r_1(x)r_1(x')
R1 <- function(p1, p2, em, val = 0) {
  if(missing(p2)) p2 <- p1
  outer(seq_len(nrow(p1)), seq_len(nrow(p2)), function(i,j) {
    r1(p1[i,]-p2[j,], em, val) - r1(p1[i,], em, val)*r1(p2[j,], em, val)
  })
}

## Adjustment functions
# Analytic update for the expectation of the boundary emulator
b_exp <- function(x, em, analytic, val = 0) {
  xK <- (x |> dplyr::mutate(across(1, ~val)))
  return(em$get_exp(x) + r1(x, em, val = val) * (apply(xK, 1, analytic) - em$get_exp(xK)))
}

# Analytic update for the (co)variance of the boundary emulator
b_cov <- function(x, xp = NULL, full = TRUE, em, val = 0) {
  xK <- (x |> dplyr::mutate(across(1, ~val)))
  if (is.null(xp)) {
    xp <- x
    xpK <- xK
  }
  else xpK <- (xp |> dplyr::mutate(across(1, ~val)))
  if (full) return(
    R1(x, xp, em, val)*em$get_cov(xK, xpK, full = TRUE)
  )
  return(
    purrr::map_dbl(seq_len(nrow(x)), ~R1(x[.,], xp[.,], em, val))*em$get_cov(xK, xpK)
  )
}

# Function for obtaining the boundary-bulk adjusted emulator expectation
bb_exp <- function(x, b_em, pre_em, dat, binv, bmod, analytic, val = 0) {
  if (nrow(x) > 1000) {
    output_vals <- c()
    for (i in 1:ceiling(nrow(x)/1000)) {
      px <- x[(1000*(i-1)+1):min(1000*i, nrow(x)),]
      output_vals <- c(output_vals, bb_exp(px, b_em, pre_em, dat, binv, bmod, analytic, val))
    }
    return(output_vals)
  }
  b_em$get_exp(x) + b_cov(x, dat, pre_em, full = TRUE, val = val) %*%
    binv %*% bmod
}

## Function to obtain the boundary-bulk adjusted emulator (co)variance
bb_cov <- function(x, xp = NULL, full = TRUE, b_em, pre_em, dat, binv, analytic, val = 0) {
  if (nrow(x) > 1000 && is.null(xp)) {
    output_vals <- c()
    for (i in 1:ceiling(nrow(x)/1000)) {
      px <- x[(1000*(i-1)+1):min(1000*i, nrow(x)),]
      output_vals <- c(output_vals, bb_cov(px, xp, full, b_em, pre_em, dat, binv, analytic, val))
    }
    return(output_vals)
  }
  cov_mat <- b_cov(x, dat, pre_em, full = TRUE, val = val)
  if (is.null(xp))
    cov_mat_p <- cov_mat
  else
    cov_mat_p <- b_cov(xp, dat, pre_em, full = TRUE, val = val)
  if (full)
    return(b_em$get_cov(x, xp, full = TRUE) - cov_mat %*% binv %*% t(cov_mat_p))
  return(b_em$get_cov(x) - diag(cov_mat %*% binv %*% t(cov_mat_p)))
}

## Function for the implausibility of the boundary-bulk emulator
bb_imp <- function(x, z, cutoff = NULL, b_em, pre_em, v_em, dat, binv, bmod, analytic, val = 0) {
  bb_exp_diff <- abs(bb_exp(x, b_em, pre_em, dat, binv, bmod, analytic, val) - z$val)
  bb_norm <- bb_cov(x, xp = NULL, full = FALSE, b_em, pre_em, dat, binv, analytic, val) + z$sigma^2 +
    v_em$get_exp(x)/reps
  imps <- bb_exp_diff/sqrt(bb_norm)
  if (is.null(cutoff)) return(imps)
  return(imps <= cutoff)
}

