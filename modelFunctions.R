############################
# Model-Specific Functions #
############################

source("baseFunctions.R")

## Shorthand function for r_1(x)
r1 <- function(pts, em, val = 0, index = 1) {
  if (!is.data.frame(pts)) pts <- data.frame(pts)
  pts <- pts |> dplyr::mutate(across(seq_len(length(pts))[-index], ~val))
  ref_pt <- data.frame(matrix(rep(0, length(pts)), nrow = 1)) |> setNames(names(pts))
  ref_pt[,index] <- val
  return(
    em$get_cov(pts, ref_pt, full = TRUE, check_neg = FALSE)/em$u_sigma^2
  )
}
## Shorthand function for R_1(x,x')=r_1(x-x')-r_1(x)r_1(x')
R1 <- function(p1, p2, em, val = 0, index = 1) {
  if(missing(p2)) p2 <- p1
  outer(seq_len(nrow(p1)), seq_len(nrow(p2)), function(i,j) {
    r1(p1[i,]-p2[j,], em, val, index) - r1(p1[i,], em, val, index)*r1(p2[j,], em, val, index)
  })
}

## SIR Model - Single Boundary
SIR_functions <- list(
  analytic_mean = function(params, t, index, init_vals) {
    params <- unlist(params, use.names = FALSE)
    gamma <- params[2]
    if (length(params) < 3)
      omega <- 0.24
    else
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
  },
  analytic_sd = function(params, t, index1, index2, init_vals) {
    params <- unlist(params, use.names = FALSE)
    gamma <- params[2]
    if (length(params) < 3)
      omega <- 0.24
    else
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
  },
  b_exp = function(x, em, analytic, vals = c(0), indices = c(1)) {
    xK <- (x |> dplyr::mutate(across(indices[1], ~vals[1])))
    return(em$get_exp(x) + r1(x, em, val = vals[1], index = indices[1]) * (apply(xK, 1, analytic) - em$get_exp(xK)))
  },
  b_cov = function(x, xp = NULL, full = TRUE, em, vals = c(0), indices = c(1)) {
    xK <- (x |> dplyr::mutate(across(indices[1], ~vals[1])))
    if (is.null(xp)) {
      xp <- x
      xpK <- xK
    }
    else xpK <- (xp |> dplyr::mutate(across(indices[1], ~vals[1])))
    if (full) return(
      R1(x, xp, em, vals[1], indices[1])*em$get_cov(xK, xpK, full = TRUE)
    )
    return(
      purrr::map_dbl(seq_len(nrow(x)), ~R1(x[.,], xp[.,], em, vals[1], indices[1]))*em$get_cov(xK, xpK)
    )
  }
)

SEIR_functions <- list(
  analytic_mean = function(params, t, index, init_vals) {
    if (index != 3) stop("No analytic solution calculated for outputs other than I.")
    I0 <- init_vals[3]
    E0 <- init_vals[2]
    muep_exp <- exp(-(params[2]+params[4])*t)
    lim_comb <- params[5]+params[6]-params[4]
    if (lim_comb == 0) return((params[4]*t*E0+I0)*muep_exp)
    mualga_exp <- exp(-(params[2]+params[5]+params[6])*t)
    algaep_exp <- exp(-lim_comb*t)
    return(params[4]*E0/lim_comb * muep_exp * (1 - algaep_exp) + I0 * mualga_exp)
  },
  analytic_sd = function(params, t, index1, index2, init_vals) {
    if (missing(index1)) index1 <- 3
    if (missing(index2)) index2 <- 3
    if (index1 != 3 || index2 != 3) stop("No analytic solution calculated for outputs other than I.")
    I0 <- init_vals[3]
    E0 <- init_vals[2]
    muep_exp <- exp(-(params[2]+params[4])*t)
    lim_comb <- params[5]+params[6]-params[4]
    if (lim_comb == 0) return(I0 * muep_exp * (1 - muep_exp) + E0 * params[4] * t * muep_exp * (1 - params[4] * t * muep_exp))
    algaep_exp <- exp(-lim_comb*t)
    mualga_exp <- exp(-(params[2]+params[5]+params[6])*t)
    comp1 <- params[4]/lim_comb * muep_exp * (1 - algaep_exp)
    return(E0 * comp1 * (1-comp1) + I0 * mualga_exp * (1 - mualga_exp))
  },
  b_exp = function(x, em, analytic, vals = c(0,0), indices = c(3,4)) {
    xK <- x |> dplyr::mutate(across(indices[1], ~vals[1]))
    xL <- x |> dplyr::mutate(across(indices[2], ~vals[2]))
    xKL <- xK |> dplyr::mutate(across(indices[2], ~vals[2]))
    return(
      em$get_exp(x, check_neg = FALSE) +
        r1(x, em, vals[1], indices[1]) * (apply(xK, 1, analytic) - em$get_exp(xK, check_neg = FALSE)) +
        r1(x, em, vals[2], indices[2]) * (apply(xL, 1, analytic) - em$get_exp(xL, check_neg = FALSE)) -
        r1(x, em, vals[1], indices[1]) * r1(x, em, vals[2], indices[2]) * (apply(xKL, 1, analytic) - em$get_exp(xKL, check_neg = FALSE))
    )
  },
  b_cov = function(x, xp = NULL, full = TRUE, em, vals = c(0,0), indices = c(3,4)) {
    xKL <- x |> dplyr::mutate(across(indices[1], ~vals[1])) |> dplyr::mutate(across(indices[2], ~vals[2]))
    if (is.null(xp)) {
      xp <- x
      xpKL <- xKL
    }
    else xpKL <- xp |> dplyr::mutate(across(indices[1], ~vals[1])) |> mutate(across(indices[2], ~vals[2]))
    if (full) return(
      R1(x, xp, em, vals[1], indices[1]) * R1(x, xp, em, vals[2], indices[2]) * em$get_cov(xKL, xpKL, full = TRUE, check_neg = FALSE)
    )
    return(
      purrr::map_dbl(seq_len(nrow(x)), ~R1(x[.,], xp[.,], em, vals[1], indices[1]) * R1(x[.,], xp[.,], em, vals[2], indices[2]))*em$get_cov(xKL, xpKL, check_neg = FALSE)
    )
  }
)


# Function for obtaining the boundary-bulk adjusted emulator expectation
bb_exp <- function(x, b_em, pre_em, dat, binv, bmod, bcov, analytic, vals = c(0), indices = c(1)) {
  if (nrow(x) > 1000) {
    output_vals <- c()
    for (i in 1:ceiling(nrow(x)/1000)) {
      px <- x[(1000*(i-1)+1):min(1000*i, nrow(x)),]
      output_vals <- c(output_vals, bb_exp(px, b_em, pre_em, dat, binv, bmod, analytic, vals, indices))
    }
    return(output_vals)
  }
  b_em$get_exp(x) + bcov(x, dat, pre_em, full = TRUE, vals = vals, indices = indices) %*%
    binv %*% bmod
}

## Function to obtain the boundary-bulk adjusted emulator (co)variance
bb_cov <- function(x, xp = NULL, full = TRUE, b_em, pre_em, dat, binv, bcov, analytic, vals = c(0), indices = c(1)) {
  if (nrow(x) > 1000 && is.null(xp)) {
    output_vals <- c()
    for (i in 1:ceiling(nrow(x)/1000)) {
      px <- x[(1000*(i-1)+1):min(1000*i, nrow(x)),]
      output_vals <- c(output_vals, bb_cov(px, xp, full, b_em, pre_em, dat, binv, analytic, vals, indices))
    }
    return(output_vals)
  }
  cov_mat <- bcov(x, dat, pre_em, full = TRUE, vals = vals, indices = indices)
  if (is.null(xp))
    cov_mat_p <- cov_mat
  else
    cov_mat_p <- bcov(xp, dat, pre_em, full = TRUE, vals = vals, indices = indices)
  if (full)
    return(b_em$get_cov(x, xp, full = TRUE) - cov_mat %*% binv %*% t(cov_mat_p))
  return(b_em$get_cov(x) - diag(cov_mat %*% binv %*% t(cov_mat_p)))
}

## Function for the implausibility of the boundary-bulk emulator
bb_imp <- function(x, z, cutoff = NULL, b_em, pre_em, v_em, dat, binv, bmod, bcov, analytic, vals = c(0), indices = c(1)) {
  bb_exp_diff <- abs(bb_exp(x, b_em, pre_em, dat, binv, bmod, bcov, analytic, vals, indices) - z$val)
  bb_norm <- bb_cov(x, xp = NULL, full = FALSE, b_em, pre_em, dat, binv, bcov, analytic, vals, indices) + z$sigma^2 +
    v_em$get_exp(x)/reps
  imps <- bb_exp_diff/sqrt(bb_norm)
  if (is.null(cutoff)) return(imps)
  return(imps <= cutoff)
}

bb_data <- list(
  bb_exp = bb_exp,
  bb_cov = bb_cov,
  bb_imp = bb_imp
)
