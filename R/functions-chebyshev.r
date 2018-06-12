
#' @title Scaling of X-Axis to [-1,1]
#' @param x_axis Arbitrary x-axis
#' @param x_val Single value in range of x_axis
#' @return Depending on input parameters: Either scaled x-axis or scaled scalar 
#'   value
#' @description
#' \code{cheb_scale} scales an arbitrary x-axis to the range [-1, 1], so it can 
#'   be used for a Chebyshev polynomial, which is defined on this span.
#' @export
cheb_scale <- function(x_axis, x_val = NA) {
  if ( is.na(x_val) == TRUE ) {
    x_cheb_scaled <- 2 * (x_axis - min(x_axis)) / (max(x_axis) - min(x_axis)) - 1
  }
  if ( is.na(x_val) == FALSE ) {
    x_cheb_scaled <- 2 * (x_val - min(x_axis)) / (max(x_axis) - min(x_axis)) - 1
  }
  return(x_cheb_scaled)
}


#' @title Rescaling of X-Axis to Original X-Axis
#' @param x_cheb Scaled x-axis
#' @param x_axis Unscaled x-axis
#' @return Rescaled x-axis
#' @description
#' \code{cheb_rescale} rescales (single or multiple) values in the range [-1,1] 
#' to the original x-axis.
#' @export
cheb_rescale <- function(x_cheb, x_axis) {
  x_rescaled <- (1/2 * (x_cheb + 1) * (max(x_axis) - min(x_axis))) + x_axis[1]
  return(x_rescaled)
}


#' @title Generating Chebyshev Polynomials of first kind
#' @param x_axis Arbitrary x-axis
#' @param n Order of the polynomial
#' @return Chebyshev polynomial of the first kind
#' @description
#' \code{cheb_1st} generates Chebyshev polynomials of the first kind of the 
#'   order n from an arbitrary x-axis. Recurrence formula is taken from 
#'   Bronstein p.xxx.
#' @export
cheb_1st <- function(x_axis, n){
  x_cheb <- cheb_scale(x_axis)
  m <- n + 1
  cheb_t_0 <- 1;  cheb_t_1 <- x_cheb;
  cheb_t <- cbind(cheb_t_0, cheb_t_1)
  if (n >= 2) {
    for (i in 3:m) {
      cheb_t_i <- 2 * x_cheb * cheb_t[,(i - 1)] - cheb_t[,(i - 2)]
      cheb_t <- cbind(cheb_t, cheb_t_i)
      rm(cheb_t_i)
    }
  }
  return(cheb_t)
}


#' @title Generating Chebyshev Polynomials of second kind
#' @param x_axis Arbitrary x-axis
#' @param n Order of the polynomial
#' @return Chebyshev polynomials of the second kind
#' @description
#' \code{cheb_2nd} generates Chebyshev polynomials of the second kind of the 
#'   order n from an arbitrary x-axis. Recurrence formula is taken from 
#'   Bronstein p.xxx.
#' @export
cheb_2nd <- function(x_axis, n){
  x_cheb <- cheb_scale(x_axis)
  m <- n + 1
  cheb_u_0 <- 1; cheb_u_1 <-  2*x_cheb
  cheb_u <- cbind(cheb_u_0, cheb_u_1)
  if (n >= 2) {
    for (i in 3:m) {
      cheb_u_i <- 2 * x_cheb * cheb_u[,(i - 1)] - cheb_u[,(i - 2)]
      cheb_u <- cbind(cheb_u, cheb_u_i)
      rm(cheb_u_i)
    }
  }
  return(cheb_u)
}


#' @title Calculation of Values of the model fit
#' @param x_axis Arbitrary x-axis or single value
#' @param cheb_coeff Coefficients of the Least-Squares-Fit with Chebyshev
#'   polynomials
#' @return Filtered model (or single value)
#' @description
#' \code{cheb_model} calculates filtered model values of arbitrary 
#' @export
cheb_model_filter <- function(x_axis, cheb_coeff) {
  n <- length(cheb_coeff) - 1
  cheb_t <- cheb_1st(x_axis, n)
  cheb_model <- cheb_t %*% cheb_coeff
  return(cheb_model)
}


#' @title Calculation of the values of the first derivation
#' @param x_axis Arbitrary x-axis
#' @param cheb_coeff Chebyshev-coefficients of the least squares model fit
#' @return Values of the first derivation of the filtered 
#'   model
#' @description
#' \code{cheb_deriv_1st} calculates values of the first derivation of the 
#'   filtered model from the Chebyshev coefficients of the least squares fit on 
#'   the original x-axis or single x-values.
#' @export
cheb_deriv_1st <- function(x_axis, cheb_coeff) {
  # if (length(x_axis) != 0) { ### ueberpruefen, ob noetig
  n <- length(cheb_coeff) - 1
  # m <- n + 1
  cheb_u <- cheb_2nd(x_axis, n)
  
  # recurrence formula for first derivation of first kind chebyshev polynomials
  # bronstein pp. 
  # dT/dx = n * U_(n-1)
  
  cheb_t_deriv_1st <- if (length(x_axis) == 1) (1:n)*t(cheb_u[,1:n]) else t((1:n) * t(cheb_u[,1:n]))
  cheb_t_deriv_1st <- cbind(0, cheb_t_deriv_1st)
  cheb_model_deriv_1st <- cheb_t_deriv_1st %*% cheb_coeff
  
  return(cheb_model_deriv_1st)
  # }
}


#' @title Calculation of the values of the second derivation
#' @param x_axis Arbitrary X-Axis
#' @param cheb_coeff Chebyshev-Coefficients
#' @return Second Derivation of the filtered model
#' @description
#' \code{cheb_deriv_2nd} calculates values of the second derivation of the 
#'   filtered model from the Chebyshev coefficients of the least squares fit on 
#'   the original x-axis or single x-values.
#' @export
cheb_deriv_2nd <- function(x_axis, cheb_coeff) {
  n <- length(cheb_coeff) - 1
  #  m <- n + 1
  l <- length(x_axis)
  cheb_t <- cheb_1st(x_axis, n)
  cheb_u <- cheb_2nd(x_axis, n)
  x_cheb <- cheb_scale(x_axis)

  # recurrence formula for first derivation of first kind chebyshev polynomials
  # bronstein pp. 
  # dT/dx = n * U_(n-1)
  
  cheb_t_deriv_2nd <- t((0:n) * t(t((0:n + 1) * t(cheb_t )) - cheb_u)) / (x_cheb ** 2 - 1)
  cheb_t_deriv_2nd[1,] <- (-1) ** (0:n) * ((0:n) ** 4 - (0:n) ** 2) / (3)
  cheb_t_deriv_2nd[l,] <- ((0:n) ** 4 - (0:n) ** 2) / (3)
  cheb_model_deriv_2nd <- cheb_t_deriv_2nd %*% cheb_coeff
  
  return(cheb_model_deriv_2nd)
}


#' @title Curve Fitting with Chebyshev Polynomials
#' @param d Data series to be fitted
#' @param x_axis Arbitrary X-Axis
#' @param n Order of the polynomial
#' @return Values of the filtered model
#' @description
#' \code{cheb_fit} fits a Chebyshev polynomial of order n to a arbitrary data 
#'   or time series using the method of least squares.
#' @export
cheb_fit <- function(d, x_axis, n){
  ## Find missing values
  mss_ind  <- which(is.na(d))

  # case discrimination for missing values, that lead to errors.
  if (length(mss_ind) >= 1) {
    cheb_t <- cheb_1st(x_axis[-mss_ind], n)
  } else {
    cheb_t <- cheb_1st(x_axis, n)
  }

  ## model calculations
  # chebyshev polynomial coefficients
  if (length(mss_ind) >= 1) {
    cheb_coeff <- solve(t(cheb_t) %*% cheb_t) %*% t(cheb_t) %*% d[-mss_ind]
  } else {
    cheb_coeff <- solve(t(cheb_t) %*% cheb_t) %*% t(cheb_t) %*% d
  }

  x_cheb <- cheb_scale(x_axis)
  cheb_model <- cheb_model_filter(x_cheb, cheb_coeff)

  return(cheb_model)
}


#' @title Curve Fitting with Chebyshev Polynomials and Finding of its Roots
#' @param d Data series to be fitted
#' @param x_axis Arbitrary x-axis
#' @param n Order of the polynomial
#' @param roots_bound_lower Lower boundary of the extreme value search
#' @param roots_bound_upper Upper boundary of the extreme value search
#' @return list of calculated parameters like Chebyshev coefficients,
#'   values of the filtered model, values of the first and seconda derivation,
#'   postitions and values of extremes
#' @description
#' \code{cheb_fit_extr} fits a Chebyshev polynomial of order n to a data series 
#'   using the method of least squares and locates the extrema of the fit using 
#'   the Newton-Raphson-algorithm.
#' @export
#' @importFrom rootSolve uniroot.all
cheb_find_extr <- function(d, x_axis, n, 
                           roots_bound_lower = NA, roots_bound_upper = NA){

  x_cheb <- cheb_scale(x_axis)
  cheb_t <- cheb_1st(x_axis, n)

  ## model calulations
  # least squares formula // ###
  cheb_coeff <- solve(t(cheb_t) %*% cheb_t) %*% t(cheb_t) %*% d

  # newton raphson method // package rootSolve
  extr <- uniroot.all(cheb_deriv_1st, cheb_coeff = cheb_coeff, lower = -1, upper = 1)
  extr_x <- if (length(extr) != 0) cheb_rescale(extr, x_axis = x_axis)

  if (!is.na(roots_bound_lower) | !is.na(roots_bound_upper)) {
    ind <- which(extr_x >= roots_bound_lower & extr_x <= roots_bound_upper)
    extr <- extr[ind]
    extr_x <- extr_x[ind]
  }
  extr.y <- if (length(extr_x) != 0) cheb_model_filter(x_axis = extr, cheb_coeff = cheb_coeff)
  extr_deriv_2nd <- if (length(extr) != 0) cheb_deriv_2nd(x_axis = extr, cheb_coeff = cheb_coeff)

  cheb_list <- list(extr_x = extr_x,
                    extr.y = extr.y,
                    extr_deriv_2nd = extr_deriv_2nd)
  return(cheb_list)
}


#' @title Routine to find maximum values of a curve
#' @param d Data series to be fitted
#' @param x_axis Arbitrary x-axis
#' @param n Order of the polynomial
#' @param maxima_bound_lower Lower boundary of the maximum value search
#' @param maxima_bound_upper Upper boundary of the maximum value search
#' @return list of calculated parameters like Chebyshev coefficients,
#'   values of the filtered model, values of the first and seconda derivation,
#'   postitions and values of maxima
#' @description
#' \code{cheb_find_max} fits a Chebyshev polynomial of order n to a data series 
#'   using the method of least squares, locates the extrema of the fit using 
#'   the Newton-Raphson-algorithm of the first derivation and checks for maxima 
#'   using the second derivative.
#' @export
#' @importFrom rootSolve uniroot.all
cheb_find_max <- function(d, x_axis, n, maxima_bound_lower = NA, maxima_bound_upper = NA){
  x_cheb <- cheb_scale(x_axis)
  cheb_t <- cheb_1st(x_axis, n)

  ## model calculations
  # method least squares // formula: #### 
  cheb_coeff <- solve(t(cheb_t) %*% cheb_t) %*% t(cheb_t) %*% d

  # newton raphson method // package rootSolve
  extr <- uniroot.all(cheb_deriv_1st, cheb_coeff = cheb_coeff, lower = -1, upper = 1)
  extr_x <- cheb_rescale(extr, x_axis = x_axis)

  if (!is.na(maxima_bound_lower) | !is.na(maxima_bound_upper)) {
    ind <- which(extr_x >= maxima_bound_lower & extr_x <= maxima_bound_upper)
    extr <- extr[ind]
    extr_x <- extr_x[ind]
  }

  # check for extreme values
  if (length(extr) != 0) {
    extr_deriv_2nd <- cheb_deriv_2nd(x_axis = extr, cheb_coeff = cheb_coeff)

    ind_max <- which(extr_deriv_2nd < 0)
    if (length(ind_max) != 0) {
      max_x <- extr_x[ind_max]
      max_y <- cheb_model_filter(x_axis = extr[ind_max], cheb_coeff = cheb_coeff)
    } else {
      max_x <- NA; max_y <- NA
    }
  } else {
    max_x <- NA; max_y <- NA;
  }
  max_list <- list(max_x = max_x,
                   max_y = max_y)
  return(max_list)
}
