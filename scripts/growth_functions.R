# Thes functions are from the growthcurver package by Kathleen Sprouffske which is no longer available on CRAN and thus the function is included here directly.

#' Number of Cells at Time t
#'
#' This function gives the number of cells or absorbance (N) at time t when
#' the parameters to the logistic equation are K, N0, and r.
#' @param k       The carrying capacity
#' @param n0      The initial population size (absorbance or individuals)
#' @param r       The exponential "growth rate"
#' @param t       The time at which you want to know N
#' @return        The number of cells, or N, at time t
#' @export
NAtT <- function(k, n0, r, t) {
  return( k / (1 + ((k - n0) / n0) * exp(-r * t)))
}


# Slope At Time t
#
# This function gives the growth rate (or slope of the curve) when
# the parameters to the logistic equation are K, N0, and r.
# @param k       The carrying capacity
# @param n0      The initial population size (absorbance or individuals)
# @param r       The exponential "growth rate"
# @param t       The time at which you want to know the growth rate
# @return        The number of cells, or N, at time t
SlopeAtT <- function(k, n0, r, t) {
  n <- NAtT(k, n0, r, t)
  return(r * n * (k - n) / k)
}

# Fastest Doubling Time
#
# This function gives you the maximum doubling time (DT) assuming exponential
# growth.
# @param r       The exponential "growth rate"
# @return        The maximum doubling time
MaxDt <- function(r) {
  return(log(2) / r)
}


# Doubling (Generation) Time at Time t
#
# This function gives you the doubling time (DT) at time t when the parameters
# of the logistic equation are K, N0, and r.
# @param k       A single integer specifying the carrying capacity
# @param n0      The initial population size
#                (in either absorbance or individuals)
# @param r       The exponential "growth rate"
# @param t       The time at which you want to know the doubling time
# @return        The doubling time at time t
DtAtT <- function(k, n0, r, t) {
  n_t <- NAtT(k, n0, r, t)
  n_halft <- 0.5 * n_t
  return(( 1 / r) * log((n_t * (k - n_halft)) / ((k - n_t) * n_halft)))
}


# Time at Inflection Point
#
# This function returns the time of the inflection point
# of the logistic equation with parameters K, N0, and r.
# @param k       A single integer specifying the carrying capacity
# @param n0      The initial population size
#                (in either absorbance or individuals)
# @param r       The exponential "growth rate"
# @return        The time of the inflection point, which occurrs when the
#                the population size N reaches half its maximum value, K
TAtInflection <- function(k, n0, r) {
  if (n0 == 0) {
    warning("Initial population size (n0) cannot be 0.")
    return(0)
  }
  t_inflection <- log(abs(k - n0) / n0) / r
  return(t_inflection)
}


# Area Under the Logistic Curve
#
# This function gives you the area under the curve from time t_min to t_max,
# when the parameters of the logistic equation are K, N0, and r. This value
# essentially combines the lag phase, growth rate, and carrying capacity
# into a single value.
# @param k       The carrying capacity
# @param n0      The initial population size
#                (in either absorbance or number of individuals)
# @param r       The exponential "growth rate"
# @param t_min   The time from which you want to know the area under the curve
#                (e.g., from time = t_min to t_max)
# @param t_max   The time to which you want to know the area under the curve
#                (e.g., from time = t_min to t_max)
# @return        The area under the curve for logistic equation with the
#                given parameters, for the specificed time range
AreaUnderCurve <- function(k, n0, r, t_min = 0, t_max) {
  auc_l <- stats::integrate(function(x) NAtT(k, n0, r, x), t_min, t_max)
  return(auc_l)
}



# Area Under the Empirical Curve
#
# This function returns the empirical "area under the curve". It uses the input
# data to do so (rather than using the logistic fit).
# @param data_t    A vector of timepoints (data_n must also
#                  be provided and be the same length).
# @param data_n    A vector of cell counts or absorbance readings.
# @param t_trim    Add up the area under the curve from the beginning to
#                  t_trim. Defaults to 0, which means don't trim.
# @return          The area under the curve
EmpiricalAreaUnderCurve <- function(data_t, data_n, t_trim = 0) {
  # make sure that both inputs are vectors
  if (!is.vector(data_t) | !is.vector(data_n)) {
    stop("Error: The input data (data_t and data_n) must be vectors.")
  }
  if (!is.numeric(data_t) | !is.numeric(data_n)) {
    stop("Error: The input data (data_t and data_n) must be numeric.")
  }
  if (t_trim > 0) {
    idx_to_keep <- data_t <= t_trim                # keep the early measurements
  }
  else {
    idx_to_keep < rep(TRUE, length(data_t))       # keep all measurements
  }
  
  x <- data_t[idx_to_keep]
  y <- data_n[idx_to_keep]
  n <- length(x)
  
  auc_e <- sum((x[2:n] - x[1:n-1]) * (y[2:n] + y[1:n-1]) /  2)
  return(auc_e)
}

#' Fits a logistic curve to data.
#'
#' This function fits a logistic curve to the supplied data, namely
#' n(t) = K / (1 + ( (K - N0) / N0) * exp(-r * t), where
#' N(t) is the number of cells (or density) at time t,
#' K is the carrying capacity,
#' N0 is the initial cell count or density, and
#' r is the "growth rate".
#' @param data_t    A vector of timepoints (data_n must also
#'                  be provided and be the same length).
#' @param data_n    A vector of cell counts or absorbance readings.
#' @return          An object of class nls.
#' @keywords        growth curves
FitLogistic <- function(data_t, data_n) {
  
  # make sure that the inputs are valid
  if (!is.vector(data_t) | !is.vector(data_n)) {
    stop("Error: The input data (data_t and data_n) must be vectors.")
  }
  if (!is.numeric(data_t) |!is.numeric(data_n)) {
    stop("Error: The input data (data_t and data_n) must be numeric.")
  }
  if (length(data_t) != length(data_n)) {
    stop("Error: The input data (data_t and data_n) must have the same length.")
  }
  
  # put together data
  d <- data.frame(cbind(data_t, data_n))
  names(d) <- c("t", "n")
  
  # make some guesses for the initial parameter values
  k_init <- max(data_n)   # carrying capacity is near the max
  n0_init <- min(data_n[data_n > 0])  # init population size is near the min
  
  # make an initial estimate for r
  glm_mod <- stats::glm(n / k_init ~ t,
                        family = stats::quasibinomial("logit"),
                        data = d)
  
  r_init <- stats::coef(glm_mod)[[2]]   # slope
  if (r_init <= 0) {
    # the slope should only be positive for a growing culture, so default
    # to something small
    r_init <- 0.001
  }
  
  suppressWarnings(
    nls_mod <- tryCatch(
      minpack.lm::nlsLM(n ~ k / (1 + ( (k - n0) / n0) * exp(-r * t)),
                        start = list(k = k_init,
                                     n0 = n0_init,
                                     r = r_init),
                        control = list(maxiter = 500),
                        lower = c(stats::median(data_n), 0, 0),
                        upper = c(Inf, max(data_n), Inf),
                        data = d),
      error = function(e) {
        stop("Error: Growthcurver FitLogistic cannot fit data.")
      }
    )
  )
  return(nls_mod)
}

#' Summarize Growth Curves
#'
#' This function finds the parameters that describe the input data's growth.
#' It does so by fitting the logistic curve to your growth curve measurements.
#'
#' The logistic curve equation is
#' \deqn{N_t = \frac{N_0 K} {N_0 + (K-N_0)e^{-rt}}}{N_t = N_0 K / (N_0 + (K - N_0) exp(-rt))}
#' where \eqn{N_t} is the number
#' of cells (or the absorbance reading) at time t, \eqn{N_0} is the initial
#' cell count (or absorbance reading), K is the carrying capacity, and r is the
#' growth rate.
#'
#' The fitness proxies returned are the parameters of the logistic equation
#' and the area under the curve (a measure that integrates the effects
#' of \eqn{N_0}, K, and r). See \code{\link{gcfit}} for more documentation on these.
#' @param data_t     A vector of timepoints (data_n must also
#'                   be provided and be the same length).
#' @param data_n     A vector of cell counts or absorbance readings.
#' @param t_trim     Measurements taken after this time should not be included
#'                   in fitting the curve. If stationary phase is variable,
#'                   this may give you a better fit. A value of 0 means no
#'                   trimming. Defaults to no trimming (0).
#' @param bg_correct The background correction method to use. No background
#'                   correction is performed for the default "none". Specifying
#'                   "min" subtracts the smallest value in a column from all the
#'                   rows in that column, and specifying "blank" subtracts
#'                   the values from the blank vector from the data_n vector.
#' @param blank      A vector of absorbance readings from a blank well
#'                   (typically contains only media) used for background
#'                   correction. The corresponding blank value is subtracted
#'                   from the data_n vector for each timepoint. Defaults to NA.
#' @return           An object of type gcfit containing the "fitness" proxies,
#'                   as well as the input data and the fitted model.
#' @seealso
#' See the accompanying Vignette for an example of how to use and interpret
#' SummarizeGrowth. \url{bit.ly/1p7w6dJ}.
#'
#' See also \code{\link{gcfit}}.
#' @examples
#' # We can check that the parameters that are found are the same
#' # as we use to generate fake experimental data. To do so, let's first
#' # generate the "experimental" data using the logistic equation,
#' # e.g., absorbance readings from a single well in a plate reader over time.
#'
#' k_in <- 0.5   # the initial carrying capacity
#' n0_in <- 1e-5 # the initial absorbance reading
#' r_in <- 1.2   # the initial growth rate
#' N <- 50       # the number of "measurements" collected during the growth
#'               # curve experiment
#'
#' data_t <- 0:N * 24 / N   # the times the measurements were made (in hours)
#' data_n <- NAtT(k = k_in, n0 = n0_in, r = r_in, t = data_t) # the measurements
#'
#' # Now summarize the "experimental" growth data that we just generated
#' gc <- SummarizeGrowth(data_t, data_n)
#'
#' # Get the possible metrics for fitness proxies
#' gc$vals$r           # growth rate is a common choice for fitness
#' gc$vals$t_gen       # doubling time, or generation time, is also common
#' gc$vals$k
#' gc$vals$n0
#' gc$vals$auc_l
#' gc$vals$auc_e
#' gc$vals$t_mid
#'
#' # Compare the data with the fit visually by plotting it
#' plot(gc)
#'
#' @export
SummarizeGrowth <- function(data_t, data_n, t_trim = 0,
                            bg_correct = "min", blank = NA) {

  # make sure that both inputs are vectors
  if (is.list(data_t) == TRUE) {
    tryCatch( {data_t <- unlist(data_t)},
              error = function(e) {stop("data_t is not a vector")}
            )
  }
  if (is.list(data_n) == TRUE) {
    tryCatch( {data_n <- unlist(data_n)},
              error = function(e) {stop("data_n is not a vector")}
            )
  }
  stopifnot(is.vector(data_t))
  stopifnot(is.vector(data_n))

  # make sure that the inputs are valid
  if (!is.numeric(data_t) |!is.numeric(data_n)) {
    stop("Error: The input data (data_t and data_n) must be numeric.")
  }
  if (length(data_t) != length(data_n)) {
    stop("Error: The input data (data_t and data_n) must have the same number of rows")
  }

  # make sure that the background correction method is valid
  if (!bg_correct %in% c("none", "min", "blank")) {
    stop(paste0(bg_correct, "is not a valid option for bg_correct"))
  }

  # check for correctness of the blank (background correction) vector
  if (bg_correct == "blank") {
    if (is.list(blank) == TRUE) {
      tryCatch( {blank <- unlist(blank)},
                error = function(e) {stop("blank is not a vector")}
      )
      stopifnot(is.vector(blank))

      if (!is.numeric(blank)) {
        stop("Error: The blank data must be numeric.")
      }
      if (length(blank) != length(data_n)) {
        stop("Error: The input data (data_n) and the background correction data (blank) must have the same number of rows.")
      }
    }
  }

  # check t_trim parameter and set dependent variables
  if (t_trim > 0) {
    t_max <- t_trim
    data_n <- data_n[data_t < t_trim]
    data_t <- data_t[data_t < t_trim]
    if (bg_correct == "blank") {
      blank <- blank[data_t < t_trim]
    }
  } else {
    t_max <- max(data_t, na.rm=TRUE)
  }

  #do the background correction
  if (bg_correct == "blank") {
    data_n <- data_n - blank
    data_n[data_n < 0] <- 0    # ensure readings are at least 0
  }
  else if (bg_correct == "min") {
    data_n <- data_n - min(data_n)
  }

  tryCatch(
    # code block
    {log_mod = FitLogistic(data_t, data_n)},
    # error handling block
    error = function(e) {}
    )

  # the data did not fit a logistic model
  if (exists("log_mod") == FALSE) {
    vals <- c(0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0,
              0, 0, 0, 0)
  }
  else {
    p <- summary(log_mod)$coefficients
    k <- p[1]
    k_se <- p[4]
    k_p <- p[10]
    n0 <- p[2]
    n0_se <- p[5]
    n0_p <- p[11]
    r <- p[3]
    r_se <- p[6]
    r_p <- p[12]

    # get the inflection point, DT, auc, sigma, df
    t_inflection <- TAtInflection(k, n0, r)

    DT <- MaxDt(r)
    sigma <- summary(log_mod)$sigma
    df <- summary(log_mod)$df[2]

    auc_l <- AreaUnderCurve(k, n0, r, 0, t_max)$value
    auc_e <- EmpiricalAreaUnderCurve(data_t, data_n, t_max)
    vals <- c(k, k_se, k_p, n0, n0_se, n0_p,
              r, r_se, r_p, sigma, df,
              t_inflection, DT, auc_l, auc_e)
  }

  val_names <- c("k", "k_se", "k_p", "n0", "n0_se", "n0_p",
                 "r", "r_se", "r_p", "sigma", "df",
                 "t_mid", "t_gen", "auc_l", "auc_e")
  vals <- stats::setNames(as.list(vals), val_names)
  class(vals) <- "gcvals"

  if (exists("log_mod") == FALSE) {
    vals$note <- "cannot fit data"
    log_mod <- ""
  }
  else if (k < n0) {
    vals$note <- "questionable fit (k < n0)"
  }
  else if (t_inflection < 0) {
    vals$note <- "questionable fit"
  }
  else {
    vals$note <- ""
  }

  ret <- list("vals" = vals, "model" = log_mod,
              "data"=list("t" = data_t, "N" = data_n))
  class(ret) <- "gcfit"
  return(ret)
}


