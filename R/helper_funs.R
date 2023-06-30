#----------------------------------------------------------------------------
#' Simulate a transformed linear model
#'
#' Generate training data (X, y) and testing data (X_test, y_test)
#' for a transformed linear model. The covariates are correlated
#' Gaussian variables. Half of the true regression coefficients
#' are zero and the other half are one. There are multiple options
#' for the transformation, which define the support of the data (see below).
#'
#' @param n number of observations in the training data
#' @param p number of covariates
#' @param g_type type of transformation; must be one of
#' \code{beta}, \code{step}, or \code{box-cox}
#' @param n_test number of observations in the testing data
#' @param heterosked logical; if TRUE, simulate the latent data with heteroskedasticity
#' @param lambda Box-Cox parameter (only applies for \code{g_type = 'box-cox'})
#' @return a list with the following elements:
#' \itemize{
#' \item \code{y}: the response variable in the training data
#' \item \code{X}: the covariates in the training data
#' \item \code{y_test}: the response variable in the testing data
#' \item \code{X_test}: the covariates in the testing data
#' \item \code{beta_true}: the true regression coefficients
#' \item \code{g_true}: the true transformation, evaluated at y
#' }
#'
#' @details The transformations vary in complexity and support
#' for the observed data, and include the following options:
#' \code{beta} yields marginally Beta(0.1, 0.5) data
#' supported on [0,1]; \code{step} generates a locally-linear
#' inverse transformation and produces positive data; and \code{box-cox}
#' refers to the signed Box-Cox family indexed by \code{lambda},
#' which generates real-valued data with examples including identity,
#' square-root, and log transformations.
#'
#' @examples
#' # Simulate data:
#' dat = simulate_tlm(n = 100, p = 5, g_type = 'beta')
#' names(dat) # what is returned
#' hist(dat$y, breaks = 25) # marginal distribution
#'
#' @importFrom stats arima.sim approxfun pnorm qbeta rexp binomial coef dnorm ecdf fitted optim poly qnorm quantile rbeta rgamma rnorm runif
#' @export
simulate_tlm = function(n, p,
                        g_type = 'beta',
                        n_test = 1000,
                        heterosked = FALSE,
                        lambda = 1){
  #----------------------------------------------------------------------------
  # Checks:
  if(!is.element(g_type, c("beta", "step", "box-cox")))
    stop("The transformation must be one of 'beta', 'step', or 'box-cox'")

  if(p <= 1)
    stop("Include at least p=2 covariates")
  #----------------------------------------------------------------------------
  # Simulate a design matrix with correlated predictors:
  ar1 = 0.75
  X = cbind(1,
            t(apply(matrix(0, nrow = n, ncol = p), 1, function(x)
              arima.sim(n = p, list(ar = ar1), sd = sqrt(1-ar1^2)))))

  # Shuffle the columns:
  ind_shuff = sample(1:p)
  X[,-1] = X[,-1][,ind_shuff]

  # Covariates on the testing data:
  X_test = cbind(1,
                 t(apply(matrix(0, nrow = n_test, ncol = p), 1, function(x)
                   arima.sim(n = p, list(ar = ar1), sd = sqrt(1-ar1^2)))))
  # Match the column shuffling:
  X_test[,-1] = X_test[,-1][,ind_shuff]
  #----------------------------------------------------------------------------
  # True transformation:
  if(g_type =='beta'){
    g_inv_true = function(z)
      qbeta(pnorm(scale(z)),
            #(0.1, .5) and (0.5, 0.5) work well
            shape1 = 0.1, shape2 = 0.5)
  }
  if(g_type == 'step'){
    g_steps = rexp(n = 10)
    g_inv_true = function(z)
      approxfun(seq(-3, 3, length.out=10),
                cumsum(g_steps), rule = 2)(scale(z))
  }
  if(g_type == 'box-cox'){
    g_inv_true = function(z) g_inv_bc(scale(z),
                                      lambda = lambda)
  }

  # True coefficients:
  beta_true = c(0, # intercept is not needed
                rep(1, ceiling(p/2)),
                rep(0, floor(p/2)))
  #----------------------------------------------------------------------------
  # Simulate the training and testing observations:

  # Latent data (transformed scale):
  if(heterosked){
    z = X %*% beta_true*(1 + rnorm(n = n))
    z_test = X_test %*% beta_true*(1 + rnorm(n = n))
  } else {
    z = X %*% beta_true + rnorm(n = n)
    z_test = X_test %*% beta_true + rnorm(n = n_test)
  }

  # Observed data:
  y = g_inv_true(z)
  y_test = g_inv_true(z_test)

  return(list(
    y = y, X = X,
    y_test = y_test, X_test = X_test,
    beta_true = beta_true,
    g_true = z
  ))
}
#----------------------------------------------------------------------------
#' Plot point and interval predictions on testing data
#'
#' Given posterior predictive samples at \code{X_test},
#' plot the point and interval estimates and compare
#' to the actual testing data \code{y_test}.
#'
#' @param post_ypred \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' @param y_test \code{n_test} testing points
#' @param alpha_level alpha-level for prediction intervals
#' @return plot of the testing data, point and interval predictions,
#' and a summary of the empirical coverage
#'
#' @examples
#' \donttest{
#' # Simulate some data:
#' dat = simulate_tlm(n = 100, p = 5, g_type = 'step')
#'
#' # Fit a semiparametric Bayesian linear model:
#' fit = sblm(y = dat$y, X = dat$X, X_test = dat$X_test)
#'
#' # Evaluate posterior predictive means and intervals on the testing data:
#' plot_pptest(fit$post_ypred, dat$y_test,
#'             alpha_level = 0.10) # coverage should be about 90%
#' }
#' @importFrom graphics abline arrows lines
#' @export
plot_pptest = function(post_ypred,
                       y_test,
                       alpha_level = 0.10

    ){

  # Fitted values:
  y_hat = colMeans(post_ypred)

  # 100*(1 - alpha_level)% PI:
  pi_y = t(apply(post_ypred, 2, quantile,
                 c(alpha_level/2, 1 - alpha_level/2)))

  # Plot the testing values, fitted, values, and intervals:
  plot(y_test, y_test, type='n',
       ylim = range(pi_y, y_test),
       xlab = 'y_test', ylab = 'y_hat',
       main = paste('Prediction intervals: testing data'))
  abline(0,1) # reference line

  # PIs:
  suppressWarnings(
    arrows(y_test, pi_y[,1], y_test, pi_y[,2],
           length=0.15, angle=90, code=3, col='gray', lwd=2)
  ) # warnings occur when the intervals are very small

  # Fitted:
  lines(y_test, y_hat, type='p', pch=2) # plot the means

  # Coverage of PIs on testing data:
  mean((pi_y[,1] <= y_test)*(pi_y[,2] >= y_test))
}
#----------------------------------------------------------------------------
#' Box-Cox transformation
#'
#' Evaluate the Box-Cox transformation, which is a scaled power transformation
#' to preserve continuity in the index \code{lambda} at zero. Negative values are
#' permitted.
#'
#' @param t argument(s) at which to evaluate the function
#' @param lambda Box-Cox parameter
#' @return The evaluation(s) of the Box-Cox function at the given input(s) \code{t}.
#'
#' @note Special cases include
#' the identity transformation (\code{lambda = 1}),
#' the square-root transformation (\code{lambda = 1/2}),
#' and the log transformation (\code{lambda = 0}).
#'
#' @examples
#' # Log-transformation:
#' g_bc(1:5, lambda = 0); log(1:5)
#'
#' # Square-root transformation: note the shift and scaling
#' g_bc(1:5, lambda = 1/2); sqrt(1:5)
#'
#' @export
g_bc = function(t, lambda) {
  if(lambda == 0) {
    # (Signed) log-transformation
    sign(t)*log(abs(t))
  } else {
    # (Signed) Box-Cox-transformation
    (sign(t)*abs(t)^lambda - 1)/lambda
  }
}
#----------------------------------------------------------------------------
#' Inverse Box-Cox transformation
#'
#' Evaluate the inverse Box-Cox transformation. Negative values are permitted.
#'
#' @param s argument(s) at which to evaluate the function
#' @param lambda Box-Cox parameter
#' @return The evaluation(s) of the inverse Box-Cox function at the given input(s) \code{s}.
#'
#' @note Special cases include
#' the identity transformation (\code{lambda = 1}),
#' the square-root transformation (\code{lambda = 1/2}),
#' and the log transformation (\code{lambda = 0}).
#'
#'#' @examples
#' # (Inverse) log-transformation:
#' g_inv_bc(1:5, lambda = 0); exp(1:5)
#'
#' # (Inverse) square-root transformation: note the shift and scaling
#' g_inv_bc(1:5, lambda = 1/2); (1:5)^2
#'
#' @export
g_inv_bc = function(s, lambda) {
  if(lambda == 0) {
    # Inverse log-transformation
    exp(s)
  } else {
    # Inverse (signed) Box-Cox-transformation
    sign(lambda*s + 1)*abs(lambda*s+1)^(1/lambda)
  }
}
#' Rank-based estimation of the linear regression coefficients
#'
#' For a transformed Gaussian linear model, compute point estimates
#' of the regression coefficients. This  approach uses the ranks of the
#' data and does not require the transformation, but must expand the
#' sample size to \code{n^2} and thus can be slow.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors (should not include an intercept!)
#' @return the estimated linear coefficients
#'
#' @examples
#' # Simulate some data:
#' dat = simulate_tlm(n = 200, p = 10, g_type = 'step')
#'
#' # Point estimates for the linear coefficients:
#' theta_hat = suppressWarnings(
#'   rank_approx(y = dat$y,
#'               X = dat$X[,-1]) # remove intercept
#' ) # warnings occur from glm.fit (fitted probabilities 0 or 1)
#'
#' # Check: correlation with true coefficients
#' cor(dat$beta_true[-1], # excluding the intercept
#'     theta_hat)
#'
#' @importFrom stats glm
#' @export
rank_approx = function(y, X){

  # Data dimensions:
  n = length(y); p = ncol(X)

  # Create the binary response and X-matrix:
  r_stack = numeric(n^2)
  X_stack  = matrix(nrow = n^2, ncol = p)
  ind = 1
  for(i in 1:n){
    for(j in 1:n){
      r_stack[ind] = I(y[i] < y[j])  # flip from chenget al
      X_stack[ind,] = -(X[i,] - X[j,])/sqrt(2)
      ind = ind + 1
    }
  }

  # Fit the probit regression:
  fit = glm(r_stack ~ X_stack - 1,
            binomial(link = "probit"))

  # Return the estimates:
  coef(fit)

  #Sigma_hat = vcov(fit) # is this meaningful?
}
#----------------------------------------------------------------------------
#' Estimate the remaining time in the MCMC based on previous samples
#' @param nsi Current iteration
#' @param timer0 Initial timer value, returned from \code{proc.time()[3]}
#' @param nsims Total number of simulations
#' @param nrep Print the estimated time remaining every \code{nrep} iterations
#' @return Table of summary statistics using the function \code{summary}
computeTimeRemaining = function(nsi, timer0, nsims, nrep=1000){

  # Only print occasionally:
  if(nsi%%nrep == 0 || nsi==100) {
    # Current time:
    timer = proc.time()[3]

    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)

    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec

    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(round(secRemaining/60, 2), "minutes remaining"))
      } else print(paste(round(secRemaining), "seconds remaining"))
    }
  }
}
