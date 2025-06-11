#---------------------------------------------------------------
# Source functions for semiparametric Bayesian regression
#---------------------------------------------------------------
#' Semiparametric Bayesian linear model
#'
#' Monte Carlo sampling for Bayesian linear regression with an
#' unknown (nonparametric) transformation. A g-prior is assumed
#' for the regression coefficients.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors (no intercept)
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param psi prior variance (g-prior)
#' @param laplace_approx logical; if TRUE, use a normal approximation
#' to the posterior in the definition of the transformation;
#' otherwise the prior is used
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param approx_g logical; if TRUE, apply large-sample
#' approximation for the transformation
#' @param nsave number of Monte Carlo simulations
#' @param ngrid number of grid points for inverse approximations
#' @param verbose logical; if TRUE, print time remaining
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{X_test}
#' \item \code{post_theta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_ypred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{model}: the model fit (here, \code{sblm})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed linear model using Monte Carlo (not MCMC) sampling.
#' The transformation is modeled as unknown and learned jointly
#' with the regression coefficients (unless \code{approx_g = TRUE}, which then uses
#' a point approximation). This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' The results are typically unchanged whether \code{laplace_approx} is TRUE/FALSE;
#' setting it to TRUE may reduce sensitivity to the prior, while setting it to FALSE
#' may speed up computations for very large datasets. By default, \code{fixedX} is
#' set to FALSE for smaller datasets (\code{n < 500}) and TRUE for larger datasets (\code{n >= 500}).
#'
#' @note The location (intercept) and scale (\code{sigma_epsilon}) are
#' not identified, so any intercepts in \code{X} and \code{X_test} will
#' be removed. The model-fitting *does* include an internal location-scale
#' adjustment, but the function only outputs inferential summaries for the
#' identifiable parameters.
#'
#' @examples
#' \donttest{
#' # Simulate some data:
#' dat = simulate_tlm(n = 100, p = 5, g_type = 'step')
#' y = dat$y; X = dat$X # training data
#' y_test = dat$y_test; X_test = dat$X_test # testing data
#'
#' hist(y, breaks = 25) # marginal distribution
#'
#' # Fit the semiparametric Bayesian linear model:
#' fit = sblm(y = y, X = X, X_test = X_test)
#' names(fit) # what is returned
#'
#' # Note: this is Monte Carlo sampling...no need for MCMC diagnostics!
#'
#' # Evaluate posterior predictive means and intervals on the testing data:
#' plot_pptest(fit$post_ypred, y_test,
#'             alpha_level = 0.10) # coverage should be about 90%
#'
#' # Check: correlation with true coefficients
#' cor(dat$beta_true, coef(fit))
#'
#' # Summarize the transformation:
#' y0 = sort(unique(y)) # posterior draws of g are evaluated at the unique y observations
#' plot(y0, fit$post_g[1,], type='n', ylim = range(fit$post_g),
#'      xlab = 'y', ylab = 'g(y)', main = "Posterior draws of the transformation")
#' temp = sapply(1:nrow(fit$post_g), function(s)
#'   lines(y0, fit$post_g[s,], col='gray')) # posterior draws
#' lines(y0, colMeans(fit$post_g), lwd = 3) # posterior mean
#' lines(y, dat$g_true, type='p', pch=2) # true transformation
#' legend('bottomright', c('Truth'), pch = 2) # annotate the true transformation
#'
#' # Posterior predictive checks on testing data: empirical CDF
#' y0 = sort(unique(y_test))
#' plot(y0, y0, type='n', ylim = c(0,1),
#'      xlab='y', ylab='F_y', main = 'Posterior predictive ECDF')
#' temp = sapply(1:nrow(fit$post_ypred), function(s)
#'   lines(y0, ecdf(fit$post_ypred[s,])(y0), # ECDF of posterior predictive draws
#'         col='gray', type ='s'))
#' lines(y0, ecdf(y_test)(y0),  # ECDF of testing data
#'      col='black', type = 's', lwd = 3)
#' }
#' @export
sblm = function(y, X, X_test = X,
                psi = length(y),
                laplace_approx = TRUE,
                fixedX = (length(y) >= 500),
                approx_g = FALSE,
                nsave = 1000,
                ngrid = 100,
                verbose = TRUE){

  # For testing:
  # X_test = X; psi = length(y); laplace_approx = TRUE; fixedX = FALSE; approx_g = FALSE; nsave = 1000; verbose = TRUE; ngrid = 100

  # Initial checks:
  if(!is.matrix(X)) stop("X must be a matrix (rows = observations, columns = variables)")
  if(!is.matrix(X_test)) stop("X_test must be a matrix (rows = observations, columns = variables)")
  if(ncol(X) != ncol(X_test)) stop('X_test and X must have the same number of columns (variables)')
  if(nrow(X) != length(y)) stop('the length of y must equal nrow(X)')

  # Data dimensions:
  n = length(y) # number of observations
  p = ncol(X) # number of variables
  n_test = nrow(X_test) # number of testing data points

  # Make sure that X does not include an intercept
  #   Note: one will be added later for location (+ scale) adjustment,
  #   but this is not part of the target parameter 'theta'
  is_intercept = apply(X, 2, function(x) all(x == 1))
  if(any(is_intercept)){

    # Quick check: need more than just an intercept here!
    if(p==1) stop('X must include at least one (non-intercept) covariate')

    # Remove the intercept:
    X = X[, !is_intercept, drop=FALSE]

    # Update dimensions:
    p = ncol(X)
  }

  # For X_test: remove intercepts anywhere, then make X_test[,1] an intercept
  #   (This will match X later)
  is_intercept_test = apply(X_test, 2, function(x) all(x == 1))
  X_test = cbind(1, X_test[,!is_intercept_test])

  # Requirement for the g-prior:
  if(p >= n) stop('The g-prior requires p < n')

  # Hyperparameters for Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)

  # Evaluate at the unique y-values:
  y0 = sort(unique(y))
  Fy_eval = Fy(y0)

  # Matrix quantities for the initialization:
  XtXinv = chol2inv(chol(crossprod(X))) # inverse of (X'X)
  xt_Sigma_x = sapply(1:n, function(i)
    crossprod(X[i,], XtXinv)%*%X[i,])

  # Grid of values for the CDF of z (based on the prior)
  z_grid = sort(unique(
    sapply(range(psi*xt_Sigma_x), function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = ngrid),
            mean = 0, # assuming prior mean zero
            sd = sqrt(1 + xtemp))
    })
  ))

  # Define the moments of the CDF of z:
  if(laplace_approx){
    # Use a normal approximation for the posterior of theta

    # Recurring terms:
    Sigma_hat_unscaled = psi/(1+psi)*XtXinv # unscaled covariance (w/o sigma)
    xt_Sigma_hat_unscaled_x = sapply(1:n, function(i)
      crossprod(X[i,], Sigma_hat_unscaled)%*%X[i,])

    # First pass: fix Fz() = qnorm(), initialize coefficients
    z = qnorm(Fy(y))
    theta_hat = Sigma_hat_unscaled%*%crossprod(X, z) # point estimate

    # Alternative rank-based approaches (SLOW! but similar...)
    # theta_hat = rank_approx(y, X) # alternative rank-based approach
    # theta_hat = coef(Rfit::rfit(y  ~ X)) # faster, but less accurate
    # that function also comes with a sd estimate

    # Second pass: update g(), then update coefficients
    # Moments of Z|X:
    mu_z = X%*%theta_hat
    sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

    # CDF of z:
    Fz_eval = Fz_fun(z = z_grid,
                     weights = rep(1/n, n),
                     mean_vec = mu_z,
                     sd_vec = sigma_z)

    # Check: update the grid if needed
    zcon = contract_grid(z = z_grid,
                         Fz = Fz_eval,
                         lower = 0.001, upper =  0.999)
    z_grid = zcon$z; Fz_eval = zcon$Fz

    # Transformation:
    g = g_fun(y = y0,
              Fy_eval = Fy_eval,
              z = z_grid,
              Fz_eval = Fz_eval)

    # Updated coefficients:
    z = g(y) # update latent data
    theta_hat = Sigma_hat_unscaled%*%crossprod(X, z) # updated coefficients

    # Moments of Z|X:
    mu_z = X%*%theta_hat
    #sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x) # no need to update

  } else {

    # Prior mean is zero:
    mu_z = rep(0, n)

    # Marginal SD based on prior:
    sigma_z = sqrt(1 + psi*xt_Sigma_x)
  }

  # Define the CDF of z:
  Fz_eval = Fz_fun(z = z_grid,
                   weights = rep(1/n, n),
                   mean_vec = mu_z,
                   sd_vec = sigma_z)

  # Check: update the grid if needed
  zcon = contract_grid(z = z_grid,
                       Fz = Fz_eval,
                       lower = 0.001, upper =  0.999)
  z_grid = zcon$z; Fz_eval = zcon$Fz

  # Compute the transformation:
  g = g_fun(y = y0, Fy_eval = Fy_eval,
            z = z_grid, Fz_eval = Fz_eval)

  # Latent data:
  z = g(y)

  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)
  #----------------------------------------------------------------------------
  # Add an intercept:
  X = cbind(1, X) # update to include intercept

  # One-time computing cost to sample (theta, sigma) efficiently:
  chXtX_psi = sqrt((1+psi)/(psi))*chol(crossprod(X))
  Xtz = crossprod(X, z) # this changes with z = g(y) (unless approx_g = TRUE)
  #----------------------------------------------------------------------------
  # Store MC output:
  post_theta = array(NA, c(nsave, p))
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))

  # Run the MC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nsave){

    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(!approx_g){

      # Bayesian bootstrap for the CDFs

      # BB CDF of y:
      Fy_eval = bb(y)(y0) # applied to y, evaluated at y0

      # BB CDF of z (only need if X is random)
      if(!fixedX){

        # Dirichlet(1) weights for x:
        weights_x = rgamma(n = n, shape = 1)
        weights_x  = weights_x/sum(weights_x)

        # BB CDF of z:
        Fz_eval = Fz_fun(z = z_grid,
                         weights = weights_x,
                         mean_vec = mu_z,
                         sd_vec = sigma_z)

      } # otherwise, the weights are 1/n (how we initialized)

      # Compute the transformation:
      g = g_fun(y = y0, Fy_eval = Fy_eval,
                z = z_grid, Fz_eval = Fz_eval)

      # Update z:
      z = g(y)
      Xtz = crossprod(X, z)

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = y_grid)
    }
    #----------------------------------------------------------------------------
    # Recurring term for Blocks 2-3:
    fsolve_theta = forwardsolve(t(chXtX_psi), Xtz)

    # Block 2: sample the scale adjustment (SD)
    SSR_psi = sum(z^2) - crossprod(Xtz, backsolve(chXtX_psi, fsolve_theta)) # SSR_psi = sum(z^2) - psi/(psi+1)*crossprod(Xtz, XtXinv%*%Xtz)
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = a_sigma + n/2,
                                  rate = b_sigma + SSR_psi/2))
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    theta = backsolve(chXtX_psi/sigma_epsilon,
                      fsolve_theta/sigma_epsilon + rnorm(p+1))
    # Previously:
    # ch_Q = chol(1/sigma_epsilon^2*(1+psi)/(psi)*XtX)
    # ell_theta = 1/sigma_epsilon^2*Xtz
    # theta = backsolve(ch_Q, forwardsolve(t(ch_Q), ell_theta) + rnorm(p + 1))
    #----------------------------------------------------------------------------
    # Store the MC:

    # Posterior samples of the model parameters:
    post_theta[nsi,] = theta[-1]/sigma_epsilon # omit the intercept & undo scaling (not identified)

    # Predictive samples of ytilde:
    ztilde = X_test%*%theta + sigma_epsilon*rnorm(n = n_test)
    post_ypred[nsi,] = g_inv(ztilde)

    # Posterior samples of the transformation:
    post_g[nsi,] = (g(y0) - theta[1])/sigma_epsilon # undo location/scale (not identified)
    #----------------------------------------------------------------------------
    if(verbose) computeTimeRemaining(nsi, timer0, nsave)
  }
  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer0
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  return(list(
    coefficients = colMeans(post_theta),
    fitted.values = colMeans(post_ypred),
    post_theta = post_theta,
    post_ypred = post_ypred,
    post_g = post_g,
    model = 'sblm', y = y, X = X[,-1], X_test = X_test[,-1], psi = psi, laplace_approx = laplace_approx, fixedX = fixedX, approx_g = approx_g))
}
#---------------------------------------------------------------
#' Semiparametric Bayesian spline model
#'
#' Monte Carlo sampling for Bayesian spline regression with an
#' unknown (nonparametric) transformation. Cubic B-splines are
#' used with a prior that penalizes roughness.
#'
#' @param y \code{n x 1} response vector
#' @param x \code{n x 1} vector of observation points; if NULL, assume equally-spaced on [0,1]
#' @param x_test \code{n_test x 1} vector of testing points; if NULL, assume equal to \code{x}
#' @param psi prior variance (inverse smoothing parameter); if NULL,
#' sample this parameter
#' @param laplace_approx logical; if TRUE, use a normal approximation
#' to the posterior in the definition of the transformation;
#' otherwise the prior is used
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param approx_g logical; if TRUE, apply large-sample
#' approximation for the transformation
#' @param nsave number of Monte Carlo simulations
#' @param ngrid number of grid points for inverse approximations
#' @param verbose logical; if TRUE, print time remaining
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{x_test}
#' \item \code{post_theta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_ypred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at \code{x_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{model}: the model fit (here, \code{sbsm})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed spline regression model using Monte Carlo (not MCMC) sampling.
#' The transformation is modeled as unknown and learned jointly
#' with the regression function (unless \code{approx_g = TRUE}, which then uses
#' a point approximation). This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' The results are typically unchanged whether \code{laplace_approx} is TRUE/FALSE;
#' setting it to TRUE may reduce sensitivity to the prior, while setting it to FALSE
#' may speed up computations for very large datasets. By default, \code{fixedX} is
#' set to FALSE for smaller datasets (\code{n < 500}) and TRUE for larger datasets (\code{n >= 500}).
#'
#' @examples
#' \donttest{
#' # Simulate some data:
#' n = 200 # sample size
#' x = sort(runif(n)) # observation points
#'
#' # Transform a noisy, periodic function:
#' y = g_inv_bc(
#'   sin(2*pi*x) + sin(4*pi*x) + rnorm(n),
#'              lambda = .5) # Signed square-root transformation
#'
#' # Fit the semiparametric Bayesian spline model:
#' fit = sbsm(y = y, x = x)
#' names(fit) # what is returned
#'
#' # Note: this is Monte Carlo sampling...no need for MCMC diagnostics!
#'
#' # Plot the model predictions (point and interval estimates):
#' pi_y = t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
#' plot(x, y, type='n', ylim = range(pi_y,y),
#'      xlab = 'x', ylab = 'y', main = paste('Fitted values and prediction intervals'))
#' polygon(c(x, rev(x)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
#' lines(x, y, type='p') # observed points
#' lines(x, fitted(fit), lwd = 3) # fitted curve
#' }
# #' @importFrom spikeSlabGAM sm
#' @export
sbsm = function(y, x = NULL,
                x_test = x,
                psi = NULL,
                laplace_approx = TRUE,
                fixedX = (length(y) >= 500),
                approx_g = FALSE,
                nsave = 1000,
                ngrid = 100,
                verbose = TRUE){

  # For testing:
  # psi = length(y); laplace_approx = TRUE; fixedX = FALSE; approx_g = FALSE; nsave = 1000; verbose = TRUE; ngrid = 100

  # Library required here:
  if (!requireNamespace("spikeSlabGAM", quietly = TRUE)) {
    stop(
      "Package \"spikeSlabGAM\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Data dimensions:
  n = length(y)

  # Observation points:
  if(is.null(x)) x = seq(0, 1, length=n)

  # Initial checks:
  if(length(x) != n) stop('x and y must have the same number of observations')

  # Recale to [0,1]:
  x = (x - min(x))/(max(x) - min(x))
  x_test = (x_test - min(x_test))/(max(x_test) - min(x_test))

  # Hyperparameters for Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001
  #----------------------------------------------------------------------------
  # Orthogonalized P-spline and related quantities:
  X = cbind(poly(x, 1), spikeSlabGAM::sm(x)) # omit intercept to initialize
  #X = X/sqrt(sum(diag(crossprod(X))))
  diagXtX = colSums(X^2)
  p = length(diagXtX)

  # Recurring term:
  xt_Sigma_x = rowSums(X^2) # sapply(1:n, function(i) sum(X[i,]^2/diagXtX))

  # Smoothing parameter:
  if(is.null(psi)){

    # Flag to sample:
    sample_psi = TRUE

    # Initialize:
    psi = n

  } else sample_psi = FALSE
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)

  # Evaluate at the unique y-values:
  y0 = sort(unique(y))
  Fy_eval = Fy(y0)

  # Grid of values for the CDF of z:
  z_grid = sort(unique(
    sapply(range(psi*xt_Sigma_x), function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = ngrid),
            mean = 0, # assuming prior mean zero
            sd = sqrt(1 + xtemp))
    })
  ))

  # Define the moments of the CDF of z:
  if(laplace_approx){
    # Use a normal approximation for the posterior of theta

    # Recurring terms:
    diag_Sigma_hat_unscaled = 1/(diagXtX + 1/psi)
    xt_Sigma_hat_unscaled_x = colSums(t(X^2)*diag_Sigma_hat_unscaled)
    #xt_Sigma_hat_unscaled_x = sapply(1:n, function(i)
    #  crossprod(X[i,], diag(diag_Sigma_hat_unscaled))%*%X[i,])

    # First pass: fix Fz() = qnorm(), initialize coefficients
    z = qnorm(Fy(y))
    theta_hat = diag_Sigma_hat_unscaled*crossprod(X, z) # point estimate

    # Second pass: update g(), then update coefficients
    # Moments of Z|X:
    mu_z = X%*%theta_hat
    sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

    # CDF of z:
    Fz_eval = Fz_fun(z = z_grid,
                     weights = rep(1/n, n),
                     mean_vec = mu_z,
                     sd_vec = sigma_z)

    # Check: update the grid if needed
    zcon = contract_grid(z = z_grid,
                         Fz = Fz_eval,
                         lower = 0.001, upper =  0.999)
    z_grid = zcon$z; Fz_eval = zcon$Fz

    # Transformation:
    g = g_fun(y = y0,
              Fy_eval = Fy_eval,
              z = z_grid,
              Fz_eval = Fz_eval)

    # Updated coefficients:
    z = g(y) # update latent data
    theta_hat = diag_Sigma_hat_unscaled*crossprod(X, z) # updated coefficients

    # Moments of Z|X:
    mu_z = X%*%theta_hat
    #sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

  } else {

    # Prior mean is zero:
    mu_z = rep(0, n)

    # Marginal SD based on prior:
    sigma_z = sqrt(1 + psi*xt_Sigma_x)
  }

  # Define the CDF of z:
  Fz_eval = Fz_fun(z = z_grid,
                   weights = rep(1/n, n),
                   mean_vec = mu_z,
                   sd_vec = sigma_z)

  # Check: update the grid if needed
  zcon = contract_grid(z = z_grid,
                       Fz = Fz_eval,
                       lower = 0.001, upper =  0.999)
  z_grid = zcon$z; Fz_eval = zcon$Fz

  # Compute the transformation:
  g = g_fun(y = y0, Fy_eval = Fy_eval,
            z = z_grid, Fz_eval = Fz_eval)

  # Latent data:
  z = g(y)

  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)
  #----------------------------------------------------------------------------
  # Add an intercept:
  X = cbind(1/sqrt(n), X)
  diagXtX = c(1, diagXtX)  # update
  Xtz = crossprod(X, z) # changes w/ z (unless approx_g = TRUE)
  #----------------------------------------------------------------------------
  # Store MC output:
  post_theta = array(NA, c(nsave, p))
  post_ypred = array(NA, c(nsave, length(x_test)))
  post_g = array(NA, c(nsave, length(y0)))
  post_psi = rep(NA, nsave)

  # Run the MC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nsave){
    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(!approx_g){

      # Bayesian bootstrap for the CDFs

      # BB CDF of y:
      Fy_eval = bb(y)(y0) # applied to y, evaluated at y0

      # BB CDF of z (only need if X is random)
      if(!fixedX){

        # Dirichlet(1) weights for x:
        weights_x = rgamma(n = n, shape = 1)
        weights_x  = weights_x/sum(weights_x)

        # CDF of Z evaluated on the grid:
        Fz_eval = Fz_fun(z = z_grid,
                         weights = weights_x,
                         mean_vec = mu_z,
                         sd_vec = sigma_z)
      } # otherwise, the weights are 1/n (how we initialized)

      # Compute the transformation:
      g = g_fun(y = y0, Fy_eval = Fy_eval,
                z = z_grid, Fz_eval = Fz_eval)

      # Update z:
      z = g(y)
      Xtz = crossprod(X, z)

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = y_grid)
    }
    #----------------------------------------------------------------------------
    # Block 2: sample the scale adjustment (SD)
    # SSR_psi = sum(z^2) - crossprod(z, X%*%solve(crossprod(X) + diag(1/psi, p+1))%*%crossprod(X,z))
    SSR_psi = sum(z^2) - crossprod(1/sqrt(diagXtX + 1/psi)*Xtz)
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = a_sigma + n/2,
                                  rate = b_sigma + SSR_psi/2))
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    Q_theta = 1/sigma_epsilon^2*(diagXtX + 1/psi)
    ell_theta = 1/sigma_epsilon^2*Xtz
    theta = rnorm(n = p+1,
                 mean = Q_theta^-1*ell_theta,
                 sd = sqrt(Q_theta^-1))
    #----------------------------------------------------------------------------
    # Block 4: sample the smoothing parameter
    if(sample_psi){
      psi = 1/rgamma(n = 1,
                     shape = 0.01 + (p+1)/2,
                     rate = 0.01 + sum(theta^2)/(2*sigma_epsilon^2))
    }
    #----------------------------------------------------------------------------
    # Store the MC:

    # Posterior samples of the model parameters:
    post_theta[nsi,] = theta[-1]/sigma_epsilon # omit the intercept & undo scaling (not identified)

    # Predictive samples of ytilde:
    #   Note: it's easier/faster to just smooth for the testing points
    #   (the orthogonalized basis is a pain to recompute)
    ztilde = stats::spline(x = x, y = X%*%theta, xout = x_test)$y +
      sigma_epsilon*rnorm(n = length(x_test))
    post_ypred[nsi,] = g_inv(ztilde)

    # Posterior samples of the transformation:
    post_g[nsi,] = (g(y0) - theta[1])/sigma_epsilon # undo location/scale (not identified)

    # Posterior samples of the prior variance (inverse smoothing parameter)
    post_psi[nsi] = psi
    #----------------------------------------------------------------------------
    if(verbose) computeTimeRemaining(nsi, timer0, nsave)
  }
  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer0
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  return(list(
    coefficients = colMeans(post_theta),
    fitted.values = colMeans(post_ypred),
    post_theta = post_theta,
    post_ypred = post_ypred,
    post_g = post_g,  post_psi = post_psi,
    model = 'sbsm', y = y, X = X[,-1], sample_psi = sample_psi, laplace_approx = laplace_approx, fixedX = fixedX, approx_g = approx_g))
}
#---------------------------------------------------------------
#' Semiparametric Bayesian Gaussian processes
#'
#' Monte Carlo sampling for Bayesian Gaussian process regression with an
#' unknown (nonparametric) transformation.
#'
#' @param y \code{n x 1} response vector
#' @param locs \code{n x d} matrix of locations
#' @param X \code{n x p} design matrix; if unspecified, use intercept only
#' @param covfun_name string name of a covariance function; see ?GpGp
#' @param locs_test \code{n_test x d} matrix of locations
#' at which predictions are needed; default is \code{locs}
#' @param X_test \code{n_test x p} design matrix for test data;
#' default is \code{X}
#' @param nn number of nearest neighbors to use; default is 30
#' (larger values improve the approximation but increase computing cost)
#' @param emp_bayes logical; if TRUE, use a (faster!) empirical Bayes
#' approach for estimating the mean function
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param approx_g logical; if TRUE, apply large-sample
#' approximation for the transformation
#' @param nsave number of Monte Carlo simulations
#' @param ngrid number of grid points for inverse approximations
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the estimated regression coefficients
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{locs_test}
#' \item \code{fit_gp} the fitted \code{GpGp_fit} object, which includes
#' covariance parameter estimates and other model information
#' \item \code{post_ypred}: \code{nsave x ntest} samples
#' from the posterior predictive distribution at \code{locs_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{model}: the model fit (here, \code{sbgp})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides Bayesian inference for a
#' transformed Gaussian process model using Monte Carlo (not MCMC) sampling.
#' The transformation is modeled as unknown and learned jointly
#' with the regression function (unless \code{approx_g = TRUE}, which then uses
#' a point approximation). This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' The results are typically unchanged whether \code{laplace_approx} is TRUE/FALSE;
#' setting it to TRUE may reduce sensitivity to the prior, while setting it to FALSE
#' may speed up computations for very large datasets. For computational efficiency,
#' the Gaussian process parameters are fixed at point estimates, and the latent Gaussian
#' process is only sampled when \code{emp_bayes = FALSE}. However, the uncertainty
#' from this term is often negligible compared to the observation errors, and the
#' transformation serves as an additional layer of robustness. By default, \code{fixedX} is
#' set to FALSE for smaller datasets (\code{n < 500}) and TRUE for larger datasets (\code{n >= 500}).
#'
#' @examples
#' \donttest{
#' # Simulate some data:
#' n = 200 # sample size
#' x = seq(0, 1, length = n) # observation points
#'
#' # Transform a noisy, periodic function:
#' y = g_inv_bc(
#'   sin(2*pi*x) + sin(4*pi*x) + rnorm(n),
#'              lambda = .5) # Signed square-root transformation
#'
#' # Package we use for fast computing w/ Gaussian processes:
#' library(GpGp)
#'
#' # Fit the semiparametric Bayesian Gaussian process:
#' fit = sbgp(y = y, locs = x)
#' names(fit) # what is returned
#' coef(fit) # estimated regression coefficients (here, just an intercept)
#' class(fit$fit_gp) # the GpGp object is also returned
#'
#' # Plot the model predictions (point and interval estimates):
#' pi_y = t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
#' plot(x, y, type='n', ylim = range(pi_y,y),
#'      xlab = 'x', ylab = 'y', main = paste('Fitted values and prediction intervals'))
#' polygon(c(x, rev(x)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
#' lines(x, y, type='p') # observed points
#' lines(x, fitted(fit), lwd = 3) # fitted curve
#' }
# #' @import GpGp fields
#' @export
sbgp = function(y, locs,
                X = NULL,
                covfun_name = "matern_isotropic",
                locs_test = locs,
                X_test = NULL,
                nn = 30,
                emp_bayes = TRUE,
                fixedX = (length(y) >= 500),
                approx_g = FALSE,
                nsave = 1000,
                ngrid = 100){

  # For testing:
  # locs = x; X = matrix(1, nrow = length(y)); covfun_name = "matern_isotropic"; locs_test = locs; X_test = X; nn = 30; emp_bayes = TRUE;fixedX = FALSE; approx_g = FALSE; nsave = 1000; ngrid = 100

  # Library required here:
  if (!requireNamespace("GpGp", quietly = TRUE)) {
    stop(
      "Package \"GpGp\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Library required here:
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop(
      "Package \"fields\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Data dimensions:
  y = as.matrix(y); n = length(y);
  locs = as.matrix(locs); d = ncol(locs)

  # Testing data:
  locs_test = as.matrix(locs_test); n_test = nrow(locs_test)

  # Covariates:
  if(is.null(X)) X = matrix(1, nrow = n)
  if(is.null(X_test)){ # supply our own testing matrix
    if(isTRUE(all.equal(locs, locs_test))){
      # If the training and testing points are the same,
      # then we input the same design matrix for the testing:
      X_test = X
    } else {
      # Otherwise, use an intercept-only with the correct dimensions:
      X_test = matrix(1, n_test)
    }
  }

  # And check:
  X = as.matrix(X); p = ncol(X)
  X_test = as.matrix(X_test)

  # Some checks needed for locs, locs_test, X, X_test
  if(nrow(locs) != n || nrow(X) != n || nrow(X_test) != n_test ||
     ncol(X_test) != p || ncol(locs_test) != d){
    stop('Check input dimensions!')
  }

  # To avoid errors for small n:
  nn = min(nn, n-1)

  # This is a temporary hack needed for sampling w/ one-dimensional inputs:
  if(!emp_bayes && d==1){
    aug = 1e-6*rnorm(n)
    locs = cbind(locs, aug)
    locs_test = cbind(locs_test, aug)
  }

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)
  #----------------------------------------------------------------------------
  print('Initial GP fit...')
  # Initial GP fit:
  z = qnorm(Fy(y))
  fit_gp = GpGp::fit_model(y = z,
                           locs = locs,
                           X = X,
                           covfun_name = covfun_name,
                           m_seq = nn,
                           silent = TRUE)

  # Fitted values for observed data:
  mu_z = GpGp::predictions(fit = fit_gp,
                            locs_pred = locs,
                            X_pred  = X)
  # SD of latent term:
  # Nugget variance (NOTE: this works for most (but not all!) covariance functions in GpGp)
  sigma_epsilon = sqrt(fit_gp$covparms[1]*fit_gp$covparms[length(fit_gp$covparms)])
  if(n < 1000){

    # Compute the covariance matrix, but remove the nugget:
    K_theta = do.call(covfun_name, list(fit_gp$covparms, locs ))
    #K_theta = match.fun(covfun_name)(fit_gp$covparms, locs)
    diag(K_theta) = diag(K_theta) - sigma_epsilon^2

    # Posterior covariance for mu requires inverses:
    K_theta_inv = chol2inv(chol(K_theta))
    Sigma_mu = chol2inv(chol(K_theta_inv + diag(1/sigma_epsilon^2, n)))

    # Extract the diagonal elements for z:
    sigma_z = sqrt(sigma_epsilon^2 + diag(Sigma_mu))

  } else {

    # Ignore the uncertainty in the regression function,
    # which is likely small when n is very large (also much faster...)
    sigma_z = rep(sigma_epsilon, n)

  }
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Evaluate CDF of Y at the unique y-values:
  y0 = sort(unique(y)); Fy_eval = Fy(y0)

  # Grid of values for the CDF of z:
  z_grid = sort(unique(
    sapply(range(mu_z), function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = ngrid),
            mean = xtemp,
            sd = sigma_epsilon)
    })
  ))

  # Evaluate the CDF of z:
  Fz_eval = Fz_fun(z = z_grid,
                   weights = rep(1/n, n),
                   mean_vec = mu_z,
                   sd_vec = sigma_z)

  # Compute the transformation:
  g = g_fun(y = y0, Fy_eval = Fy_eval,
            z = z_grid, Fz_eval = Fz_eval)

  # Latent data:
  z = g(y)

  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)
  #----------------------------------------------------------------------------
  print('Updated GP fit...')
  # Now update the GP coefficients:
  fit_gp = GpGp::fit_model(y = z,
                           locs = locs,
                           X = X,
                           covfun_name = covfun_name,
                           start_parms = fit_gp$covparms,
                           m_seq = nn,
                           silent = TRUE)

  # Fitted values for observed data:
  mu_z = GpGp::predictions(fit = fit_gp,
                           locs_pred = locs,
                           X_pred  = X,
                           m = nn)

  # Fitted values for testing data:
  if(isTRUE(all.equal(X, X_test)) &&
       isTRUE(all.equal(locs, locs_test))){
    # If the testing and training data are identical,
    # then there is no need to apply a separate predict function
    z_test = mu_z
  } else {
    z_test = GpGp::predictions(fit = fit_gp,
                               locs_pred = locs_test,
                               X_pred  = X_test,
                               m = nn)
  }

  # SD of latent term:
  # Nugget variance (NOTE: this works for most (but not all!) covariance functions in GpGp!)
  sigma_epsilon = sqrt(fit_gp$covparms[1]*fit_gp$covparms[length(fit_gp$covparms)])
  if(n < 1000){

    # Compute the covariance matrix, but remove the nugget:
    K_theta = do.call(covfun_name, list(fit_gp$covparms, locs ))
    #K_theta = match.fun(covfun_name)(fit_gp$covparms, locs)
    diag(K_theta) = diag(K_theta) - sigma_epsilon^2

    # Posterior covariance for mu requires inverses:
    K_theta_inv = chol2inv(chol(K_theta))
    Sigma_mu = chol2inv(chol(K_theta_inv + diag(1/sigma_epsilon^2, n)))

    # Extract the diagonal elements for z:
    sigma_z = sqrt(sigma_epsilon^2 + diag(Sigma_mu))

  } else {
    # Ignore the uncertainty in the regression function,
    # which is likely small when n is very large (this is also faster...)
    sigma_z = rep(sigma_epsilon, n)
  }

  # Estimated coefficients:
  theta = fit_gp$betahat
  #----------------------------------------------------------------------------
  # Store MC output:
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))

  print('Sampling...')

  # Run the MC:
  for(nsi in 1:nsave){
    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(!approx_g){

      # Bayesian bootstrap for the CDFs

      # BB CDF of y:
      Fy_eval = bb(y)(y0) # applied to y, evaluated at y0

      # BB CDF of z (only need if X is random)
      if(!fixedX){

        # Dirichlet(1) weights for x:
        weights_x = rgamma(n = n, shape = 1)
        weights_x  = weights_x/sum(weights_x)

        # CDF of Z evaluated on the grid:
        Fz_eval = Fz_fun(z = z_grid,
                         weights = weights_x,
                         mean_vec = mu_z,
                         sd_vec = sigma_z)
      } # otherwise, the weights are 1/n (how we initialized)

      # Compute the transformation:
      g = g_fun(y = y0, Fy_eval = Fy_eval,
                z = z_grid, Fz_eval = Fz_eval)

      # Update z:
      z = g(y)

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = y_grid)
    }
    #----------------------------------------------------------------------------
    # w/o empirical Bayes approximation:
    #   update fitted curve (based on z)
    #   & sample sigma_epsilon
    if(!emp_bayes){
      # Block 2: update the fitted GP
      fit_gp = GpGp::fit_model(y = z,
                               locs = locs, X = X, covfun_name = covfun_name, start_parms = fit_gp$covparms, m_seq = nn, silent = TRUE)
      f_hat = GpGp::predictions(fit = fit_gp,
                                locs_pred = locs,
                                X_pred = X,
                                m = nn)

      # Block 3: sample the scale adjustment (SD)
      sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                    shape = .001 + n/2,
                                    rate = .001 + sum((z - f_hat)^2)/2))
    }
    #----------------------------------------------------------------------------
    # Store the MC:

    # Predictive samples of ytilde:
    if(emp_bayes){
      ztilde = z_test + sigma_epsilon*rnorm(n = n_test)
    } else {
      ztilde = GpGp::cond_sim(fit = fit_gp,
                              locs_pred = locs_test,
                              X_pred = X_test,
                              m = nn)
    }
    post_ypred[nsi,] = g_inv(ztilde)

    # Posterior samples of the transformation:
    post_g[nsi,] = (g(y0) - theta[1])/sigma_epsilon # undo location/scale (not identified)
    #----------------------------------------------------------------------------
  }
  print('Done!')

  return(list(
    coefficients = theta,
    fitted.values = colMeans(post_ypred),
    fit_gp = fit_gp,
    post_ypred = post_ypred,
    post_g = post_g,
    model = 'sbgp', y = y, X = X, X_test = X_test, nn = nn, emp_bayes = emp_bayes, fixedX = fixedX, approx_g = approx_g))
}
#' Semiparametric Bayesian quantile regression
#'
#' MCMC sampling for Bayesian quantile regression with an
#' unknown (nonparametric) transformation. Like in traditional Bayesian
#' quantile regression, an asymmetric Laplace distribution is assumed
#' for the errors, so the regression models targets the specified quantile.
#' However, these models are often woefully inadequate for describing
#' observed data. We introduce a nonparametric transformation to
#' improve model adequacy while still providing inference for the
#' regression coefficients and the specified quantile. A g-prior is assumed
#' for the regression coefficients.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors (no intercept)
#' @param tau the target quantile (between zero and one)
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param psi prior variance (g-prior)
#' @param laplace_approx logical; if TRUE, use a normal approximation
#' to the posterior in the definition of the transformation;
#' otherwise the prior is used
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param approx_g logical; if TRUE, apply large-sample
#' approximation for the transformation
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param ngrid number of grid points for inverse approximations
#' @param verbose logical; if TRUE, print time remaining
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{fitted.values} the estimated \code{tau}th quantile at test points \code{X_test}
#' \item \code{post_theta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_ypred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_qtau}: \code{nsave x n_test} samples of the \code{tau}th conditional quantile at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{model}: the model fit (here, \code{sbqr})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed quantile linear model.
#' The transformation is modeled as unknown and learned jointly
#' with the regression coefficients (unless \code{approx_g} = TRUE, which then uses
#' a point approximation). This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' The results are typically unchanged whether \code{laplace_approx} is TRUE/FALSE;
#' setting it to TRUE may reduce sensitivity to the prior, while setting it to FALSE
#' may speed up computations for very large datasets. Similarly, treating the
#' covariates as fixed (\code{fixedX = TRUE}) can substantially improve
#' computing efficiency, so we make this the default.
#'
#' @note The location (intercept) is not identified, so any intercepts
#' in \code{X} and \code{X_test} will be removed. The model-fitting *does*
#' include an internal location-scale adjustment, but the function only outputs
#' inferential summaries for the identifiable parameters.
#'
#' @examples
#' \donttest{
#' # Simulate some heteroskedastic data (no transformation):
#' dat = simulate_tlm(n = 200, p = 10, g_type = 'box-cox', heterosked = TRUE, lambda = 1)
#' y = dat$y; X = dat$X # training data
#' y_test = dat$y_test; X_test = dat$X_test # testing data
#'
#' # Target this quantile:
#' tau = 0.05
#'
#' # Fit the semiparametric Bayesian quantile regression model:
#' fit = sbqr(y = y, X = X, tau = tau, X_test = X_test)
#' names(fit) # what is returned
#'
#' # Posterior predictive checks on testing data: empirical CDF
#' y0 = sort(unique(y_test))
#' plot(y0, y0, type='n', ylim = c(0,1),
#'      xlab='y', ylab='F_y', main = 'Posterior predictive ECDF')
#' temp = sapply(1:nrow(fit$post_ypred), function(s)
#'   lines(y0, ecdf(fit$post_ypred[s,])(y0), # ECDF of posterior predictive draws
#'         col='gray', type ='s'))
#' lines(y0, ecdf(y_test)(y0),  # ECDF of testing data
#'      col='black', type = 's', lwd = 3)
#' }
#'
# #' @importFrom quantreg rq
# #' @importFrom statmod rinvgauss
#' @export
sbqr = function(y, X, tau = 0.5,
                X_test = X,
                psi = length(y),
                laplace_approx = TRUE,
                fixedX = TRUE,
                approx_g = FALSE,
                nsave = 1000,
                nburn = 100,
                ngrid = 100,
                verbose = TRUE){

  # For testing:
  # X_test = X; tau = 0.10; psi = length(y); laplace_approx = TRUE; fixedX=FALSE; approx_g = FALSE; nsave = 1000; nburn = 100; verbose = TRUE; ngrid = 100

  # Library required here:
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop(
      "Package \"quantreg\" must be installed to use this function.",
      call. = FALSE
    )
  }
  # Library required here:
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop(
      "Package \"statmod\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Initial checks:
  if(!is.matrix(X)) stop("X must be a matrix (rows = observations, columns = variables)")
  if(!is.matrix(X_test)) stop("X_test must be a matrix (rows = observations, columns = variables)")
  if(ncol(X) != ncol(X_test)) stop('X_test and X must have the same number of columns (variables)')
  if(nrow(X) != length(y)) stop('the length of y must equal nrow(X)')

  # Here: intercept is 1st column of X, X_test
  #   (this is not part of the target parameter theta)
  is_intercept = apply(X, 2, function(x) all(x == 1))
  X = cbind(1, X[,!is_intercept])
  is_intercept_test = apply(X_test, 2, function(x) all(x == 1))
  X_test = cbind(1, X_test[,!is_intercept_test])

  # Data dimensions:
  n = length(y) # number of observations
  p = ncol(X) - 1 # number of variables (excluding intercept)
  n_test = nrow(X_test) # number of testing data points

  # Requirement for the g-prior:
  if(p >= n) stop('The g-prior requires p < n')
  #----------------------------------------------------------------------------
  # Recurring terms:
  a_tau = (1-2*tau)/(tau*(1-tau))
  b_tau = sqrt(2/(tau*(1-tau)))

  # Key matrix quantities:
  XtX = crossprod(X)
  XtXinv = chol2inv(chol(XtX))
  xt_Sigma_x = sapply(1:n, function(i)
    crossprod(X[i,], XtXinv)%*%X[i,])
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)

  # Evaluate at the unique y-values:
  y0 = sort(unique(y))
  Fy_eval = Fy(y0)

  # Grid of values for the CDF of z:
  z_grid = sort(unique(
    sapply(range(psi*xt_Sigma_x), function(xtemp){
      qnorm(seq(0.001, 0.999, length.out = ngrid),
            mean = 0 + a_tau, # assuming prior mean zero
            sd = sqrt(b_tau^2 + xtemp))
    })
  ))

  # Define the moments of the CDF of z (BEFORE parameter expansions)
  if(laplace_approx){
    # Use a normal approximation for the posterior of theta

    # First pass: fix Fz() = qnorm(), initialize coefficients
    z = qnorm(Fy(y))
    fit = quantreg::rq(z ~ X - 1, tau = tau)
    Sigma_hat = summary(fit, se = 'boot', covariance = TRUE)$cov
    xt_Sigma_hat_x = sapply(1:n, function(i)
      crossprod(X[i,], Sigma_hat)%*%X[i,])

    # Moments of X%*%theta:
    mu_z = fitted(fit)
    sigma_z = sqrt(xt_Sigma_hat_x)

    # CDF of z, using Monte Carlo:
    Fz_eval = rowMeans(sapply(rexp(n = 100), function(xi_s){
      Fz_fun(z = z_grid,
             weights = rep(1/n, n),
             mean_vec = mu_z + a_tau*xi_s,
             sd_vec = sqrt(sigma_z^2 + xi_s*b_tau^2))
    }))

    # Check: update the grid if needed
    zcon = contract_grid(z = z_grid,
                         Fz = Fz_eval,
                         lower = 0.001, upper =  0.999)
    z_grid = zcon$z; Fz_eval = zcon$Fz

    # Transformation:
    g = g_fun(y = y0,
              Fy_eval = Fy_eval,
              z = z_grid,
              Fz_eval = Fz_eval)

    # Updated coefficients:
    z = g(y) # update latent data
    fit = quantreg::rq(z ~ X - 1, tau = tau)
    Sigma_hat = summary(fit, se = 'boot', covariance = TRUE)$cov
    xt_Sigma_hat_x = sapply(1:n, function(i)
      crossprod(X[i,], Sigma_hat)%*%X[i,])

    # Moments of X%*%theta:
    mu_z = fitted(fit)
    sigma_z = sqrt(xt_Sigma_hat_x)

  } else {

    # Using the prior
    mu_z = rep(0, n)
    sigma_z = sqrt(psi*xt_Sigma_x)
  }

  # CDF of z, using Monte Carlo:
  Fz_eval = rowMeans(sapply(rexp(n = 100), function(xi_s){
    Fz_fun(z = z_grid,
           weights = rep(1/n, n),
           mean_vec = mu_z + a_tau*xi_s,
           sd_vec = sqrt(sigma_z^2 + xi_s*b_tau^2))
  }))

  # Check: update the grid if needed
  zcon = contract_grid(z = z_grid,
                       Fz = Fz_eval,
                       lower = 0.001, upper =  0.999)
  z_grid = zcon$z; Fz_eval = zcon$Fz

  # Compute the transformation:
  g = g_fun(y = y0, Fy_eval = Fy_eval,
            z = z_grid, Fz_eval = Fz_eval)

  # Latent data:
  z = g(y)

  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)
  #----------------------------------------------------------------------------
  # Initialize the coefficients:
  theta = chol2inv(chol(XtX))%*%crossprod(X,z) # initial coefs
  #----------------------------------------------------------------------------
  # Store MC output:
  post_theta = array(NA, c(nsave, p))
  post_ypred = post_qtau = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:(nburn + nsave)){

    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(!approx_g){

      # Bayesian bootstrap for the CDFs

      # BB CDF of y:
      Fy_eval = bb(y)(y0) # applied to y, evaluated at y0

      # BB CDF of z (only need if X is random)
      if(!fixedX){

        # Dirichlet(1) weights for x:
        weights_x = rgamma(n = n, shape = 1)
        weights_x  = weights_x/sum(weights_x)

        # CDF of Z evaluated on the grid
        #   (MC averages are good even for few replicates!)
        Fz_eval = rowMeans(sapply(rexp(n = 10), function(xi_s){
          Fz_fun(z = z_grid,
                 weights = weights_x,
                 mean_vec = mu_z + a_tau*xi_s,
                 sd_vec = sqrt(sigma_z^2 + xi_s*b_tau^2))
        }))
      } # otherwise, the weights are 1/n (how we initialized)

      # Compute the transformation:
      g = g_fun(y = y0, Fy_eval = Fy_eval,
                z = z_grid, Fz_eval = Fz_eval)

      # Update z:
      z = g(y)

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = y_grid)
    }
    #----------------------------------------------------------------------------
    # Block 2: parameter expansion
    xi = 1/statmod::rinvgauss(n = n,
                              mean = sqrt((2 + a_tau^2/b_tau^2)/((z - X%*%theta)^2/b_tau^2)),
                              shape = 2 + a_tau^2/b_tau^2)
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    Q_theta = crossprod(X/sqrt(b_tau^2*xi)) + 1/psi*XtX # t(X)%*%diag(1/(b_tau^2*xi))%*%X + 1/psi*XtX
    ell_theta = crossprod(X/(b_tau^2*xi), z - a_tau*xi) # t(X)%*%diag(1/(b_tau^2*xi))%*%(z - a_tau*xi)
    ch_Q = chol(Q_theta)
    theta = backsolve(ch_Q,
                      forwardsolve(t(ch_Q), ell_theta) +
                        rnorm(p + 1))
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){
      # Posterior samples of the model parameters:
      post_theta[nsi - nburn,] = theta[-1]

      # Quantile at the testing point:
      post_qtau[nsi - nburn,] = g_inv(X_test%*%theta)

      # Predictive samples of ytilde:
      xi_test = rexp(n = n_test, rate = 1)
      ztilde = X_test%*%theta + a_tau*xi_test + b_tau*sqrt(xi_test)*rnorm(n = n_test)
      post_ypred[nsi - nburn,] = g_inv(ztilde)

      # Posterior samples of the transformation:
      post_g[nsi - nburn,] = g(y0) - theta[1] # undo location (not identified; no scale here)
    }
    #----------------------------------------------------------------------------
    if(verbose) computeTimeRemaining(nsi, timer0, nsave)
  }
  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer0
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  return(list(
    coefficients = colMeans(post_theta),
    fitted.values = colMeans(post_qtau),
    post_theta = post_theta,
    post_ypred = post_ypred,
    post_qtau = post_qtau,
    post_g = post_g,
    model = 'sbqr', y = y, X = X[,-1], X_test = X_test[,-1], psi = psi, laplace_approx = laplace_approx, fixedX = fixedX, approx_g = approx_g, tau = tau))
}
#---------------------------------------------------------------
#' Post-processing with importance sampling
#'
#' Given Monte Carlo draws from the surrogate posterior,
#' apply sampling importance reweighting (SIR) to
#' correct for the true model likelihood.
#'
#' @param fit a fitted model object that includes
#' \itemize{
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{post_theta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_ypred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{model}: the model fit (\code{sblm} or \code{sbsm})
#' }
#' @param sir_frac fraction of draws to sample for SIR
#' @param nsims_prior number of draws from the prior
#' @param marg_x logical; if TRUE, compute the weights marginal over the covariates
#' @param verbose logical; if TRUE, print time remaining
#' @return the fitted model object with the posterior draws subsampled
#' based on the SIR adjustment
#'
#' @details The Monte Carlo sampling for \code{\link{sblm}} and
#' \code{\link{sbsm}} uses a surrogate likelihood for posterior inference,
#' which enables much faster and easier computing. SIR provides a correction for
#' the actual (specified) likelihood. However, this correction
#' step typically does not produce any noticeable
#' discrepancies, even for small sample sizes.
#'
#' @note SIR sampling is done *without* replacement, so \code{sir_frac}
#' is typically between 0.1 and 0.5. The \code{nsims_priors} draws
#' are used to approximate a prior expectation, but larger values
#' can significantly slow down this function. The importance weights
#' can be computed conditionally (\code{marg_x = FALSE}) or unconditionally
#' (\code{marg_x = TRUE}) on the covariates, corresponding to whether
#' or not the covariates are marginalized out in the likelihood. The
#' conditional version is much faster.
#'
#' @examples
#' \donttest{
#' # Simulate some data:
#' dat = simulate_tlm(n = 50, p = 5, g_type = 'step')
#' y = dat$y; X = dat$X # training data
#' y_test = dat$y_test; X_test = dat$X_test # testing data
#'
#' hist(y, breaks = 10) # marginal distribution
#'
#' # Fit the semiparametric Bayesian linear model:
#' fit = sblm(y = y, X = X, X_test = X_test)
#' names(fit) # what is returned
#'
#' # Update with SIR:
#' fit_sir = sir_adjust(fit)
#' names(fit_sir) # what is returned
#'
#' # Prediction: unadjusted vs. adjusted?
#'
#' # Point estimates:
#' y_hat = fitted(fit)
#' y_hat_sir = fitted(fit_sir)
#' cor(y_hat, y_hat_sir) # similar
#'
#' # Interval estimates:
#' pi_y = t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
#' pi_y_sir = t(apply(fit_sir$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
#'
#' # PI overlap (%):
#' overlaps = 100*sapply(1:length(y_test), function(i){
#'   # innermost part
#'   (min(pi_y[i,2], pi_y_sir[i,2]) - max(pi_y[i,1], pi_y_sir[i,1]))/
#'     # outermost part
#'     (max(pi_y[i,2], pi_y_sir[i,2]) - min(pi_y[i,1], pi_y_sir[i,1]))
#' })
#' summary(overlaps) # mostly close to 100%
#'
#' # Coverage of PIs on testing data (should be ~ 90%)
#' mean((pi_y[,1] <= y_test)*(pi_y[,2] >= y_test)) # unadjusted
#' mean((pi_y_sir[,1] <= y_test)*(pi_y_sir[,2] >= y_test)) # adjusted
#'
#' # Plot together with testing data:
#' plot(y_test, y_test, type='n', ylim = range(pi_y, pi_y_sir, y_test),
#'      xlab = 'y_test', ylab = 'y_hat', main = paste('Prediction intervals: testing data'))
#' abline(0,1) # reference line
#' suppressWarnings(
#'   arrows(y_test, pi_y[,1], y_test, pi_y[,2],
#'          length=0.15, angle=90, code=3, col='gray', lwd=2)
#' ) # plot the PIs (unadjusted)
#' suppressWarnings(
#'   arrows(y_test, pi_y_sir[,1], y_test, pi_y_sir[,2],
#'          length=0.15, angle=90, code=3, col='darkgray', lwd=2)
#' ) # plot the PIs (adjusted)
#' lines(y_test, y_hat, type='p', pch=2) # plot the means (unadjusted)
#' lines(y_test, y_hat_sir, type='p', pch=3) # plot the means (adjusted)
#' }
# #' @importFrom MASS mvrnorm
#' @export
sir_adjust = function(fit,
                      sir_frac = 0.3,
                      nsims_prior = 100,
                      marg_x = FALSE,
                      verbose = TRUE){

  # Checks:
  if(fit$model == 'sbgp'){
    warning('sbgp uses an empirical bayes approximation for the parameters,
    so SIR sampling is not needed; returning the original fit...')
    return(fit)
  }
  if(fit$model != 'sblm' &&
     fit$model != 'sbsm')
    stop('Currently implemented for sblm and sbsm only')

  if(sir_frac <= 0 || sir_frac >= 1){
    stop('sir_frac must be between zero and one')
  }

  # Extract some recurring terms:
  nsave = nrow(fit$post_ypred) # number of draws
  y = fit$y; n = length(y);
  X = fit$X; p = ncol(X)

  # Useful for matching:
  ind_y = match(y, sort(unique(y)))

  # Not identified:
  sigma_epsilon = 1 # fit$sigma_epsilon

  # Sample from the prior:

  # Linear model case:
  if(fit$model == 'sblm'){

    # Library required here:
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop(
        "Package \"MASS\" must be installed to use this function.",
        call. = FALSE
      )
    }

    # Matrix quantities (note: these *could* be passed in from fit)
    XtXinv = chol2inv(chol(crossprod(X)))
    xt_Sigma_x = sapply(1:n, function(i)
      crossprod(X[i,], XtXinv)%*%X[i,])

    # Sample:
    prior_theta = MASS::mvrnorm(n = nsims_prior,
                                mu = rep(0,p),
                                Sigma = sigma_epsilon^2*fit$psi*XtXinv)
  }

  # Spline case:
  if(fit$model == 'sbsm'){

    # Matrix quantities (note: these *could* be passed in from fit)
    xt_Sigma_x = rowSums(X^2)

    # Sample:
    prior_theta = matrix(rnorm(n = p*nsims_prior,
                               mean = 0,
                               sd = sigma_epsilon*sqrt(fit$psi)),
                         nrow = nsims_prior)
  }

  # Compute the log-importance weights:
  post_logw = rep(NA, nsave)
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nsave){

    # Extract g(y), properly matched:
    z = fit$post_g[nsi, ind_y]

    # Dirichlet(1) weights for x:
    weights_x = rgamma(n = n, shape = 1)
    weights_x  = weights_x/sum(weights_x)

    # The numerator can be unconditional or conditional on x
    if(marg_x){
      # Log-numerator, marginal over x:
      log_numer = log(mean(sapply(1:nsims_prior, function(s){
        prod(sapply(1:n, function(i){
          sum(weights_x*dnorm(z[i],
                              mean = X%*%prior_theta[s,],
                              sd = sigma_epsilon))
        }))
      })))
    } else {
      # Log-numerator, conditional on x:
      # log_numer = log(mean(sapply(1:nsims_prior, function(s){
      #   exp(sum(dnorm(z,
      #                 mean = X%*%prior_theta[s,],
      #                 sd = sigma_epsilon,
      #                 log=TRUE)))
      # })))
      # For numerical stability...
      #   this needs to be log-transformed before averaging
      #   but we can rescale for convenience...
      temp = sapply(1:nsims_prior, function(s){
        sum(dnorm(z,
                  mean = X%*%prior_theta[s,],
                  sd = sigma_epsilon,
                  log=TRUE))
      })
      log_numer = log(mean(exp(temp - max(temp)))) # normalize, then exponentiate, then average, then log
    }

    # Log-denominator
    log_denom = sum(log(
      sapply(1:n, function(i){
        sum(weights_x*dnorm(z[i],
                            mean = 0,
                            sd = sigma_epsilon*sqrt(1 + fit$psi*xt_Sigma_x)))
      })
    ))

    # Store the log-weights:
    post_logw[nsi] = log_numer - log_denom

    # Check the timing:
    if(verbose) computeTimeRemaining(nsi, timer0, nsave)
  }
  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer0
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  # Indices to keep:
  ind_sir = sample(1:nsave,
                   size = ceiling(sir_frac*nsave),
                   prob = exp(post_logw - max(post_logw)), # normalize
                   replace = FALSE)
  # Update:
  fit$post_theta = fit$post_theta[ind_sir,];
  fit$coefficients = colMeans(fit$post_theta)
  fit$post_ypred = fit$post_ypred[ind_sir,]
  fit$fitted.values = colMeans(fit$post_ypred)
  fit$post_g = fit$post_g[ind_sir,]

  return(fit)
}
#---------------------------------------------------------------
#' Compute the transformation
#'
#' Given the CDFs of z and y, compute a smoothed function
#' to evaluate the transformation
#'
#' @param y vector of points at which the CDF of y is evaluated
#' @param Fy_eval CDF of y evaluated at \code{y}
#' @param z vector of points at which the CDF of z is evaluated
#' @param Fz_eval CDF of z evaluated at \code{z}
#' @return A smooth monotone function which can be used for evaluations of the transformation.
#' @importFrom stats splinefun
g_fun = function(y, Fy_eval, z, Fz_eval){

  # Quick checks:
  if(length(y) != length(Fy_eval)){
    stop('length of y must equal length of Fy_eval')
  }
  if(length(z) != length(Fz_eval)){
    stop('length of z must equal length of Fz_eval')
  }

  # Remove duplicates:
  z_unique = which(!duplicated(Fz_eval))
  Fz_eval = Fz_eval[z_unique]; z = z[z_unique]
  y_unique = which(!duplicated(Fy_eval))
  Fy_eval = Fy_eval[y_unique]; y = y[y_unique]

  # Inverse the CDF of z:
  Fz_inv = function(s) stats::spline(Fz_eval, z, method = "hyman",
                                     xout = s)$y
  # Ref: https://stats.stackexchange.com/questions/390931/compute-quantile-function-from-a-mixture-of-normal-distribution/390936#390936

  # Check:
  # plot(z, Fz_inv(Fz_eval)); abline(0,1)

  # Initial transformation:
  g0 = Fz_inv(Fy_eval)

  # Make sure we have only finite values of g0 (infinite values occur for F_y = 0 or F_y = 1)
  y = y[which(is.finite(g0))]; g0 = g0[which(is.finite(g0))]

  # Return the smoothed (monotone) transformation:
  return(stats::splinefun(y, g0, method = 'monoH.FC'))
}
#----------------------------------------------------------------------------
#' Approximate inverse transformation
#'
#' Compute the inverse function of a transformation \code{g} based on a grid search.
#'
#' @param g the transformation function
#' @param t_grid grid of arguments at which to evaluate the transformation function
#' @return A function which can be used for evaluations of the
#' (approximate) inverse transformation function.
g_inv_approx = function(g, t_grid) {

  # Evaluate g() on the grid:
  g_grid = g(t_grid)

  # Approximate inverse function:
  function(s) {
    sapply(s, function(si)
      t_grid[which.min(abs(si - g_grid))])
  }
}
#---------------------------------------------------------------
#' Compute the latent data CDF
#'
#' Assuming a Gaussian latent data distribution (given x),
#' compute the CDF on a grid of points
#'
#' @param z vector of points at which the CDF of z is evaluated
#' @param weights \code{n}-dimensional vector of weights; if NULL,
#' assume 1/n
#' @param mean_vec \code{n}-dimensional vector of means; if NULL,
#' assume mean zero
#' @param sd_vec \code{n}-dimensional vector of standard deviations
#' @return CDF of z evaluated at \code{z}
Fz_fun = function(z,
                  weights = NULL,
                  mean_vec = NULL,
                  sd_vec){
  # Dimension:
  n = length(sd_vec)

  # Checks:
  if(is.null(weights)) weights = rep(1/n, n)
  if(is.null(mean_vec)) mean_vec = rep(0, n)

  # Return:
  rowSums(sapply(1:n, function(i){
    weights[i]*pnorm(z,
                mean = mean_vec[i],
                sd = sd_vec[i])
  }))
}
#' Grid contraction
#'
#' Contract the grid if the evaluation points exceed some threshold.
#' This removes the corresponding z values.
#' We can add points back to achieve the same (approximate) length.
#'
#' @param z grid points (ordered)
#' @param Fz function evaluated at those grid points
#' @param lower lower threshold at which to check Fz
#' @param upper upper threshold at which to check Fz
#' @param add_back logical; if true, expand the grid to
#' (about) the original size
#' @param monotone logical; if true, enforce monotonicity
#' on the expanded grid
#' @return a list containing the grid points \code{z} and the (interpolated) function
#' \code{Fz} at those points
contract_grid = function(z, Fz, lower, upper, add_back = TRUE, monotone = TRUE){

  # Identify the offending indexes:
  bad_ind = which(Fz < lower | Fz > upper)

  if(length(bad_ind) > 0){

    # Original:
    z0 = z

    # Original length
    ngrid = length(z)

    # Subset:
    z = z[-bad_ind]

    # Add back, if desired:
    if(add_back){
      z = sort(unique(c(z,
                        seq(min(z),
                            max(z),
                            length.out = ngrid - length(z) + 2))))
    }
    if(monotone){
      Fz = stats::spline(x = z0[-bad_ind], y = Fz[-bad_ind], method = 'hyman', xout = z)$y
    } else {
      Fz = stats::spline(x = z0[-bad_ind], y = Fz[-bad_ind], xout = z)$y
    }
  }
  list(z = z,
       Fz = Fz)
}
