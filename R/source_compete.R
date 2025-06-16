#---------------------------------------------------------------
# Source functions for *competing methods* for semiparametric Bayesian analysis
#---------------------------------------------------------------
#' Bayesian linear model with a Box-Cox transformation
#'
#' MCMC sampling for Bayesian linear regression with a
#' (known or unknown) Box-Cox transformation. A g-prior is assumed
#' for the regression coefficients.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors (no intercept)
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param psi prior variance (g-prior)
#' @param lambda Box-Cox transformation; if NULL, estimate this parameter
#' @param sample_lambda logical; if TRUE, sample lambda, otherwise
#' use the fixed value of lambda above or the MLE (if lambda unspecified)
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' \item \code{post_lambda} \code{nsave} posterior samples of lambda
#' \item \code{post_sigma} \code{nsave} posterior samples of sigma
#' \item \code{model}: the model fit (here, \code{blm_bc})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed linear model via MCMC sampling. The transformation is
#' parametric from the Box-Cox family, which has one parameter \code{lambda}.
#' That parameter may be fixed in advanced or learned from the data.
#'
#' @note Box-Cox transformations may be useful in some cases, but
#' in general we recommend the nonparametric transformation (with
#' Monte Carlo, not MCMC sampling) in \code{\link{sblm}}.
#'
#' @note An intercept is automatically added to \code{X} and
#' \code{X_test}. The coefficients reported do *not* include
#' this intercept parameter, since it is not identified
#' under more general transformation models (e.g., \code{\link{sblm}}).
#'
#' @examples
#' # Simulate some data:
#' dat = simulate_tlm(n = 100, p = 5, g_type = 'step')
#' y = dat$y; X = dat$X # training data
#' y_test = dat$y_test; X_test = dat$X_test # testing data
#'
#' hist(y, breaks = 25) # marginal distribution
#'
#' # Fit the Bayesian linear model with a Box-Cox transformation:
#' fit = blm_bc(y = y, X = X, X_test = X_test)
#' names(fit) # what is returned
#' round(quantile(fit$post_lambda), 3) # summary of unknown Box-Cox parameter
#'
#' @export
blm_bc = function(y, X, X_test = X,
                   psi = length(y),
                   lambda = NULL,
                   sample_lambda = TRUE,
                   nsave = 1000,
                   nburn = 1000,
                   nskip = 0,
                   verbose = TRUE){

  # For testing:
  # X_test = X; psi = length(y); sample_lambda  = FALSE; lambda = NULL; nsave = 1000; nburn = 1000; nskip = 0; verbose = TRUE

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
  y0 = sort(unique(y))
  XtX = crossprod(X)
  XtXinv = chol2inv(chol(XtX))
  #----------------------------------------------------------------------------
  # Initialize the parameters:

  # The scale is not well-identified for unknown lambda, so initialize at one
  sigma_epsilon = 1

  # Optimize the Box-Cox parameter:
  if(is.null(lambda)){

    # Compute the hat matrix:
    XtXinv = chol2inv(chol(XtX))
    XtXinvXt = tcrossprod(XtXinv, X)
    H = X%*%XtXinvXt # hat matrix

    # Marginal precision for z:
    Sigma_z_inv  = chol2inv(chol(
      sigma_epsilon^2*(diag(n) + psi*H)
    ))

    # MLE:
    opt = optim(par = 0.5,
                fn = function(l_bc){
                  z = g_bc(y, lambda = l_bc)
                  crossprod(z, Sigma_z_inv)%*%(z)
                }, method  = "L-BFGS-B", lower = 0.01, upper = 1
    )
    if(opt$convergence  ==  0){
      # Successful optimization
      lambda = opt$par
    } else {
      warning('Optimization failed...setting lambda = 1/2')
      lambda = 1/2
    }
  }

  # Latent data:
  z = g_bc(y, lambda = lambda)
  Xtz = crossprod(X, z) # only changes  if sample_lambda

  # Coefficients:
  theta = chol2inv(chol(XtX))%*%Xtz
  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_theta = array(NA, c(nsave, p))
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))
  post_lambda = post_sigma = rep(NA, nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(sample_lambda){
      # Sample lambda:
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           # Likelihood for lambda:
                           sum(dnorm(g_bc(y, lambda = l_bc),
                                     mean = X%*%theta,
                                     sd = sigma_epsilon, log = TRUE)) +
                             # This is the prior on lambda, truncated to [0, 2]
                             dnorm(l_bc, mean = 1/2, sd = 1/2, log = TRUE)
                         },
                         w = 1/2, m = 50, lower = 0, upper = 2)

      # Update z:
      z = g_bc(y, lambda = lambda)
      Xtz = crossprod(X, z)

    }
    #----------------------------------------------------------------------------
    # Block 2: sample the error SD
    SSR_psi = sum(z^2) - psi/(psi+1)*crossprod(Xtz, XtXinv%*%Xtz)
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = .001 + n/2,
                                  rate = .001 + SSR_psi/2))
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    ch_Q = chol(1/sigma_epsilon^2*(1+psi)/(psi)*XtX)
    ell_theta = 1/sigma_epsilon^2*Xtz
    theta = backsolve(ch_Q,
                      forwardsolve(t(ch_Q), ell_theta) +
                        rnorm(p+1))
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Posterior samples of the model parameters:
        post_theta[isave,] = theta[-1]

        # Predictive samples of ytilde:
        ztilde = X_test%*%theta + sigma_epsilon*rnorm(n = n_test)
        post_ypred[isave,] = g_inv_bc(ztilde, lambda = lambda)

        # Posterior samples of the transformation:
        post_g[isave,] = g_bc(y0, lambda = lambda)
        post_lambda[isave] = lambda

        # Posterior samples of the error SD:
        post_sigma[isave] = sigma_epsilon

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
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
    post_g = post_g, post_lambda = post_lambda, post_sigma = post_sigma,
    model = 'blm_bc', y = y, X = X[,-1], X_test = X_test[,-1], psi = psi))
}
#' Bayesian spline model with a Box-Cox transformation
#'
#' MCMC sampling for Bayesian spline regression with a
#' (known or unknown) Box-Cox transformation.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param x \code{n x 1} vector of observation points; if NULL, assume equally-spaced on [0,1]
#' @param x_test \code{n_test x 1} vector of testing points; if NULL, assume equal to \code{x}
#' @param psi prior variance (inverse smoothing parameter); if NULL,
#' sample this parameter
#' @param lambda Box-Cox transformation; if NULL, estimate this parameter
#' @param sample_lambda logical; if TRUE, sample lambda, otherwise
#' use the fixed value of lambda above or the MLE (if lambda unspecified)
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' \item \code{post_lambda} \code{nsave} posterior samples of lambda
#' \item \code{model}: the model fit (here, \code{sbsm_bc})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed spline model via MCMC sampling. The transformation is
#' parametric from the Box-Cox family, which has one parameter \code{lambda}.
#' That parameter may be fixed in advanced or learned from the data.
#'
#' @note Box-Cox transformations may be useful in some cases, but
#' in general we recommend the nonparametric transformation (with
#' Monte Carlo, not MCMC sampling) in \code{\link{sbsm}}.
#'
#' @examples
#' # Simulate some data:
#' n = 100 # sample size
#' x = sort(runif(n)) # observation points
#'
#' # Transform a noisy, periodic function:
#' y = g_inv_bc(
#'   sin(2*pi*x) + sin(4*pi*x) + rnorm(n, sd = .5),
#'              lambda = .5) # Signed square-root transformation
#'
#' # Fit the Bayesian spline model with a Box-Cox transformation:
#' fit = bsm_bc(y = y, x = x)
#' names(fit) # what is returned
#' round(quantile(fit$post_lambda), 3) # summary of unknown Box-Cox parameter
#'
#' # Plot the model predictions (point and interval estimates):
#' pi_y = t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
#' plot(x, y, type='n', ylim = range(pi_y,y),
#'      xlab = 'x', ylab = 'y', main = paste('Fitted values and prediction intervals'))
#' polygon(c(x, rev(x)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
#' lines(x, y, type='p')
#' lines(x, fitted(fit), lwd = 3)
#'
# #' @importFrom spikeSlabGAM sm
#' @export
bsm_bc = function(y, x = NULL,
                   x_test = x,
                   psi = NULL,
                   lambda = NULL,
                   sample_lambda = TRUE,
                   nsave = 1000,
                   nburn = 1000,
                   nskip = 0,
                   verbose = TRUE){

  # For testing:
  # psi = length(y); sample_lambda  = FALSE; lambda = NULL; nsave = 1000; nburn = 1000; nskip = 0; verbose = TRUE

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
  #----------------------------------------------------------------------------
  # Orthogonalized P-spline and related quantities:
  X = cbind(1/sqrt(n), poly(x, 1), spikeSlabGAM::sm(x))
  X = X/sqrt(sum(diag(crossprod(X))))
  diagXtX = colSums(X^2)
  p = length(diagXtX)

  if(is.null(psi)){
    sample_psi = TRUE
    psi = n # initialized
  } else sample_psi = FALSE

  # Recurring terms:
  y0 = sort(unique(y))
  #----------------------------------------------------------------------------
  # Initialize the parameters:

  # The scale is not well-identified for unknown lambda:
  sigma_epsilon = 1

  # Optimize the Box-Cox parameter:
  if(is.null(lambda)){

    # Key term:
    XXt = tcrossprod(X)

    # Marginal precision for z:
    Sigma_z_inv  = chol2inv(chol(
      sigma_epsilon^2*(diag(n) + psi*XXt)
    ))

    # MLE:
    opt = optim(par = 0.5,
                fn = function(l_bc){
                  z = g_bc(y, lambda = l_bc)
                  crossprod(z, Sigma_z_inv)%*%(z)
                }, method  = "L-BFGS-B", lower = 0.01, upper = 1
    )
    if(opt$convergence  ==  0){
      # Successful optimization
      lambda = opt$par
    } else {
      warning('Optimization failed...setting lambda = 1/2')
      lambda = 1/2
    }
  }

  # Latent data:
  z = g_bc(y, lambda = lambda)
  Xtz = crossprod(X, z) # only changes  if sample_lambda

  # Coefficients:
  Q_theta = 1/sigma_epsilon^2*diagXtX + 1/(sigma_epsilon^2*psi)
  ell_theta = 1/sigma_epsilon^2*Xtz # only changes if g changes
  theta = Q_theta^-1*ell_theta
  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_theta = array(NA, c(nsave, p))
  post_ypred = array(NA, c(nsave, length(x_test)))
  post_g = array(NA, c(nsave, length(y0)))
  post_lambda = rep(NA, nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(sample_lambda){
      # Sample lambda:
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           # Likelihood for lambda:
                           sum(dnorm(g_bc(y, lambda = l_bc),
                                     mean = X%*%theta,
                                     sd = sigma_epsilon, log = TRUE)) +
                             # This is the prior on lambda, truncated to [0, 2]
                             dnorm(l_bc, mean = 1/2, sd = 1/2, log = TRUE)
                         },
                         w = 1/2, m = 50, lower = 0, upper = 2)

      # Update z:
      z = g_bc(y, lambda = lambda)
      Xtz = crossprod(X, z)
    }
    # Block 2: sample the scale adjustment (SD)
    # SSR_psi = sum(z^2) - crossprod(z, X%*%solve(crossprod(X) + diag(1/psi, p))%*%crossprod(X,z))
    SSR_psi = sum(z^2) - crossprod(1/sqrt(diagXtX + 1/psi)*Xtz)
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = .001 + n/2,
                                  rate = .001 + SSR_psi/2))

    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    Q_theta = 1/sigma_epsilon^2*diagXtX + 1/(sigma_epsilon^2*psi)
    ell_theta = 1/sigma_epsilon^2*Xtz
    theta = rnorm(n = p,
                  mean = Q_theta^-1*ell_theta,
                  sd = sqrt(Q_theta^-1))
    #----------------------------------------------------------------------------
    # Block 4: sample the smoothing parameter
    if(sample_psi){
      psi = 1/rgamma(n = 1,
                     shape = 0.01 + p/2,
                     rate = 0.01 + sum(theta^2)/(2*sigma_epsilon^2))
    }
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Posterior samples of the model parameters:
        post_theta[isave,] = theta

        # Predictive samples of ytilde:
        ztilde = splinefun(x, X%*%theta)(x_test) +
          sigma_epsilon*rnorm(n = length(x_test))
        post_ypred[isave,] = g_inv_bc(ztilde, lambda = lambda)

        # Posterior samples of the transformation:
        post_g[isave,] = g_bc(y0, lambda = lambda)
        post_lambda[isave] = lambda

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
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
    post_g = post_g, post_lambda = post_lambda,
    model = 'sbsm_bc', y = y, X = X, psi = psi))
}
#' Bayesian Gaussian processes with a Box-Cox transformation
#'
#' MCMC sampling for Bayesian Gaussian process regression with a
#' (known or unknown) Box-Cox transformation.
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
#' @param lambda Box-Cox transformation; if NULL, estimate this parameter
#' @param sample_lambda logical; if TRUE, sample lambda, otherwise
#' use the fixed value of lambda above or the MLE (if lambda unspecified)
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{locs_test}
#' \item \code{fit_gp} the fitted \code{GpGp_fit} object, which includes
#' covariance parameter estimates and other model information
#' \item \code{post_ypred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at \code{locs_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{post_lambda} \code{nsave} posterior samples of lambda
#' \item \code{model}: the model fit (here, \code{bgp_bc})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides Bayesian inference for
#' transformed Gaussian processes. The transformation is
#' parametric from the Box-Cox family, which has one parameter \code{lambda}.
#' That parameter may be fixed in advanced or learned from the data.
#' For computational efficiency, the Gaussian process parameters are
#' fixed at point estimates, and the latent Gaussian process is only sampled
#' when \code{emp_bayes} = FALSE.
#'
#' @note Box-Cox transformations may be useful in some cases, but
#' in general we recommend the nonparametric transformation (with
#' Monte Carlo, not MCMC sampling) in \code{\link{sbgp}}.
#'
#' @examples
#' \donttest{
#' # Simulate some data:
#' n = 200 # sample size
#' x = seq(0, 1, length = n) # observation points
#'
#' # Transform a noisy, periodic function:
#' y = g_inv_bc(
#'   sin(2*pi*x) + sin(4*pi*x) + rnorm(n, sd = .5),
#'              lambda = .5) # Signed square-root transformation
#'
#' # Package we use for fast computing w/ Gaussian processes:
#' library(GpGp)
#'
#' # Fit a Bayesian Gaussian process with Box-Cox transformation:
#' fit = bgp_bc(y = y, locs = x)
#' names(fit) # what is returned
#' coef(fit) # estimated regression coefficients (here, just an intercept)
#' class(fit$fit_gp) # the GpGp object is also returned
#' round(quantile(fit$post_lambda), 3) # summary of unknown Box-Cox parameter
#'
#' # Plot the model predictions (point and interval estimates):
#' pi_y = t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
#' plot(x, y, type='n', ylim = range(pi_y,y),
#'      xlab = 'x', ylab = 'y', main = paste('Fitted values and prediction intervals'))
#' polygon(c(x, rev(x)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
#' lines(x, y, type='p')
#' lines(x, fitted(fit), lwd = 3)
#' }
#'
# #' @import GpGp fields
#' @export
bgp_bc = function(y, locs,
                   X = NULL,
                   covfun_name = "matern_isotropic",
                   locs_test = locs,
                   X_test = NULL,
                   nn = 30,
                   emp_bayes = TRUE,
                   lambda = NULL,
                   sample_lambda = TRUE,
                   nsave = 1000,
                   nburn = 1000,
                   nskip = 0){

  # For testing:
  # X = matrix(1, nrow = length(y)); covfun_name = "matern_isotropic"; locs_test = locs; X_test = X; nn = 30; sample_lambda  = FALSE; lambda = NULL; nsave = 1000; nburn = 1000; nskip = 0;

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
  n = length(y);
  locs = as.matrix(locs); d = ncol(locs)

  # Testing data:
  locs_test = as.matrix(locs_test); n_test = nrow(locs_test)

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

  # Recurring term
  y0 = sort(unique(y))
  #----------------------------------------------------------------------------
  print('Initial GP fit...')
  # Initial GP fit:
  fit_gp = GpGp::fit_model(y = y,
                           locs = locs,
                           X = X,
                           covfun_name = covfun_name,
                           m_seq = nn,
                           silent = TRUE)

  # Fitted values for observed data:
  z_hat = GpGp::predictions(fit = fit_gp,
                            locs_pred = locs,
                            X_pred  = X)

  # First and last define the nugget variance:
  # NOTE: this works for most (but not all!) covariance functions in GpGp!
  sigma_epsilon = sqrt(fit_gp$covparms[1]*
                         fit_gp$covparms[length(fit_gp$covparms)])


  # Optimize the Box-Cox parameter:
  if(is.null(lambda)){

    # conditional MLE:
    opt = optim(par = 0.5,
                fn = function(l_bc){
                  z = g_bc(y, lambda = l_bc)
                  mean((z - z_hat)^2)
                }, method  = "L-BFGS-B", lower = 0.01, upper = 1
    )
    if(opt$convergence  ==  0){
      # Successful optimization
      lambda = opt$par
    } else {
      warning('Optimization failed...setting lambda = 1/2')
      lambda = 1/2
    }
  }

  # Latent data:
  z = g_bc(y, lambda = lambda)
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
  z_hat = GpGp::predictions(fit = fit_gp,
                            locs_pred = locs,
                            X_pred  = X,
                            m = nn)

  # Fitted values for testing data:
  if(isTRUE(all.equal(X, X_test)) &&
     isTRUE(all.equal(locs, locs_test))){
    # If the testing and training data are identical,
    # then there is no need to apply a separate predict function
    z_test = z_hat
  } else {
    z_test = GpGp::predictions(fit = fit_gp,
                               locs_pred = locs_test,
                               X_pred  = X_test,
                               m = nn)
  }

  # Nuggest variance:
  sigma_epsilon = sqrt(fit_gp$covparms[1]*
                         fit_gp$covparms[length(fit_gp$covparms)])

  # Estimated coefficients:
  theta = fit_gp$betahat
  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))
  post_lambda = rep(NA, nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Sample the transformation
    if(sample_lambda){
      # Sample lambda:
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           # Likelihood for lambda:
                           sum(dnorm(g_bc(y, lambda = l_bc),
                                     mean = z_hat,
                                     sd = sigma_epsilon, log = TRUE)) +
                             # This is the prior on lambda, truncated to [0, 2]
                             dnorm(l_bc, mean = 1/2, sd = 1/2, log = TRUE)
                         },
                         w = 1/2, m = 50, lower = 0, upper = 2)

      # Update z:
      z = g_bc(y, lambda = lambda)
    }
    #----------------------------------------------------------------------------
    # Sample the error SD
    # sigma_epsilon = 1/sqrt(rgamma(n = 1,
    #                               shape = .001 + n/2,
    #                               rate = .001 +
    #                                 sum((z - z_hat)^2)/2
    # ))
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Predictive samples of ytilde:
        if(emp_bayes){
          ztilde = z_test + sigma_epsilon*rnorm(n = n_test)
        } else {
          ztilde = GpGp::cond_sim(fit = fit_gp,
                                  locs_pred = locs_test,
                                  X_pred = X_test,
                                  m = nn)
        }
        post_ypred[isave,] = g_inv_bc(ztilde, lambda = lambda)

        # Posterior samples of the transformation:
        post_g[isave,] = g_bc(y0, lambda = lambda)
        post_lambda[isave] = lambda

        # And reset the skip counter:
        skipcount = 0
      }
    }
  }

  return(list(
    coefficients = theta,
    fitted.values = colMeans(post_ypred),
    fit_gp = fit_gp,
    post_ypred = post_ypred,
    post_g = post_g, post_lambda = post_lambda,
    model = 'bgp_bc', y = y, X = X))
}
#' Bayesian quantile regression
#'
#' MCMC sampling for Bayesian quantile regression.
#' An asymmetric Laplace distribution is assumed for the errors,
#' so the regression models targets the specified quantile.
#' A g-prior is assumed for the regression coefficients.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors (no intercept)
#' @param tau the target quantile (between zero and one)
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param psi prior variance (g-prior)
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' \item \code{model}: the model fit (here, \code{bqr})
#' }
#' as well as the arguments passed
#'
#' @note The asymmetric Laplace distribution is advantageous because
#' it links the regression model (\code{X\%*\%theta}) to a pre-specified
#' quantile (\code{tau}). However, it is often a poor model for
#' observed data, and the semiparametric version \code{\link{sbqr}} is
#' recommended in general.
#'
#' @note An intercept is automatically added to \code{X} and
#' \code{X_test}. The coefficients reported do *not* include
#' this intercept parameter.
#'
#' @examples
#' # Simulate some heteroskedastic data (no transformation):
#' dat = simulate_tlm(n = 100, p = 5, g_type = 'box-cox', heterosked = TRUE, lambda = 1)
#' y = dat$y; X = dat$X # training data
#' y_test = dat$y_test; X_test = dat$X_test # testing data
#'
#' # Target this quantile:
#' tau = 0.05
#'
#' # Fit the Bayesian quantile regression model:
#' fit = bqr(y = y, X = X, tau = tau, X_test = X_test)
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
#'
#' # The posterior predictive checks usually do not pass!
#' # try ?sbqr instead...
#'
# #' @importFrom statmod rinvgauss
#' @export
bqr = function(y, X, tau = 0.5,
               X_test = X,
               psi = length(y),
               nsave = 1000,
               nburn = 1000,
               nskip = 0,
               verbose = TRUE){

  # For testing:
  # tau = 0.5; psi = length(y); nsave = 1000; nburn = 1000; nskip = 0; verbose = TRUE

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
  y0 = sort(unique(y))
  XtX = crossprod(X)
  #----------------------------------------------------------------------------
  # Initialize the parameters:

  # Coefficients:
  theta = chol2inv(chol(XtX))%*%crossprod(X,y)
  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_theta = array(NA, c(nsave, p))
  post_ypred = post_qtau = array(NA, c(nsave, n_test))

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: parameter expansion
    xi = 1/statmod::rinvgauss(n = n,
                              mean = sqrt((2 + a_tau^2/b_tau^2)/((y - X%*%theta)^2/b_tau^2)),
                              shape = 2 + a_tau^2/b_tau^2)
    # xi = sapply(1:n, function(i){
    #   rgig(n = 1,
    #        lambda = 0.5,
    #        chi = (y[i] - X[i,]%*%theta)^2/b_tau^2,
    #        psi = 2 + a_tau^2/b_tau^2
    #   )
    # })
    #----------------------------------------------------------------------------
    # Block 2: sample the regression coefficients
    Q_theta = crossprod(X/sqrt(b_tau^2*xi)) + 1/psi*XtX # t(X)%*%diag(1/(b_tau^2*xi))%*%X + 1/psi*XtX
    ell_theta = crossprod(X/(b_tau^2*xi), y - a_tau*xi) # t(X)%*%diag(1/(b_tau^2*xi))%*%(y - a_tau*xi)
    ch_Q = chol(Q_theta)
    theta = backsolve(ch_Q,
                      forwardsolve(t(ch_Q), ell_theta) +
                        rnorm(p + 1))
    #----------------------------------------------------------------------------
    # Block 3: sample the error SD
    # sigma_epsilon = 1/sqrt(rgamma(n = 1,
    #                               shape = .001 + n/2 + p/2,
    #                               rate = .001 +
    #                                 sum((z - X%*%theta)^2)/2 +
    #                                 crossprod(theta, XtX)%*%theta/(2*psi)
    # ))
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Posterior samples of the model parameters:
        post_theta[isave,] = theta[-1]

        # Quantile at the testing point:
        post_qtau[nsi - nburn,] = X_test%*%theta

        # Predictive samples of ytilde:
        xi_test = rexp(n = n_test, rate = 1)
        post_ypred[isave,] = X_test%*%theta + a_tau*xi_test + b_tau*sqrt(xi_test)*rnorm(n = n_test)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
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
    model = 'bqr', y = y, X = X[,-1], X_test = X_test[,-1], psi = psi, tau = tau))
}
#' Bayesian linear model with a Box-Cox transformation and a horseshoe prior
#'
#' MCMC sampling for Bayesian linear regression with 1) a
#' (known or unknown) Box-Cox transformation and 2) a horseshoe prior for
#' the (possibly high-dimensional) regression coefficients.
#'
#' @param y \code{n x 1} vector of observed counts
#' @param X \code{n x p} matrix of predictors (no intercept)
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param lambda Box-Cox transformation; if NULL, estimate this parameter
#' @param sample_lambda logical; if TRUE, sample lambda, otherwise
#' use the fixed value of lambda above or the MLE (if lambda unspecified)
#' @param only_theta logical; if TRUE, only return posterior draws of the
#' regression coefficients (for speed)
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' \item \code{post_lambda}: \code{nsave} posterior samples of lambda
#' \item \code{post_sigma}: \code{nsave} posterior samples of sigma
#' \item \code{model}: the model fit (here, \code{blm_bc_hs})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed linear model via MCMC sampling. The transformation is
#' parametric from the Box-Cox family, which has one parameter \code{lambda}.
#' That parameter may be fixed in advanced or learned from the data.
#'
#' The horseshoe prior is especially useful for high-dimensional settings with
#' many (possibly correlated) covariates. This function
#' uses a fast Cholesky-forward/backward sampler when \code{p < n}
#' and the Bhattacharya et al. (<https://doi.org/10.1093/biomet/asw042>) sampler
#' when \code{p > n}. Thus, the sampler can scale linear in \code{n}
#' (for fixed/small \code{p}) or linear in \code{p} (for fixed/small \code{n}).
#'
#' @note Box-Cox transformations may be useful in some cases, but
#' in general we recommend the nonparametric transformation in \code{\link{sblm_hs}}.
#'
#' @note An intercept is automatically added to \code{X} and
#' \code{X_test}. The coefficients reported do *not* include
#' this intercept parameter, since it is not identified
#' under more general transformation models (e.g., \code{\link{sblm_hs}}).
#'
#' @examples
#' # Simulate data from a transformed (sparse) linear model:
#' dat = simulate_tlm(n = 100, p = 50, g_type = 'step', prop_sig = 0.1)
#' y = dat$y; X = dat$X # training data
#'
#' hist(y, breaks = 25) # marginal distribution
#'
#' # Fit the Bayesian linear model with a Box-Cox transformation & a horseshoe prior:
#' fit = blm_bc_hs(y = y, X = X, verbose = FALSE)
#' names(fit) # what is returned
#'
#' @export
blm_bc_hs = function(y, X, X_test = X,
                     lambda = NULL,
                     sample_lambda = TRUE,
                     only_theta = FALSE,
                     nsave = 1000,
                     nburn = 1000,
                     nskip = 0,
                     verbose = TRUE){

  # For testing:
  # X_test = X; sample_lambda  = FALSE; lambda = NULL; nsave = 1000; nburn = 1000; nskip = 0; verbose = TRUE

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

  # Here, we require that lambda is either input or sampled
  if(is.null(lambda) & !sample_lambda) stop('If sample_lambda = FALSE, the user must input a value of lambda')

  # Hyperparameters for Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001
  #----------------------------------------------------------------------------
  # Recurring terms:
  y0 = sort(unique(y))
  XtX = crossprod(X)
  #----------------------------------------------------------------------------
  # Initialize the parameters:

  # The scale is not well-identified for unknown lambda, so initialize at one
  sigma_epsilon = 1

  # Box-Cox parameter, if unspecified
  if(is.null(lambda)) lambda = 0.5

  # Latent data:
  z = g_bc(y, lambda = lambda)
  Xtz = crossprod(X, z) # only changes if sample_lambda = TRUE

  # Coefficients:
  if(p < n){
    theta = chol2inv(chol(XtX))%*%Xtz
  } else {
    theta = sampleFastGaussian(Phi = X/sigma_epsilon,
                               Ddiag = rep(100, p+1),
                               alpha = z/sigma_epsilon)
  }
  # Horseshoe prior parameters:
  #         theta_j ~ N(0, lambda_j^2)
  #         lambda_j ~ C+(0, tau)
  #         tau ~ C+(0,1)
  # Local precision parameters: eta_lambda_j = 1/lambda_j^2
  eta_lambda_j = c(10^-6, # flat prior for intercept
                   1/square_stabilize(theta[-1]))
  xi_lambda_j = rep(1, p)

  # Global precision parameters: eta_tau = 1/tau^2
  eta_tau = xi_tau = 1
  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_theta = array(NA, c(nsave, p))
  if(only_theta){
    post_ypred = post_g = post_lambda = post_sigma = NULL
  } else {
    post_ypred = array(NA, c(nsave, n_test))
    post_g = array(NA, c(nsave, length(y0)))
    post_lambda = post_sigma = rep(NA, nsave)
  }

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(sample_lambda){
      # Sample lambda:
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           # Likelihood for lambda:
                           sum(dnorm(g_bc(y, lambda = l_bc),
                                     mean = X%*%theta,
                                     sd = sigma_epsilon, log = TRUE)) +
                             # This is the prior on lambda, truncated to [0, 2]
                             dnorm(l_bc, mean = 1/2, sd = 1/2, log = TRUE)
                         },
                         w = 1/2, m = 50, lower = 0, upper = 2)

      # Update z:
      z = g_bc(y, lambda = lambda)
      Xtz = crossprod(X, z)
    }
    #----------------------------------------------------------------------------
    # Recurring terms for Blocks 2-3
    ch_Q = chol(XtX + diag(eta_lambda_j, p + 1))
    fsolve_theta = forwardsolve(t(ch_Q), Xtz)

    # Block 2: sample the error SD (unconditional on theta!)
    SSR_hs = sum(z^2) - crossprod(Xtz, backsolve(ch_Q, fsolve_theta))
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = a_sigma + n/2,
                                  rate = b_sigma + SSR_hs/2))
    #----------------------------------------------------------------------------
    # Block 3: sample the regression coefficients
    if(p >= n){
      # Fast sampler for p >= n (BHATTACHARYA et al., 2016)
      theta = sampleFastGaussian(Phi = X/sigma_epsilon,
                                 Ddiag = sigma_epsilon^2/eta_lambda_j,
                                 alpha = z/sigma_epsilon)
    } else {
      # Fast sampler for p < n (Rue, 2001)
      theta = backsolve(ch_Q/sigma_epsilon,
                        fsolve_theta/sigma_epsilon + rnorm(p+1))

      # Previously:
      # ch_Q = 1/sigma_epsilon*chol(XtX + diag(eta_lambda_j, p+1))
      # ell_theta = 1/sigma_epsilon^2*crossprod(X, z)
      # theta = backsolve(ch_Q, forwardsolve(t(ch_Q), ell_theta) + rnorm(p+1))
    }
    #----------------------------------------------------------------------------
    # Block 4: sample the horseshoe prior variance parameters

    # Precision parameters & parameter expansions (local):
    eta_lambda_j[-1] = rgamma(n = p,
                              shape = 1,
                              rate = xi_lambda_j + square_stabilize(theta[-1]/sigma_epsilon)/2
    )
    xi_lambda_j = rgamma(n = p,
                         shape = 1,
                         rate = eta_lambda_j[-1] + eta_tau)

    # Precision parameters & parameter expansions (global):
    eta_tau = rgamma(n = 1,
                     shape = p/2 + 1/2,
                     rate = sum(xi_lambda_j) + xi_tau)
    xi_tau = rgamma(n = 1,
                    shape = 1,
                    rate = eta_tau + 1)
    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Posterior samples of the model parameters:
        post_theta[isave,] = theta[-1]

        if(!only_theta){
          # Predictive samples of ytilde:
          ztilde = X_test%*%theta + sigma_epsilon*rnorm(n = n_test)
          post_ypred[isave,] = g_inv_bc(ztilde, lambda = lambda)

          # Posterior samples of the transformation:
          post_g[isave,] = g_bc(y0, lambda = lambda)
          post_lambda[isave] = lambda

          # Posterior samples of the error SD:
          post_sigma[isave] = sigma_epsilon
        }

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
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
    fitted.values = ifelse(only_theta, NA, colMeans(post_ypred)),
    post_theta = post_theta,
    post_ypred = post_ypred,
    post_g = post_g, post_lambda = post_lambda, post_sigma = post_sigma,
    model = 'blm_bc_hs', y = y, X = X[,-1], X_test = X_test[,-1], sample_lambda = sample_lambda, only_theta = only_theta))
}
