#' Semiparametric Bayesian linear model with horseshoe priors for
#' high-dimensional data
#'
#' MCMC sampling for semiparametric Bayesian linear regression with
#' 1) an unknown (nonparametric) transformation and 2) a horseshoe prior for
#' the (possibly high-dimensional) regression coefficients. Here, unlike \code{\link{sblm}},
#' Gibbs sampling is needed for the regression coefficients and the horseshoe
#' prior variance components. The transformation \code{g} is still sampled
#' unconditionally on the regression coefficients, which provides a more
#' efficient blocking within the Gibbs sampler.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors (no intercept)
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param approx_g logical; if TRUE, apply large-sample
#' approximation for the transformation
#' @param init_screen for the initial approximation, number of covariates
#' to pre-screen (necessary when \code{p > n}); if NULL, use \code{n/log(n)}
#' @param pilot_hs logical; if TRUE, use a short pilot run with a horseshoe
#' prior to estimate the marginal CDF of the latent z (otherwise, use a sparse Laplace approximation)
#' @param nsave number of MCMC simulations to save
#' @param nburn number of MCMC iterations to discard
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
#' \item \code{model}: the model fit (here, \code{sblm_hs})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed linear model with horseshoe priors using efficiently-blocked Gibbs sampling.
#' The transformation is modeled as unknown and learned jointly
#' with the regression coefficients (unless \code{approx_g} = TRUE, which then uses
#' a point approximation). This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#'
#' The horseshoe prior is especially useful for high-dimensional settings with
#' many (possibly correlated) covariates. Compared to sparse or spike-and-slab
#' alternatives (see \code{\link{sblm_ssvs}}), the horseshoe prior
#' delivers more scalable computing in \code{p}. This function
#' uses a fast Cholesky-forward/backward sampler when \code{p < n}
#' and the Bhattacharya et al. (<https://doi.org/10.1093/biomet/asw042>) sampler
#' when \code{p > n}. Thus, the sampler can scale linear in \code{n}
#' (for fixed/small \code{p}) or linear in \code{p} (for fixed/small \code{n}).
#' Empirically, the horseshoe prior performs best under sparse regimes,
#' i.e., when the number of true signals (nonzero regression coefficients)
#' is a small fraction of the total number of variables.
#'
#' To learn the transformation, \code{SeBR} infers the marginal CDF
#' of the latent data model \code{Fz} by integrating over the covariates
#' \code{X} and the coefficients \code{theta}. When \code{fixedX = TRUE}, the
#' \code{X} averaging is empirical; otherwise it uses the Bayesian bootstrap (\code{\link{bb}}).
#' By default, \code{fixedX} is set to \code{FALSE} for smaller datasets (\code{n < 500})
#' and \code{TRUE} for larger datasets. When \code{pilot_hs = TRUE},
#' the algorithm fits an initial linear regression model
#' with a horseshoe prior (\code{\link{blm_bc_hs}}) to transformed data
#' (under a preliminary point estimate of the transformation) and
#' uses that posterior distribution to integrate over \code{theta}.
#' Otherwise, this marginalization is done using a sparse Laplace approximation
#' for speed and simplicity.
#'
#' @note The location (intercept) and scale (\code{sigma_epsilon}) are
#' not identified, so any intercepts in \code{X} and \code{X_test} will
#' be removed. The model-fitting *does* include an internal location-scale
#' adjustment, but the function only outputs inferential summaries for the
#' identifiable parameters.
#'
#' @examples
#' \donttest{
#' # Simulate data from a transformed (sparse) linear model:
#' dat = simulate_tlm(n = 100, p = 50, g_type = 'step', prop_sig = 0.1)
#' y = dat$y; X = dat$X # training data
#' y_test = dat$y_test; X_test = dat$X_test # testing data
#'
#' hist(y, breaks = 25) # marginal distribution
#'
#' # Fit the semiparametric Bayesian linear model with a horseshoe prior:
#' fit = sblm_hs(y = y, X = X, X_test = X_test)
#' names(fit) # what is returned
#'
#' # Evaluate posterior predictive means and intervals on the testing data:
#' plot_pptest(fit$post_ypred, y_test,
#'             alpha_level = 0.10) # coverage should be about 90%
#'
#' # Check: correlation with true coefficients
#' cor(dat$beta_true, coef(fit))
#'
#' # Compute 95% credible intervals for the coefficients:
#' ci_theta = t(apply(fit$post_theta, 2, quantile, c(0.05/2, 1 - 0.05/2)))
#'
#' # True positive/negative rates for "selected" coefficients:
#' selected = ((ci_theta[,1] >0 | ci_theta[,2] < 0)) # intervals exclude zero
#' sigs_true = dat$beta_true != 0 # true signals
#' (TPR = sum(selected & sigs_true)/sum(sigs_true))
#' (TNR = sum(!selected & !sigs_true)/sum(!sigs_true))
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
#' @importFrom stats cor median
#' @export
sblm_hs = function(y, X, X_test = X,
                   fixedX = (length(y) >= 500),
                   approx_g = FALSE,
                   init_screen = NULL,
                   pilot_hs = FALSE,
                   nsave = 1000,
                   nburn = 1000,
                   ngrid = 100,
                   verbose = TRUE){

  # For testing:
  # fixedX = FALSE; init_screen = NULL; pilot_hs = FALSE; approx_g = FALSE; nsave = 1000; nburn = 1000; verbose = TRUE; ngrid = 100

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

  # Hyperparameters for Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001
  #----------------------------------------------------------------------------
  # Pre-screen variable for the initial approximation
  # based on marginal rank (spearman) correlations, which are invariant to g()

  # Pre-screening: how many variables to include in the initial approx?
  if(is.null(init_screen))
    init_screen = min(p, floor(n/log(n)))

  # Absolute rank correlations:
  acor = abs(cor(y, X, method = 'spearman'))

  # Find the cutoff and take the top init_screen:
  cutoff = sort(acor, decreasing = TRUE)[init_screen]
  var_screen = which(acor >= cutoff) # including the intercept
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)

  # Evaluate at the unique y-values:
  y0 = sort(unique(y))
  Fy_eval = Fy(y0)

  # Initial transformation: fix Fz() = qnorm()
  z = qnorm(Fy(y))

  # Recurring terms (for pre-screened variables):
  Sigma_hat_unscaled = chol2inv(chol(crossprod(X[,var_screen]))) # unscaled covariance (pre-screened)
  xt_Sigma_hat_unscaled_x = sapply(1:n, function(i)
    crossprod(X[i,var_screen], Sigma_hat_unscaled)%*%X[i,var_screen]) # sandwiched by X[i,] (pre-screened)

  # OLS initialize (pre-screened) coefficients
  theta_hat = Sigma_hat_unscaled%*%crossprod(X[,var_screen], z) # point estimate (pre-screened)

  # Moments of Z|X:
  mu_z = X[,var_screen]%*%theta_hat
  sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

  # Grid of values for the CDF of z:
  z_grid = sort(unique(
    sapply(range(mu_z), function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = ngrid),
            mean = xtemp,
            sd = sqrt(1 + median(xt_Sigma_hat_unscaled_x)))
    })
  ))

  # Second pass: update g(), then update coefficients

  # CDF of z:
  Fz_eval = Fz_fun(z = z_grid,
                   weights = rep(1/n, n),
                   mean_vec = mu_z,
                   sd_vec = sigma_z)

  # Check: update the grid if needed
  zcon = contract_grid(z = z_grid,
                       Fz = Fz_eval,
                       lower = 0.001, upper =  0.999)
  z_grid = zcon$z; Fz_eval = zcon$Fz; ngrid = length(z_grid)

  # Transformation:
  g = g_fun(y = y0,
            Fy_eval = Fy_eval,
            z = z_grid,
            Fz_eval = Fz_eval)

  # Update latent data
  z = g(y)

  if(pilot_hs){

    # Pilot run: estimate the posterior of theta under a horseshoe prior (fixed transformation)
    # and use this to compute Fz
    post_theta = blm_bc_hs(y = z, X = X,
                           lambda = 1, sample_lambda = FALSE, # fixed identity transformation
                           only_theta = TRUE,  nsave = 500, nburn = 500, verbose = FALSE)$post_theta # quick run; only returns post_theta

    # Initialize coefficients:
    theta = colMeans(post_theta)

    # Now store F_{Z | X = x_i}(t) for all x_i & all t in z_grid
    #   this term integrates out theta from its posterior:
    Fzx_eval = matrix(0, nrow = ngrid, ncol = n)
    post_Xtheta = tcrossprod(X, post_theta); rm(post_theta)
    zrep = rep(z_grid, times = ncol(post_Xtheta))
    Fzx_eval = sapply(1:n, function(i){
      rowMeans(matrix(pnorm(zrep,
                            mean = rep(post_Xtheta[i,], each = ngrid),
                            sd = 1), nrow = ngrid))
    })

    # Quick clean up
    rm(post_Xtheta,theta_hat,Sigma_hat_unscaled,xt_Sigma_hat_unscaled_x,mu_z,sigma_z)

    # Update the transformation based on this approximation:
    if(approx_g || fixedX){

      # Point estimate:
      Fz_eval = rowMeans(Fzx_eval)

      # Transformation:
      g = g_fun(y = y0,
                Fy_eval = Fy_eval,
                z = z_grid,
                Fz_eval = Fz_eval)

      # Update latent data
      z = g(y)
    }

  } else {

    # updated coefficients (pre-screened)
    theta_hat = Sigma_hat_unscaled%*%crossprod(X[,var_screen], z)

    # Moments of Z|X:
    mu_z = X[,var_screen]%*%theta_hat
    # sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x) # no need to update

    # Initialize coefficients,
    theta = rep(0, p) # sparse otherwise
    theta[var_screen] = theta_hat; rm(theta_hat)
  }
  #----------------------------------------------------------------------------
  # Update for the intercept:
  theta = c(0, theta) # add intercept
  X = cbind(1, X) # update to include intercept
  XtX = crossprod(X); # recurring terms
  Xtz = crossprod(X, z) # changes w/ z (unless approx_g = TRUE)
  #----------------------------------------------------------------------------
  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)
  #----------------------------------------------------------------------------
  # Now initialize the horseshoe prior parameters:
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
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))
  post_sigma = rep(NA, nsave) # for testing

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

        # CDF of Z evaluated on the grid:
        if(pilot_hs){
          Fz_eval = Fzx_eval%*%weights_x
        } else {
          Fz_eval = Fz_fun(z = z_grid,
                           weights = weights_x,
                           mean_vec = mu_z,
                           sd_vec = sigma_z)
        }
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
    # Recurring terms for Blocks 2-3
    ch_Q = chol(XtX + diag(eta_lambda_j, p + 1))
    fsolve_theta = forwardsolve(t(ch_Q), Xtz)
    #----------------------------------------------------------------------------
    # Block 2: sample the scale adjustment (SD)
    SSR_hs = sum(z^2) - crossprod(Xtz, backsolve(ch_Q, fsolve_theta)) # SSR_hs = sum(z^2) - crossprod(z, X%*%chol2inv(chol(XtX + diag(eta_lambda_j, p + 1)))%*%crossprod(X, z))
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

      # Posterior samples of the model parameters:
      post_theta[nsi - nburn,] = theta[-1]/sigma_epsilon # omit the intercept & undo scaling (not identified)

      # Predictive samples of ytilde:
      ztilde = X_test%*%theta + sigma_epsilon*rnorm(n = n_test)
      post_ypred[nsi - nburn,] = g_inv(ztilde)

      # Posterior samples of the transformation:
      post_g[nsi - nburn,] = (g(y0) - theta[1])/sigma_epsilon # undo location/scale (not identified)

      # Posterior samples of sigma:
      post_sigma[nsi - nburn] = sigma_epsilon # for testing
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nsave + nburn)
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
    post_sigma = post_sigma, # for testing
    model = 'sblm_hs', y = y, X = X[,-1], X_test = X_test[,-1], fixedX = fixedX, approx_g = approx_g, init_screen = init_screen, pilot_hs = pilot_hs))
}
#' Semiparametric Bayesian linear model with stochastic search variable selection
#'
#' MCMC sampling for semiparametric Bayesian linear regression with
#' 1) an unknown (nonparametric) transformation and 2) a sparsity prior on
#' the (possibly high-dimensional) regression coefficients. Here, unlike \code{\link{sblm}},
#' Gibbs sampling is used for the variable inclusion indicator variables
#' \code{gamma}, referred to as stochastic search variable selection (SSVS).
#' All remaining terms--including the transformation \code{g}, the regression
#' coefficients \code{theta}, and any predictive draws--are drawn directly from
#' the joint posterior (predictive) distribution.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors (no intercept)
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param psi prior variance (g-prior)
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param approx_g logical; if TRUE, apply large-sample
#' approximation for the transformation
#' @param init_screen for the initial approximation, number of covariates
#' to pre-screen (necessary when \code{p > n}); if NULL, use \code{n/log(n)}
#' @param a_pi shape1 parameter of the (Beta) prior inclusion probability
#' @param b_pi shape2 parameter of the (Beta) prior inclusion probability
#' @param nsave number of MCMC simulations to save
#' @param nburn number of MCMC iterations to discard
#' @param ngrid number of grid points for inverse approximations
#' @param verbose logical; if TRUE, print time remaining
#' @return a list with the following elements:
#' \itemize{
#' \item \code{coefficients} the posterior mean of the regression coefficients
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{X_test}
#' \item \code{selected}: the variables (columns of \code{X}) selected by the median probability model
#' \item \code{pip}: (marginal) posterior inclusion probabilities for each variable
#' \item \code{post_theta}: \code{nsave x p} samples from the posterior distribution
#' of the regression coefficients
#' \item \code{post_gamma}: \code{nsave x p} samples from the posterior distribution
#' of the variable inclusion indicators
#' \item \code{post_ypred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{model}: the model fit (here, \code{sblm_ssvs})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed linear model with sparse g-priors on the regression coefficients.
#' The transformation is modeled as unknown and learned jointly
#' with the regression coefficients (unless \code{approx_g} = TRUE, which then uses
#' a point approximation). This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' By default, \code{fixedX} is set to \code{FALSE} for smaller datasets (\code{n < 500})
#' and \code{TRUE} for larger datasets.
#'
#' The sparsity prior is especially useful for variable selection. Compared
#' to the horseshoe prior version (\code{\link{sblm_hs}}), the sparse g-prior
#' is advantageous because 1) it truly allows for sparse (i.e., exactly zero)
#' coefficients in the prior and posterior, 2) it incorporates covariate
#' dependencies via the g-prior structure, and 3) it tends to perform well
#' under both sparse and non-sparse regimes, while the horseshoe version only
#' performs well under sparse regimes. The disadvantage is that
#' SSVS does not scale nearly as well in \code{p}.
#'
#' Following Scott and Berger (<https://doi.org/10.1214/10-AOS792>),
#' we include a \code{Beta(a_pi, b_pi)} prior on the prior inclusion probability. This term
#' is then sampled with the variable inclusion indicators \code{gamma} in a
#' Gibbs sampling block. All other terms are sampled using direct Monte Carlo
#' (not MCMC) sampling.
#'
#' Alternatively, model probabilities can be computed directly
#' (by Monte Carlo, not MCMC/Gibbs sampling) using \code{\link{sblm_modelsel}}.
#'
#' @note The location (intercept) and scale (\code{sigma_epsilon}) are
#' not identified, so any intercepts in \code{X} and \code{X_test} will
#' be removed. The model-fitting *does* include an internal location-scale
#' adjustment, but the function only outputs inferential summaries for the
#' identifiable parameters.
#'
#' @examples
#' \donttest{
#' # Simulate data from a transformed (sparse) linear model:
#' dat = simulate_tlm(n = 100, p = 15, g_type = 'step')
#' y = dat$y; X = dat$X # training data
#' y_test = dat$y_test; X_test = dat$X_test # testing data
#'
#' hist(y, breaks = 25) # marginal distribution
#'
#' # Fit the semiparametric Bayesian linear model with sparsity priors:
#' fit = sblm_ssvs(y = y, X = X, X_test = X_test)
#' names(fit) # what is returned
#'
#' # Evaluate posterior predictive means and intervals on the testing data:
#' plot_pptest(fit$post_ypred, y_test,
#'             alpha_level = 0.10) # coverage should be about 90%
#'
#' # Check: correlation with true coefficients
#' cor(dat$beta_true, coef(fit))
#'
#' # Selected coefficients under median probability model:
#' fit$selected
#'
#' # True signals:
#' which(dat$beta_true != 0)
#'
#' # Summarize the transformation:
#' y0 = sort(unique(y)) # posterior draws of g are evaluated at the unique y observations
#' plot(y0, fit$post_g[1,], type='n', ylim = range(fit$post_g),
#'      xlab = 'y', ylab = 'g(y)', main = "Posterior draws of the transformation")
#' temp = sapply(1:nrow(fit$post_g), function(s)
#'   lines(y0, fit$post_g[s,], col='gray')) # posterior draws
#' lines(y0, colMeans(fit$post_g), lwd = 3) # posterior mean
#' lines(y, dat$g_true, type='p', pch=2) # true transformation
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
#' @importFrom stats cor median
#' @export
sblm_ssvs = function(y, X, X_test = X,
                     psi = length(y),
                     fixedX = (length(y) >= 500),
                     approx_g = FALSE,
                     init_screen = NULL,
                     a_pi = 1, b_pi = 1,
                     nsave = 1000,
                     nburn = 1000,
                     ngrid = 100,
                     verbose = TRUE){

  # For testing:
  # psi = length(y); init_screen = NULL; a_pi = b_pi = 1; fixedX =FALSE; approx_g = FALSE; nsave = 1000; nburn = 1000; verbose = TRUE; ngrid = 100

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
  X_test = cbind(1,
                 X_test[,!is_intercept_test])

  # Hyperparameters for Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001
  #----------------------------------------------------------------------------
  # Pre-screen variable for the initial approximation
  # based on marginal rank (spearman) correlations, which are invariant to g()

  # Pre-screening: how many variables to include in the initial approx?
  if(is.null(init_screen))
    init_screen = min(p, floor(n/log(n)))

  # Absolute rank correlations:
  acor = abs(cor(y, X, method = 'spearman'))

  # Find the cutoff and take the top init_screen:
  cutoff = sort(acor, decreasing = TRUE)[init_screen]
  gamma = 1.0*(acor >= cutoff) # selection indicator
  var_screen = which(gamma==1) # including the intercept

  # Initialize the prior inclusion probability (restrict to [.1, .9])
  pi_gamma = max(0.1, min(0.9, mean(gamma)))
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)

  # Evaluate at the unique y-values:
  y0 = sort(unique(y))
  Fy_eval = Fy(y0)

  # Initial transformation: fix Fz() = qnorm()
  z = qnorm(Fy(y))

  # Recurring terms (for pre-screened variables):
  Sigma_hat_unscaled = psi/(1+psi)*chol2inv(chol(crossprod(X[,var_screen]))) # unscaled covariance (pre-screened)
  xt_Sigma_hat_unscaled_x = sapply(1:n, function(i)
    crossprod(X[i,var_screen], Sigma_hat_unscaled)%*%X[i,var_screen]) # sandwiched by X[i,] (pre-screened)

  # OLS initialize (pre-screened) coefficients
  theta_hat = Sigma_hat_unscaled%*%crossprod(X[,var_screen], z) # point estimate (pre-screened)

  # Moments of Z|X:
  mu_z = X[,var_screen]%*%theta_hat
  sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

  # Grid of values for the CDF of z:
  z_grid = sort(unique(
    sapply(range(mu_z), function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = ngrid),
            mean = xtemp,
            sd = sqrt(1 + median(xt_Sigma_hat_unscaled_x)))
    })
  ))

  # Second pass: update g(), then update coefficients

  # CDF of z:
  Fz_eval = Fz_fun(z = z_grid,
                   weights = rep(1/n, n),
                   mean_vec = mu_z,
                   sd_vec = sigma_z)

  # Check: update the grid if needed
  zcon = contract_grid(z = z_grid,
                       Fz = Fz_eval,
                       lower = 0.001, upper =  0.999)
  z_grid = zcon$z; Fz_eval = zcon$Fz; ngrid = length(z_grid)

  # Transformation:
  g = g_fun(y = y0,
            Fy_eval = Fy_eval,
            z = z_grid,
            Fz_eval = Fz_eval)

  # Update latent data
  z = g(y)

  # updated coefficients (pre-screened)
  theta_hat = Sigma_hat_unscaled%*%crossprod(X[,var_screen], z)

  # Moments of Z|X:
  mu_z = X[,var_screen]%*%theta_hat
  # sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x) # no need to update

  # Initialize coefficients, including intercept:
  theta = rep(0, p+1) # sparse otherwise
  theta[1+var_screen] = theta_hat # '+ 1' shifts for intercept

  # Add an intercept to X and gamma (inclusion):
  X = cbind(1, X)
  gamma = c(1, gamma) # always include the intercept
  post_probs = gamma # initialize the posterior probabilities
  #----------------------------------------------------------------------------
  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)
  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_theta = array(NA, c(nsave, p))
  post_gamma = array(NA, c(nsave, p))
  pip = rep(0, p)
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))

  # Run the MC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:(nsave + nburn)){

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
    # Block 2: sample the variable inclusion indicators
    #   (in a random order)
    log_o_prior = log(pi_gamma) - log(1 - pi_gamma) # prior log-odds (recurring term)
    for(j in sample(2:(p+1), p)){ # all but the intercept

      # Sample from [gamma[j] | data, gamma[-j]]:
      gamma0 = gamma1 = gamma; # placeholders
      gamma0[j] = 0 # "inactive j"
      gamma1[j] = 1 # "active j"

      # Log-odds of selection, simplified version:
      log_oj = log_o_prior - 1/2*log(1+psi) -(n/2 + a_sigma)*(
        log(SSR_gprior(z, X[, gamma1==1, drop=FALSE], psi) + 2*b_sigma) -
          log(SSR_gprior(z, X[, gamma0==1, drop=FALSE], psi) + 2*b_sigma)
      )

      # # Log-odds of selection, detailed version (slower):
      #
      # # Number of active terms, excluding j:
      # p_gamma_nj = sum(gamma[-j])
      #
      # # Log-likelihood under "inactive j"
      # log_like0 = -p_gamma_nj/2*log(1+psi) -
      #   (n/2 + a_sigma)*log(SSR_gprior(z, X[, gamma0==1, drop=FALSE], psi) + 2*b_sigma)
      #
      # # Log-likelihood under "active j"
      # log_like1 = -(1+p_gamma_nj)/2*log(1+psi) -
      #   (n/2 + a_sigma)*log(SSR_gprior(z, X[, gamma1==1, drop=FALSE], psi) + 2*b_sigma)
      #
      # # Log-prior under "inactive j"
      # log_prior0 = log(1 - pi_gamma)
      #
      # # Log-prior under "active j"
      # log_prior1 = log(pi_gamma)
      #
      # # Log-odds of selection
      # log_oj = (log_like1 + log_prior1) - (log_like0 + log_prior0)

      # Odds:
      o_j = exp(log_oj)

      # Probabilities:
      post_probs[j] = o_j/(1 + o_j)

      # Sample:
      if(runif(1) < post_probs[j]){ #if(runif(1) < o_j/(1 + o_j)){
        gamma[j] = 1
      } else gamma[j] = 0
    }

    # Summarize the selections:
    var_screen = which(gamma == 1) # indices to be included
    p_gamma = length(var_screen) - 1 # number of included variables (not counting intercept)

    # Also sample the inclusion probabilities:
    pi_gamma = rbeta(n = 1,
                     shape1 = a_pi + p_gamma, # sum(gamma[-1]==1)
                     shape2 = b_pi + (p - p_gamma)) # sum(gamma[-1]==0)

    # Subset the covariate (X) terms for selected variables:
    X_gamma = X[,var_screen, drop=FALSE]
    #----------------------------------------------------------------------------
    # Block 3: sample the scale adjustment (SD)

    # Gamma rate requires the sum-scaled-residuals (SSR) under g-prior:
    SSR_psi = SSR_gprior(z, X_gamma, psi) # sum(z^2) - psi/(psi+1)*crossprod(z, X_gamma%*%XtXinv_gamma%*%crossprod(X_gamma, z))

    # Note: this assumes Gamma(a_sigma, b_sigma) prior on the error precision
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = a_sigma + n/2,
                                  rate = b_sigma + SSR_psi/2))
    #----------------------------------------------------------------------------
    # Block 4: sample the regression coefficients
    theta = rep(0, p + 1) # reset every time
    ch_Q = chol(1/sigma_epsilon^2*(1+psi)/(psi)*crossprod(X_gamma))
    ell_theta = 1/sigma_epsilon^2*crossprod(X_gamma, z)
    theta[var_screen] = backsolve(ch_Q,
                                  forwardsolve(t(ch_Q), ell_theta) +
                                    rnorm(p_gamma + 1)) # intercept too

    # Store the MCMC:
    if(nsi > nburn){

      # Posterior samples of the model parameters:
      post_theta[nsi - nburn,] = theta[-1]/sigma_epsilon # omit the intercept & undo scaling (not identified)

      # Posterior samples of the variable inclusion indicators:
      post_gamma[nsi - nburn,] = gamma[-1] # exclude intercept

      # Posterior inclusion probabilities:
      pip = pip + 1/nsave*post_probs[-1]  # exclude intercept

      # Predictive samples of ytilde:
      ztilde = X_test%*%theta + sigma_epsilon*rnorm(n = n_test)
      post_ypred[nsi - nburn,] = g_inv(ztilde)

      # Posterior samples of the transformation:
      post_g[nsi - nburn,] = (g(y0) - theta[1])/sigma_epsilon # undo location/scale (not identified)
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nsave + nburn)
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
    selected = which(pip > 0.5),
    pip = pip,
    post_theta = post_theta,
    post_gamma = post_gamma,
    post_ypred = post_ypred,
    post_g = post_g,
    model = 'sblm_ssvs', y = y, X = X[,-1], X_test = X_test[,-1], psi = psi, fixedX = fixedX, approx_g = approx_g, init_screen = init_screen, a_pi = a_pi, b_pi = b_pi))
}
#' Model selection for semiparametric Bayesian linear regression
#'
#' Compute model probabilities for semiparametric Bayesian linear regression with
#' 1) an unknown (nonparametric) transformation and 2) a sparsity prior on
#' the regression coefficients. The model probabilities are computed
#' using direct Monte Carlo (not MCMC) sampling.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors (no intercept)
#' @param prob_inclusion prior inclusion probability for each variable
#' @param psi prior variance (g-prior)
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param init_screen for the initial approximation, number of covariates
#' to pre-screen (necessary when \code{p > n}); if NULL, use \code{n/log(n)}
#' @param nsave number of Monte Carlo simulations
#' @param override logical; if TRUE, the user may override the default
#' cancellation of the function call when \code{p > 15}
#' @param ngrid number of grid points for inverse approximations
#' @param verbose logical; if TRUE, print time remaining
#' @return a list with the following elements:
#' \itemize{
#' \item \code{post_probs} the posterior probabilities for each model
#' \item \code{all_models}: \code{2^p x p} matrix where each row corresponds to a model
#' from \code{post_probs} and each column indicates inclusion (TRUE) or exclusion (FALSE) for that variable
#' \item \code{model}: the model fit (here, \code{sblm_modelsel})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian model selection for a
#' transformed linear model with sparse g-priors on the regression coefficients.
#' The transformation is modeled as unknown and learned jointly
#' with the model probabilities. This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' By default, \code{fixedX} is set to \code{FALSE} for smaller datasets (\code{n < 500})
#' and \code{TRUE} for larger datasets.
#'
#' Enumeration of all possible subsets is computationally demanding and
#' should be reserved only for small \code{p}. The function will exit for
#' \code{p > 15} unless \code{override = TRUE}.
#'
#' This function exclusively computes model probabilities and does not
#' provide other coefficient inference or prediction. These additions
#' would be straightforward, but are omitted to save on computing time.
#' For prediction, inference, and computation
#' with moderate to large \code{p}, use \code{\link{sblm_ssvs}}.
#'
#' @note The location (intercept) and scale (\code{sigma_epsilon}) are
#' not identified, so any intercept in \code{X} will be removed.
#' The model-fitting *does* include an internal location-scale
#' adjustment, but the model probabilities only refer to the
#' non-intercept variables in \code{X}.
#'
#' @examples
#' \donttest{
#' # Simulate data from a transformed (sparse) linear model:
#' dat = simulate_tlm(n = 100, p = 5, g_type = 'beta')
#' y = dat$y; X = dat$X
#'
#' hist(y, breaks = 25) # marginal distribution
#'
#' # Package for conveniently computing all subsets:
#' library(plyr)
#'
#' # Fit the semiparametric Bayesian linear model with model selection:
#' fit = sblm_modelsel(y = y, X = X)
#' names(fit) # what is returned
#'
#' # Summarize the probabilities of each model (by size):
#' plot(rowSums(fit$all_models), fit$post_probs,
#'      xlab = 'Model sizes', ylab = 'p(model | data)',
#'      main = 'Posterior model probabilities', pch = 2, ylim = c(0,1))
#'
#' # Highest probability model:
#' hpm = which.max(fit$post_probs)
#' fit$post_probs[hpm] # probability
#' which(fit$all_models[hpm,]) # which variables
#' which(dat$beta_true != 0) # ground truth
#' }
#'
#' @importFrom stats dbinom
#' @export
sblm_modelsel = function(y, X,
                         prob_inclusion = 0.5,
                         psi = length(y),
                         fixedX = (length(y) >= 500),
                         init_screen = NULL,
                         nsave = 1000,
                         override = FALSE,
                         ngrid = 100,
                         verbose = TRUE){

  # For testing:
  # prob_inclusion = 0.5; psi = length(y); fixedX = FALSE; init_screen = NULL;  nsave = 1000; verbose = TRUE; ngrid = 100

  # Initial checks:
  if(!is.matrix(X)) stop("X must be a matrix (rows = observations, columns = variables)")
  if(nrow(X) != length(y)) stop('the length of y must equal nrow(X)')

  # Data dimensions:
  n = length(y) # number of observations
  p = ncol(X) # number of variables

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

  # Naming:
  if(is.null(colnames(X))) colnames(X) = paste('X', 1:p, sep='')

  # Add a warning!
  if(p > 15 && !override)
    stop("This function is not recommended when p > 15!
          Try 'sblm_ssvs' or override this break by re-running 'sblm_modelsel' with 'override = TRUE'")
  #----------------------------------------------------------------------------
  # Pre-screen variable for the initial approximation
  # based on marginal rank (spearman) correlations, which are invariant to g()

  # Pre-screening: how many variables to include in the initial approx?
  if(is.null(init_screen))
    init_screen = min(p, floor(n/log(n)))

  # Absolute rank correlations:
  acor = abs(cor(y, X, method = 'spearman'))

  # Find the cutoff and take the top init_screen:
  cutoff = sort(acor, decreasing = TRUE)[init_screen]
  gamma = 1.0*(acor >= cutoff) # selection indicator
  var_screen = which(gamma==1) # including the intercept

  # Hyperparameters for Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)

  # Evaluate at the unique y-values:
  y0 = sort(unique(y))
  Fy_eval = Fy(y0)

  # Initial transformation: fix Fz() = qnorm()
  z = qnorm(Fy(y))

  # Recurring terms (for pre-screened variables):
  Sigma_hat_unscaled = psi/(1+psi)*chol2inv(chol(crossprod(X[,var_screen]))) # unscaled covariance (pre-screened)
  xt_Sigma_hat_unscaled_x = sapply(1:n, function(i)
    crossprod(X[i,var_screen], Sigma_hat_unscaled)%*%X[i,var_screen]) # sandwiched by X[i,] (pre-screened)

  # OLS initialize (pre-screened) coefficients
  theta_hat = Sigma_hat_unscaled%*%crossprod(X[,var_screen], z) # point estimate (pre-screened)

  # Moments of Z|X:
  mu_z = X[,var_screen]%*%theta_hat
  sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x)

  # Grid of values for the CDF of z:
  z_grid = sort(unique(
    sapply(range(mu_z), function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = ngrid),
            mean = xtemp,
            sd = sqrt(1 + median(xt_Sigma_hat_unscaled_x)))
    })
  ))

  # Second pass: update g(), then update coefficients

  # CDF of z:
  Fz_eval = Fz_fun(z = z_grid,
                   weights = rep(1/n, n),
                   mean_vec = mu_z,
                   sd_vec = sigma_z)

  # Check: update the grid if needed
  zcon = contract_grid(z = z_grid,
                       Fz = Fz_eval,
                       lower = 0.001, upper =  0.999)
  z_grid = zcon$z; Fz_eval = zcon$Fz; ngrid = length(z_grid)

  # Transformation:
  g = g_fun(y = y0,
            Fy_eval = Fy_eval,
            z = z_grid,
            Fz_eval = Fz_eval)

  # Update latent data
  z = g(y)

  # updated coefficients (pre-screened)
  theta_hat = Sigma_hat_unscaled%*%crossprod(X[,var_screen], z)

  # Moments of Z|X:
  mu_z = X[,var_screen]%*%theta_hat
  # sigma_z = sqrt(1 + xt_Sigma_hat_unscaled_x) # no need to update
  #----------------------------------------------------------------------------
  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)
  #----------------------------------------------------------------------------
  # Construct the subsets, forcing an intercept:
  X = cbind(1, X); colnames(X)[1] = 'Intercept' # add that intercept

  # All subsets, including an intercept term:
  all_models = all_subsets(1:(p+1))
  colnames(all_models) = colnames(X) # name

  # Remove the ones that *exclude* the intercept:
  all_models = all_models[all_models[,1],]

  # One-time cost:
  log_prior_all = dbinom(0:p, p, prob = prob_inclusion, log = TRUE)
  #----------------------------------------------------------------------------
  # Store MC output:
  post_probs_mc = array(NA, c(nsave, nrow(all_models)))

  # Run the MC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nsave){

    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
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
    #----------------------------------------------------------------------------
    # Block 2: compute the model probabilities

    # Posterior log-probabilities (unnormalized) for *all* models:
    log_probs_unnorm = rep(0, nrow(all_models))
    for(m in 1:nrow(all_models)){

      # Inclusion indicators:
      gamma = 1.0*(all_models[m,])

      # Number of variables (including intercept)
      p_gamma = sum(gamma)

      # Covariates:
      X_gamma = X[,gamma==1,drop=FALSE]

      # Log-likelihood term for this model:
      log_like_m = -p_gamma/2*log(1+psi) - (n/2 + a_sigma)*log(SSR_gprior(z, X_gamma, psi) + 2*b_sigma)

      # Log-prior term for this model:
      log_prior_m = log_prior_all[p_gamma] # dbinom(p_gamma - 1, p, prob = prob_inclusion, log = TRUE)

      # Unnormalized posterior log-probability for this model:
      log_probs_unnorm[m] = log_like_m + log_prior_m
    }
    probs_unnorm = exp(log_probs_unnorm - max(log_probs_unnorm)) # normalize, then exponentiate
    post_probs_mc[nsi,] = probs_unnorm/sum(probs_unnorm)
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
    post_probs = colMeans(post_probs_mc),
    all_models = all_models[,-1], # omit the intercept
    model = 'sblm_modelsel', y = y, X = X[,-1], prob_inclusion = prob_inclusion, psi = psi, fixedX = fixedX, init_screen = init_screen))
}
#----------------------------------------------------------------------------
#' Sample a Gaussian vector using Bhattacharya et al. (2016)
#'
#' Sample from N(mu, Sigma) where Sigma = solve(crossprod(Phi) + solve(D))
#' and mu = Sigma*crossprod(Phi, alpha):
#'
#' @param Phi \code{n x p} matrix (of predictors)
#' @param Ddiag \code{p x 1} vector of diagonal components (of prior variance)
#' @param alpha \code{n x 1} vector (of data, scaled by variance)
#' @return Draw from N(mu, Sigma), which is \code{p x 1}, and is computed in \code{O(n^2*p)}

#' @note Assumes D is diagonal, but extensions are available
#'
#' @references Bhattacharya, Chakraborty, and Mallick (2016, <https://doi.org/10.1093/biomet/asw042>)
#'
#' @export
sampleFastGaussian = function(Phi, Ddiag, alpha){

  # Dimensions:
  Phi = as.matrix(Phi); n = nrow(Phi); p = ncol(Phi)

  # Step 1:
  u = rnorm(n = p, mean = 0, sd = sqrt(Ddiag))
  delta = rnorm(n = n, mean = 0, sd = 1)

  # Step 2:
  v = Phi%*%u + delta

  # Step 3:
  w = solve(crossprod(sqrt(Ddiag)*t(Phi)) + diag(n), #Phi%*%diag(Ddiag)%*%t(Phi) + diag(n)
            alpha - v)

  # Step 4:
  theta =  u + Ddiag*crossprod(Phi, w)

  # Return theta:
  theta
}
#' Numerically stabilize the squared elements
#'
#' Given a vector to be squared, add a numeric buffer
#' for the elements very close to zero.
#'
#' @param vec vector of inputs to be squared
#' @return a vector of the same length as `vec`
#' @importFrom stats mad
square_stabilize = function(vec){
  vec^2 + any(vec^2 < 10^-16)*max(10^-8, mad(vec)/10^6)
}
#' Compute the sum-squared-residuals term under Zellner's g-prior
#'
#' These sum-squared-residuals (SSR) arise in the variance (or precision)
#' term under 1) Zellner's g-prior on the coefficients and a Gamma prior on the
#' error precision and 2) marginalization over the coefficients.
#'
#' @param y vector of response variables
#' @param X matrix of covariates; if NULL, return \code{sum(y^2)}
#' @param psi prior variance (g-prior)
#' @return a positive scalar
#' @importFrom stats lm
SSR_gprior = function(y, X = NULL, psi){

  if(is.null(X) || ncol(X) == 0) return(sum(y^2))

  # Using lm:
  y_hat = fitted(lm(y ~ X - 1))
  SSR = sum(y^2) - psi/(psi+1)*crossprod(y, y_hat)

  # Alternatively:
  #SSR_psi = sum(y^2) - psi/(psi+1)*crossprod(y, X%*%XtXinv%*%crossprod(X, y))

  return(SSR)
}
#' Compute all subsets of a set
#'
#' Given a set of variables, compute the inclusion indicators for
#' all possible subsets.
#'
#' @param set the set from which to compute all subsets (e.g., \code{1:p})
#' @return a data frame where the rows indicate the \code{2^p} different subsets
#' and the columns indicate inclusion (logical) for each element in that subset
#' @references Code adapted from <https://www.r-bloggers.com/2012/04/generating-all-subsets-of-a-set/>
# #' @examples
# #' library(plyr) # required for this function
# #' all_subsets(1:3) # all subsets from {1,2,3}
# #' @export
all_subsets <- function(set) {

  # Library required here:
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop(
      "Package \"plyr\" must be installed to use this function call.",
      call. = FALSE
    )
  }

  bin = expand.grid(plyr::rlply(length(set), c(FALSE, TRUE)))
  as.matrix(attr(plyr::mlply(bin, function(...) { set[c(...)] }), 'split_labels'))
}
