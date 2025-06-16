#----------------------------------------------------------------------------
#' Simulate a transformed linear model
#'
#' Generate training data (X, y) and testing data (X_test, y_test)
#' for a transformed linear model. The covariates are correlated
#' Gaussian variables. A user-specified proportion (\code{prop_sig})
#' of the regression coefficients are nonozero (= 1) and the rest are zero.
#' There are multiple options for the transformation, which define the support
#' of the data (see below).
#'
#' @param n number of observations in the training data
#' @param p number of covariates
#' @param g_type type of transformation; must be one of
#' \code{beta}, \code{step}, or \code{box-cox}
#' @param n_test number of observations in the testing data
#' @param heterosked logical; if TRUE, simulate the latent data with heteroskedasticity
#' @param lambda Box-Cox parameter (only applies for \code{g_type = 'box-cox'})
#' @param prop_sig proportion of signals (nonzero coefficients)
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
#' @note The design matrices \code{X} and \code{X_test}
#' do not include an intercept and there is no
#' intercept parameter in \code{beta_true}. The
#' location/scale of the data are not identified
#' in general transformed regression models, so
#' recovering them is not a goal.
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
                        lambda = 1,
                        prop_sig = 0.5){
  #----------------------------------------------------------------------------
  # Checks:
  if(!is.element(g_type, c("beta", "step", "box-cox")))
    stop("The transformation must be one of 'beta', 'step', or 'box-cox'")

  if(p <= 1)
    stop("Include at least p=2 covariates")
  #----------------------------------------------------------------------------
  # Simulate a design matrix with correlated predictors:
  ar1 = 0.75
  X = t(apply(matrix(0, nrow = n, ncol = p), 1, function(x)
              arima.sim(n = p, list(ar = ar1), sd = sqrt(1-ar1^2))))

  # Shuffle the columns:
  ind_shuff = sample(1:p)
  X = X[,ind_shuff]

  # Covariates on the testing data:
  X_test = t(apply(matrix(0, nrow = n_test, ncol = p), 1, function(x)
                   arima.sim(n = p, list(ar = ar1), sd = sqrt(1-ar1^2))))
  # Match the column shuffling:
  X_test = X_test[,ind_shuff]
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
  beta_true = c(rep(1, ceiling(p*prop_sig)),
                rep(0, p - ceiling(p*prop_sig)))
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
#' Bayesian bootstrap posterior sampler for the CDF
#'
#' Compute one Monte Carlo draw from the Bayesian bootstrap (BB)
#' posterior distribution of the cumulative distribution function (CDF).
#'
#' @param y the data from which to infer the CDF (preferably sorted)
#' @return a function that can evaluate the sampled CDF at any argument(s)
#'
#' @details Assuming the data \code{y} are iid from an unknown distribution,
#' the Bayesian bootstrap (BB) is a nonparametric model for this distribution. The
#' BB is a limiting case of a Dirichlet process prior (without
#' any hyperparameters) that admits direct Monte Carlo (not MCMC) sampling.
#'
#' This function computes one draw from the BB posterior
#' distribution for the CDF \code{Fy}.
#'
#' @note This code is inspired by \code{ggdist::weighted_ecdf}.
#'
#' @examples
#' # Simulate data:
#' y = rnorm(n = 100)
#'
#' # One draw from the BB posterior:
#' Fy = bb(y)
#'
#' class(Fy) # this is a function
#' Fy(0) # some example use (for this one draw)
#' Fy(c(.5, 1.2))
#'
#' # Plot several draws from the BB posterior distribution:
#' ys = seq(-3, 3, length.out=1000)
#' plot(ys, ys, type='n', ylim = c(0,1),
#'      main = 'Draws from BB posterior', xlab = 'y', ylab = 'F(y)')
#' for(s in 1:50) lines(ys, bb(y)(ys), col='gray')
#'
#' # Add ECDF for reference:
#' lines(ys, ecdf(y)(ys), lty=2)
#'
#' @importFrom stats rgamma approxfun
#' @export
bb = function(y){

  # Length of data:
  n = length(y)

  # Sort the y's, if unsorted
  if(is.unsorted(y)) y = y[order(y)]

  # Dirichlet(1) weights:
  weights_y = rgamma(n = n, shape = 1)
  weights_y  = weights_y/sum(weights_y)

  # Key input (with rescaling by n/(n+1) as in the paper for boundary reasons)
  sum_weights_y = n/(n+1)*cumsum(weights_y)

  # Use approxfun() for fast computing (w/ n/(n+1) rescaling as above)
  Fy = approxfun(y, sum_weights_y,
                 yleft = 0, yright = n/(n+1),
                 ties = "ordered",
                 method = "constant")
  return(Fy)
}
#' Hierarchical Bayesian bootstrap posterior sampler
#'
#' Compute one Monte Carlo draw from the hierarchical Bayesian bootstrap (HBB)
#' posterior distribution of the cumulative distribution function (CDF) for
#' each group. The common (BB) and group-specific (HBB) weights are also returned.
#'
#' @param y the data from which to infer the group-specific CDFs
#' @param groups the group assignment for each element of \code{y}
#' @param sample_alphas logical; if TRUE, sample the concentration hyperparameters
#' from their marginal posterior distribution
#' @param shape_alphas (optional) shape parameter for the Gamma prior on each \code{alphas} (if sampled)
#' @param rate_alphas (optional) rate parameter for the Gamma prior on each \code{alphas} (if sampled)
#' @param alphas (optional) vector of fixed concentration hyperparameters corresponding
#' to the unique levels in \code{groups} (used when \code{sample_alphas = FALSE})
#' @param M a positive scaling term to set a default value of \code{alphas} when
#' it is unspecified (\code{alphas = NULL}) and not sampled (\code{sample_alphas = FALSE})
#'
#' @return a list with the following elements:
#' \itemize{
#'  \item \code{Fyc}: a list of functions where each entry corresponds to a group
#'  and that group-specific function can evaluate the sampled CDF at any argument(s)
#'  \item \code{weights_y}: sampled weights from the common (BB) distribution (\code{n}-dimensional)
#'  \item \code{weights_yc}: sampled weights from each of the \code{K} groups (\code{K x n})
#'  \item \code{alphas}: the (fixed or sampled) concentration hyperparameters
#' }
#'
#' @details Assuming the data \code{y} are independent with unknown,
#' group-specific distributions, the hierarchical Bayesian bootstrap (HBB) from
#' Oganisian et al. (<https://doi.org/10.1515/ijb-2022-0051>) is a nonparametric model
#' for each distribution. The HBB includes hierarchical shrinkage across these
#' groups toward a common distribution (the \code{\link{bb}}). The HBB admits
#' direct Monte Carlo (not MCMC) sampling.
#'
#' The shrinkage toward this common distribution is determined by the concentration
#' hyperparameters \code{alphas}. Each component of \code{alphas} corresponds to
#' one of the groups. Larger values encourage more shrinkage toward
#' the common distribution, while smaller values allow more substantial deviations for that group.
#'
#' When \code{sample_alphas=TRUE}, each component of \code{alphas} is sampled from its marginal
#' posterior distribution, assuming independent Gamma(\code{shape_alphas}, \code{rate_alphas})
#' priors. This step uses a simple grid approximation to enable efficient sampling that
#' preserves joint Monte Carlo sampling with the group-specific and common distributions.
#' See \code{\link{concen_hbb}} for details. Note that diffuse priors on \code{alphas}
#' tends to produce more aggressive shrinkage toward the common distribution (complete pooling).
#' For moderate shrinkage, we use the default values \code{shape_alphas = 30*K} and \code{rate_alphas = 1}
#' where \code{K} is the number of groups.
#'
#' When \code{sample_alphas=FALSE}, these concentration hyperparameters are fixed
#' at user-specified values. That can be done by specifying \code{alphas} directly.
#' Alternatively, if \code{alphas} is left unspecified (\code{alphas = NULL}),
#' we adopt the default from Oganisian et al. which sets the \code{c}th entry to \code{M*n/nc}
#' where \code{M} is user-specified and \code{nc} is the number of observations in group \code{c}.
#' For further guidance on the choice of \code{M}:
#' \itemize{
#' \item \code{M = 0.01/K} approximates separate BB's by group (no pooling);
#' \item \code{M} between 10 and 100 gives moderate shrinkage (partial pooling); and
#' \item \code{M = 100*max(nc)} approximates a common BB (complete pooling).
#' }
#'
#' @note If supplying \code{alphas} with distinct entries, make sure that the
#' groups are ordered properly; these entries should match \code{sort(unique(groups))}.
#'
#' @references Oganisian et al. (<https://doi.org/10.1515/ijb-2022-0051>)
#'
#' @examples
#' # Sample size and number of groups:
#' n = 500
#' K = 3
#'
#' # Define the groups, then assign:
#' ugroups = paste('g', 1:K, sep='') # groups
#' groups = sample(ugroups, n, replace = TRUE) # assignments
#'
#' # Simulate the data: iid normal, then add group-specific features
#' y = rnorm(n = n) # data
#' for(g in ugroups)
#'   y[groups==g] = y[groups==g] + 3*rnorm(1) # group-specific
#'
#' # One draw from the HBB posterior of the CDF:
#' samp_hbb = hbb(y, groups)
#'
#' names(samp_hbb) # items returned
#' Fyc = samp_hbb$Fyc # list of CDFs
#' class(Fyc) # this is a list
#' class(Fyc[[1]]) # each element is a function
#'
#' c = 1 # try: vary in 1:K
#' Fyc[[c]](0) # some example use (for this one draw)
#' Fyc[[c]](c(.5, 1.2))
#'
#' # Plot several draws from the HBB posterior distribution:
#' ys = seq(min(y), max(y), length.out=1000)
#' plot(ys, ys, type='n', ylim = c(0,1),
#'      main = 'Draws from HBB posteriors', xlab = 'y', ylab = 'F_c(y)')
#' for(s in 1:50){ # some draws
#'
#'   # BB CDF:
#'   Fy = bb(y)
#'   lines(ys, Fy(ys), lwd=3) # plot CDF
#'
#'   # HBB:
#'   Fyc = hbb(y, groups)$Fyc
#'
#'   # Plot CDFs by group:
#'   for(c in 1:K) lines(ys, Fyc[[c]](ys), col=c+1, lwd=3)
#' }
#'
#' # For reference, add the ECDFs by group:
#' for(c in 1:K) lines(ys, ecdf(y[groups==ugroups[c]])(ys), lty=2)
#'
#' legend('bottomright', c('BB', paste('HBB:', ugroups)), col = 1:(K+1), lwd=3)
#'
#' @importFrom stats rgamma approxfun
#' @export
hbb = function(y, groups,
               sample_alphas = FALSE,
               shape_alphas = NULL,
               rate_alphas = NULL,
               alphas = NULL, M = 30){

  #----------------------------------------------------------------------------
  # Start with checks and basic definitions

  # Length of data:
  n = length(y)

  # Check: lengths:
  if(length(groups) != n)
    stop('y and groups must have the same length')

  # Unique groups (levels) and their number:
  ugroups = sort(unique(groups))
  K = length(ugroups)

  # Sort the y's (and groups), if unsorted
  if(is.unsorted(y)) {
    ord = order(y)
    y = y[ord]
    groups = groups[ord]
  }

  #----------------------------------------------------------------------------
  # Step 0: sample or fix the concentration hyperparameters
  if(sample_alphas){

    # Default values:
    if(is.null(shape_alphas)) shape_alphas = 30*K
    if(is.null(rate_alphas)) rate_alphas = 1

    # One draw of alphas:
    alphas = concen_hbb(groups = groups,
                           shape_alphas = shape_alphas,
                           rate_alphas = rate_alphas,
                           nsave = 1,
                           ngrid = 500)

  } else {

    # When hyperparameters are unspecified, use the default from Organisian et al:
    if(is.null(alphas)){
      nc = table(groups)[ugroups] # number of observations by group
      alphas = M*n/nc
    }

    # Check whether the alphas have the right length and
    if(length(alphas) != K)
      stop('alphas must have length equal to the number of unique groups')
  }
  #----------------------------------------------------------------------------
  # Step 1: sample the "global" CDF (using BB)

  # Dirichlet(1) weights:
  weights_y = rgamma(n = n, shape = 1)
  weights_y  = weights_y/sum(weights_y)

  # Key input (with rescaling by n/(n+1) as in the paper for boundary reasons)
  sum_weights_y = n/(n+1)*cumsum(weights_y)

  # Use approxfun() for fast computing (w/ n/(n+1) rescaling as above)
  #   Note: this is never used, so no need to compute it...
  # Fy0 = approxfun(y, sum_weights_y,
  #                 yleft = 0, yright = n/(n+1),
  #                 ties = "ordered",
  #                 method = "constant")
  #----------------------------------------------------------------------------
  # Step 2: sample the group-specific CDFs
  Fyc = vector('list', K) # list of functions...
  names(Fyc) = ugroups
  weights_yc = matrix(NA, nrow = K, ncol = n) # to store the sampled weights

  # For each group:
  for(c in 1:K){

    # Dirichlet(ac) weights:
    ac = alphas[c]*weights_y + I(groups==ugroups[c])
    weights_yc[c,] = rgamma(n = n, shape = ac)
    weights_yc[c,]  = weights_yc[c,]/sum(weights_yc[c,])

    # Key input (with rescaling by n/(n+1) as in the paper for boundary reasons)
    sum_weights_yc = n/(n+1)*cumsum(weights_yc[c,])

    # Use approxfun() for fast computing (w/ n/(n+1) rescaling as above)
    Fyc[[c]] = approxfun(y, sum_weights_yc,
                         yleft = 0, yright = n/(n+1),
                         ties = "ordered",
                         method = "constant")
  }

  return(list(
    #Fy0 = Fy0, # not needed
    Fyc = Fyc,
    weights_y = weights_y,
    weights_yc = weights_yc,
    alphas = alphas
  ))
}
#' Posterior sampling algorithm for the HBB concentration hyperparameters
#'
#' Compute Monte Carlo draws from the (marginal) posterior distribution of the
#' concentration hyperparameters of the hierarchical Bayesian bootstrap
#' (\code{\link{hbb}}). The HBB is a nonparametric model for group-specific
#' distributions; each group has a concentration parameter, where
#' larger values encourage more shrinkage toward a common distribution.
#'
#' @param groups the group assignments in the observed data
#' @param shape_alphas (optional) shape parameter for the Gamma prior
#' @param rate_alphas (optional) rate parameter for the Gamma prior
#' @param nsave (optional) number of Monte Carlo simulations
#' @param ngrid (optional) number of grid points
#' @return \code{nsave x K} samples of the concentration hyperparameters
#' corresponding to the \code{K} groups
#'
#' @details The concentration hyperparameters are assigned
#' independent Gamma(\code{shape_alphas}, \code{rate_alphas}) priors.
#' This function uses a grid approximation to the marginal posterior
#' with the goal of producing a simple algorithm. Because this is a
#' *marginal* posterior sampler, it can be used with the \code{\link{hbb}}
#' sampler (which conditions on \code{alphas}) to provide a joint
#' Monte Carlo (not MCMC) sampling algorithm for the concentration
#' hyperparameters, the group-specific CDFs, and the common CDF.

#' Note that diffuse priors on \code{alphas} tend to put posterior mass on
#' large values, which leads to more aggressive shrinkage toward the common distribution
#' (complete pooling). For moderate shrinkage, we use the default values
#' \code{shape_alphas = 30*K} and \code{rate_alphas = 1}, where \code{K} is the
#' number of groups.
#'
#' @references Oganisian et al. (<https://doi.org/10.1515/ijb-2022-0051>)
#'
#' @examples
#' # Dimensions:
#' n = 500 # number of observations
#' K = 3 # number of groups
#'
#' # Assign groups w/ unequal probabilities:
#' ugroups = paste('g', 1:K, sep='') # groups
#' groups = sample(ugroups,
#'                 size = n,
#'                 replace = TRUE,
#'                 prob = 1:K) # unequally weighted (unnormalized)
#'
#' # Summarize:
#' table(groups)/n
#'
#' # Marginal posterior sampling for alpha:
#' post_alpha = concen_hbb(groups)
#'
#' # Summarize: posterior distributions
#' for(c in 1:K) {
#'   hist(post_alpha[,c],
#'        main = paste("Concentration parameter: group", ugroups[c]),
#'        xlim = range(post_alpha))
#'   abline(v = mean(post_alpha[,c]), lwd=3) # posterior mean
#' }
#'
#' @importFrom stats dgamma
#' @export
concen_hbb = function(groups,
                         shape_alphas = NULL,
                         rate_alphas = NULL,
                         nsave = 1000,
                         ngrid = 500){

  # Unique groups (levels) and their number:
  ugroups = sort(unique(groups))
  K = length(ugroups)

  # Number of observations by group:
  nc = table(groups)[ugroups]

  # Total number of observations:
  n = length(groups)

  # Default values:
  if(is.null(shape_alphas)) shape_alphas = 30*K
  if(is.null(rate_alphas)) rate_alphas = 1

  # Establish a grid
  #   alpha = 0.01 ~ separate BBs by group
  #   alpha = 100*n ~ common BB
  # Use a regular grid for "small" values, then log-grid for "large" values
  a_grid = c(seq(0.01, 10, length.out = ceiling(ngrid/5)),
             exp(seq(log(10.01), log(100*n), length.out = ngrid - ceiling(ngrid/5))))

  # Log-likelihood for each group's alpha:
  #   NOTE: use recurring terms instead to save compute time
  # loglike_alpha = function(alphac, nc) lgamma(alphac) + nc*log(alphac) - lgamma(alphac + nc)

  # Recurring terms shared for all c = 1:K (and all simulations)
  lga = lgamma(a_grid)
  log_prior = dgamma(a_grid,
                     shape = shape_alphas,
                     rate = rate_alphas,
                     log = TRUE)

  # Storage:
  post_alphas = matrix(NA, nrow = nsave, ncol = K); colnames(post_alphas) = ugroups

  # For each group:
  for(c in 1:K){

    # Log-posterior (up to a constant) on the grid:
    log_post = log_prior +
      lga + nc[c]*log(a_grid) - lgamma(a_grid + nc[c]) # log-likelihood on grid

    # Add a constant to avoid (some) numerical issues:
    log_post = log_post + abs(max(log_post))

    # Also exclude very low probability grid points:
    sub_c = which(exp(log_post) > 10^-16)
    post_alphas[, c] = sample(x = a_grid[sub_c],
                              size = nsave,
                              replace = TRUE,
                              prob = exp(log_post[sub_c]))
  }

  # For a single draw, convert to a vector
  if(nsave==1) post_alphas = post_alphas[1,]

  return(post_alphas)
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
#' @examples
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
#' Estimate the remaining time in the algorithm
#' @param nsi current iteration
#' @param timer0 initial timer value from \code{proc.time()[3]}
#' @param nsims total number of simulations
#' @param nprints total number of printed updates
#' @return estimate of remaining time
computeTimeRemaining = function(nsi, timer0, nsims, nprints = 2){

  # Print every nrep:
  nrep = ceiling(nsims/(nprints +1)) + 1

  # Only print occasionally:
  if(nsi%%nrep == 0){ # || nsi==ninit) {
    # Current time:
    timer = proc.time()[3]

    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)

    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec

    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hr remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(ceiling(secRemaining/60), "min remaining"))
      } else print(paste(ceiling(secRemaining), "sec remaining"))
    }
  }
}
