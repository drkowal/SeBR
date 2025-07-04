% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source_sba.R
\name{sbgp}
\alias{sbgp}
\title{Semiparametric Bayesian Gaussian processes}
\usage{
sbgp(
  y,
  locs,
  X = NULL,
  covfun_name = "matern_isotropic",
  locs_test = locs,
  X_test = NULL,
  nn = 30,
  emp_bayes = TRUE,
  fixedX = (length(y) >= 500),
  approx_g = FALSE,
  nsave = 1000,
  ngrid = 100
)
}
\arguments{
\item{y}{\code{n x 1} response vector}

\item{locs}{\code{n x d} matrix of locations}

\item{X}{\code{n x p} design matrix; if unspecified, use intercept only}

\item{covfun_name}{string name of a covariance function; see ?GpGp}

\item{locs_test}{\code{n_test x d} matrix of locations
at which predictions are needed; default is \code{locs}}

\item{X_test}{\code{n_test x p} design matrix for test data;
default is \code{X}}

\item{nn}{number of nearest neighbors to use; default is 30
(larger values improve the approximation but increase computing cost)}

\item{emp_bayes}{logical; if TRUE, use a (faster!) empirical Bayes
approach for estimating the mean function}

\item{fixedX}{logical; if TRUE, treat the design as fixed (non-random) when sampling
the transformation; otherwise treat covariates as random with an unknown distribution}

\item{approx_g}{logical; if TRUE, apply large-sample
approximation for the transformation}

\item{nsave}{number of Monte Carlo simulations}

\item{ngrid}{number of grid points for inverse approximations}
}
\value{
a list with the following elements:
\itemize{
\item \code{coefficients} the estimated regression coefficients
\item \code{fitted.values} the posterior predictive mean at the test points \code{locs_test}
\item \code{fit_gp} the fitted \code{GpGp_fit} object, which includes
covariance parameter estimates and other model information
\item \code{post_ypred}: \code{nsave x ntest} samples
from the posterior predictive distribution at \code{locs_test}
\item \code{post_g}: \code{nsave} posterior samples of the transformation
evaluated at the unique \code{y} values
\item \code{model}: the model fit (here, \code{sbgp})
}
as well as the arguments passed in.
}
\description{
Monte Carlo sampling for Bayesian Gaussian process regression with an
unknown (nonparametric) transformation.
}
\details{
This function provides Bayesian inference for a
transformed Gaussian process model using Monte Carlo (not MCMC) sampling.
The transformation is modeled as unknown and learned jointly
with the regression function (unless \code{approx_g = TRUE}, which then uses
a point approximation). This model applies for real-valued data, positive data, and
compactly-supported data (the support is automatically deduced from the observed \code{y} values).
The results are typically unchanged whether \code{laplace_approx} is TRUE/FALSE;
setting it to TRUE may reduce sensitivity to the prior, while setting it to FALSE
may speed up computations for very large datasets. For computational efficiency,
the Gaussian process parameters are fixed at point estimates, and the latent Gaussian
process is only sampled when \code{emp_bayes = FALSE}. However, the uncertainty
from this term is often negligible compared to the observation errors, and the
transformation serves as an additional layer of robustness. By default, \code{fixedX} is
set to FALSE for smaller datasets (\code{n < 500}) and TRUE for larger datasets (\code{n >= 500}).
}
\examples{
\donttest{
# Simulate some data:
n = 200 # sample size
x = seq(0, 1, length = n) # observation points

# Transform a noisy, periodic function:
y = g_inv_bc(
  sin(2*pi*x) + sin(4*pi*x) + rnorm(n),
             lambda = .5) # Signed square-root transformation

# Package we use for fast computing w/ Gaussian processes:
library(GpGp)

# Fit the semiparametric Bayesian Gaussian process:
fit = sbgp(y = y, locs = x)
names(fit) # what is returned
coef(fit) # estimated regression coefficients (here, just an intercept)
class(fit$fit_gp) # the GpGp object is also returned

# Plot the model predictions (point and interval estimates):
pi_y = t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90\% PI
plot(x, y, type='n', ylim = range(pi_y,y),
     xlab = 'x', ylab = 'y', main = paste('Fitted values and prediction intervals'))
polygon(c(x, rev(x)),c(pi_y[,2], rev(pi_y[,1])),col='gray', border=NA)
lines(x, y, type='p') # observed points
lines(x, fitted(fit), lwd = 3) # fitted curve
}
}
