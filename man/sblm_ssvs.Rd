% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source_varsel.R
\name{sblm_ssvs}
\alias{sblm_ssvs}
\title{Semiparametric Bayesian linear model with stochastic search variable selection}
\usage{
sblm_ssvs(
  y,
  X,
  X_test = X,
  psi = length(y),
  fixedX = (length(y) >= 500),
  approx_g = FALSE,
  init_screen = NULL,
  a_pi = 1,
  b_pi = 1,
  nsave = 1000,
  nburn = 1000,
  ngrid = 100,
  verbose = TRUE
)
}
\arguments{
\item{y}{\code{n x 1} response vector}

\item{X}{\code{n x p} matrix of predictors (no intercept)}

\item{X_test}{\code{n_test x p} matrix of predictors for test data;
default is the observed covariates \code{X}}

\item{psi}{prior variance (g-prior)}

\item{fixedX}{logical; if TRUE, treat the design as fixed (non-random) when sampling
the transformation; otherwise treat covariates as random with an unknown distribution}

\item{approx_g}{logical; if TRUE, apply large-sample
approximation for the transformation}

\item{init_screen}{for the initial approximation, number of covariates
to pre-screen (necessary when \code{p > n}); if NULL, use \code{n/log(n)}}

\item{a_pi}{shape1 parameter of the (Beta) prior inclusion probability}

\item{b_pi}{shape2 parameter of the (Beta) prior inclusion probability}

\item{nsave}{number of MCMC simulations to save}

\item{nburn}{number of MCMC iterations to discard}

\item{ngrid}{number of grid points for inverse approximations}

\item{verbose}{logical; if TRUE, print time remaining}
}
\value{
a list with the following elements:
\itemize{
\item \code{coefficients} the posterior mean of the regression coefficients
\item \code{fitted.values} the posterior predictive mean at the test points \code{X_test}
\item \code{selected}: the variables (columns of \code{X}) selected by the median probability model
\item \code{pip}: (marginal) posterior inclusion probabilities for each variable
\item \code{post_theta}: \code{nsave x p} samples from the posterior distribution
of the regression coefficients
\item \code{post_gamma}: \code{nsave x p} samples from the posterior distribution
of the variable inclusion indicators
\item \code{post_ypred}: \code{nsave x n_test} samples
from the posterior predictive distribution at test points \code{X_test}
\item \code{post_g}: \code{nsave} posterior samples of the transformation
evaluated at the unique \code{y} values
\item \code{model}: the model fit (here, \code{sblm_ssvs})
}
as well as the arguments passed in.
}
\description{
MCMC sampling for semiparametric Bayesian linear regression with
1) an unknown (nonparametric) transformation and 2) a sparsity prior on
the (possibly high-dimensional) regression coefficients. Here, unlike \code{\link{sblm}},
Gibbs sampling is used for the variable inclusion indicator variables
\code{gamma}, referred to as stochastic search variable selection (SSVS).
All remaining terms--including the transformation \code{g}, the regression
coefficients \code{theta}, and any predictive draws--are drawn directly from
the joint posterior (predictive) distribution.
}
\details{
This function provides fully Bayesian inference for a
transformed linear model with sparse g-priors on the regression coefficients.
The transformation is modeled as unknown and learned jointly
with the regression coefficients (unless \code{approx_g} = TRUE, which then uses
a point approximation). This model applies for real-valued data, positive data, and
compactly-supported data (the support is automatically deduced from the observed \code{y} values).
By default, \code{fixedX} is set to \code{FALSE} for smaller datasets (\code{n < 500})
and \code{TRUE} for larger datasets.

The sparsity prior is especially useful for variable selection. Compared
to the horseshoe prior version (\code{\link{sblm_hs}}), the sparse g-prior
is advantageous because 1) it truly allows for sparse (i.e., exactly zero)
coefficients in the prior and posterior, 2) it incorporates covariate
dependencies via the g-prior structure, and 3) it tends to perform well
under both sparse and non-sparse regimes, while the horseshoe version only
performs well under sparse regimes. The disadvantage is that
SSVS does not scale nearly as well in \code{p}.

Following Scott and Berger (<https://doi.org/10.1214/10-AOS792>),
we include a \code{Beta(a_pi, b_pi)} prior on the prior inclusion probability. This term
is then sampled with the variable inclusion indicators \code{gamma} in a
Gibbs sampling block. All other terms are sampled using direct Monte Carlo
(not MCMC) sampling.

Alternatively, model probabilities can be computed directly
(by Monte Carlo, not MCMC/Gibbs sampling) using \code{\link{sblm_modelsel}}.
}
\note{
The location (intercept) and scale (\code{sigma_epsilon}) are
not identified, so any intercepts in \code{X} and \code{X_test} will
be removed. The model-fitting *does* include an internal location-scale
adjustment, but the function only outputs inferential summaries for the
identifiable parameters.
}
\examples{
\donttest{
# Simulate data from a transformed (sparse) linear model:
dat = simulate_tlm(n = 100, p = 15, g_type = 'step')
y = dat$y; X = dat$X # training data
y_test = dat$y_test; X_test = dat$X_test # testing data

hist(y, breaks = 25) # marginal distribution

# Fit the semiparametric Bayesian linear model with sparsity priors:
fit = sblm_ssvs(y = y, X = X, X_test = X_test)
names(fit) # what is returned

# Evaluate posterior predictive means and intervals on the testing data:
plot_pptest(fit$post_ypred, y_test,
            alpha_level = 0.10) # coverage should be about 90\%

# Check: correlation with true coefficients
cor(dat$beta_true, coef(fit))

# Selected coefficients under median probability model:
fit$selected

# True signals:
which(dat$beta_true != 0)

# Summarize the transformation:
y0 = sort(unique(y)) # posterior draws of g are evaluated at the unique y observations
plot(y0, fit$post_g[1,], type='n', ylim = range(fit$post_g),
     xlab = 'y', ylab = 'g(y)', main = "Posterior draws of the transformation")
temp = sapply(1:nrow(fit$post_g), function(s)
  lines(y0, fit$post_g[s,], col='gray')) # posterior draws
lines(y0, colMeans(fit$post_g), lwd = 3) # posterior mean
lines(y, dat$g_true, type='p', pch=2) # true transformation

# Posterior predictive checks on testing data: empirical CDF
y0 = sort(unique(y_test))
plot(y0, y0, type='n', ylim = c(0,1),
     xlab='y', ylab='F_y', main = 'Posterior predictive ECDF')
temp = sapply(1:nrow(fit$post_ypred), function(s)
  lines(y0, ecdf(fit$post_ypred[s,])(y0), # ECDF of posterior predictive draws
        col='gray', type ='s'))
lines(y0, ecdf(y_test)(y0),  # ECDF of testing data
     col='black', type = 's', lwd = 3)
}

}
