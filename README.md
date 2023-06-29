SeBR: Semiparametric Bayesian Regression
================

**Overview.** Data transformations are a useful companion for parametric
regression models. A well-chosen or learned transformation can greatly
enhance the applicability of a given model, especially for data with
irregular marginal features (e.g., multimodality, skewness) or various
data domains (e.g., real-valued, positive, or compactly-supported data).

Given paired data $(x_i,y_i)$ for $i=1,\ldots,n$, `SeBR` implements
efficient and fully Bayesian inference for *semiparametric regression
models* that incorporate (1) an unknown data transformation

$$
g(y_i) = z_i
$$

and (2) a useful parametric regression model

$$
z_i  \stackrel{indep}{\sim} P_{Z \mid \theta, X = x_i}
$$

with unknown parameters $\theta$.

**Examples.** We focus on the following important special cases of
$P_{Z \mid \theta, X}$:

1.  The **linear model** is a natural starting point:

$$
z_i = x_i'\theta + \epsilon_i, \quad \epsilon_i \stackrel{iid}{\sim} N(0, \sigma_\epsilon^2)
$$

The transformation $g$ broadens the applicability of this useful class
of models, including for positive or compactly-supported data, while
$P_{Z \mid \theta, X=x} = N(x'\theta, \sigma_\epsilon^2)$.

2.  The **quantile regression model** replaces the Gaussian assumption
    in the linear model with an *asymmetric Laplace* distribution (ALD)

$$
z_i = x_i'\theta + \epsilon_i, \quad \epsilon_i \stackrel{iid}{\sim} ALD(\tau)
$$

to target the $\tau$th quantile of $z$ at $x$, or equivalently, the
$g^{-1}(\tau)$th quantile of $y$ at $x$. The ALD is quite often a very
poor model for real data, especially when $\tau$ is near zero or one.
The transformation $g$ offers a pathway to significantly improve the
model adequacy, while still targeting the desired quantile of the data.

3.  The **Gaussian process (GP) model** generalizes the linear model to
    include a nonparametric regression function,

$$
z_i = f_\theta(x_i) + \epsilon_i, \quad  \epsilon_i \stackrel{iid}{\sim} N(0, \sigma_\epsilon^2)
$$

where $f_\theta$ is a GP and $\theta$ parameterizes the mean and
covariance functions. Although GPs offer substantial flexibility for the
regression function $f_\theta$, this model may be inadequate when $y$
has irregular marginal features or a restricted domain (e.g., positive
or compact).

**Challenges:** The goal is to provide fully Bayesian posterior
inference for the unknowns $(g, \theta)$ and posterior predictive
inference for future/unobserved data $\tilde y(x)$. We prefer a model
and algorithm that offer both (i) flexible modeling of $g$ and (ii)
efficient posterior and predictive computations.

**Innovations:** Our approach (<https://arxiv.org/abs/2306.05498>)
specifies a *nonparametric* model for $g$, yet also provides *Monte
Carlo* (not MCMC) sampling for the posterior and predictive
distributions. As a result, we control the approximation accuracy via
the number of simulations, but do *not* require the lengthy runs,
burn-in periods, convergence diagnostics, or inefficiency factors that
accompany MCMC. The Monte Carlo sampling is typically quite fast.

# Using `SeBR`

The package `SeBR` is installed and loaded as follows:

``` r
# install.packages("devtools")
# devtools::install_github("drkowal/SeBR")
library(SeBR) 
```

The main functions in `SeBR` are:

- `sblm()`: Monte Carlo sampling for posterior and predictive inference
  with the *semiparametric Bayesian linear model*;

- `sbsm()`: Monte Carlo sampling for posterior and predictive inference
  with the *semiparametric Bayesian spline model*, which replaces the
  linear model with a spline for nonlinear modeling of
  $x \in \mathbb{R}$;

- `sbqr()`: blocked Gibbs sampling for posterior and predictive
  inference with the *semiparametric Bayesian quantile regression*; and

- `sbgp()`: Monte Carlo sampling for predictive inference with the
  *semiparametric Bayesian Gaussian process model*.

Each function returns a point estimate of $\theta$ (`coefficients`),
point predictions at some specified testing points (`fitted.values`),
posterior samples of the transformation $g$ (`post_g`), and posterior
predictive samples of $\tilde y(x)$ at the testing points
(`post_ypred`), as well as other function-specific quantities (e.g.,
posterior draws of $\theta$, `post_theta`). The calls `coef()` and
`fitted()` extract the point estimates and point predictions,
respectively.

**Note:** The package also includes Box-Cox variants of these functions,
i.e., restricting $g$ to the (signed) Box-Cox parametric family
$g(t; \lambda) = \{\mbox{sign}(t) \vert t \vert^\lambda - 1\}/\lambda$
with known or unknown $\lambda$. The parametric transformation is less
flexible, especially for irregular marginals or restricted domains, and
requires MCMC sampling. These functions (e.g., `blm_bc()`, etc.) are
primarily for benchmarking.

Detailed documentation and examples are available at
<https://drkowal.github.io/SeBR/>.
