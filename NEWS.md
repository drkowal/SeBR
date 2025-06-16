# SeBR 1.1.0


## Improvements to previous functionality
* Added `bb()` to sample from the Bayesian bootstrap (BB) posterior more efficiently.
* Added a `fixedX` case for when the covariates are fixed (not random), which also improves computing time for all semiparametric regression functions.
* Since location (intercept) and scale (error standard deviation) are not identifiable in the general transformed regression model, these are no longer reported as coefficients/parameters.
* The posterior draws of the transformation `post_g` now report `(g - intercept)/scale` instead of `g`, which properly corresponds to the transformation under the location-scale identified model. Now, `post_g` can be compared directly to the "true" transformations from simulated data without any further location-scale matching. 

## Fewer dependencies
* `fields` and `GpGp` are only needed for `sbgp()` and `bgp_bc()`.
* `plyr` is only needed for `sblm_modelsel()`.
* `statmod` is only needed for `sbqr()` and `bqr()`.
* `quantreg` is only needed for `sbqr()`.  
* `spikeSlabGAM` is only needed for `sbsm()` and `bsm_bc()`.

## New functions
* Added `sblm_hs()` for semiparametric regression with horseshoe priors. 
* Added `blm_bc_hs()` for Box-Cox transformed regression with horseshoe priors. 
* Added `sblm_ssvs()` for stochastic search variable selection
for semiparametric regression with sparsity priors. 
* Added `sblm_modelsel()` for model/variable selection for semiparametric regression with sparsity priors.
* Added `hbb()` function to sample from the hierarchical BB (HBB) posterior. `concen_hbb()` samples from the marginal posterior distribution of the HBB concentration parameters.
