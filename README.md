# The Variance-Delinkage trade-off

This repo contains scripts to reproduce the figures for the paper "The Shrinkage-Delinkage Trade-off: An Analysis of Factorized Gaussian Approximations for Variational Inference":
* Numerical experiments on Gaussian targets (Figures 1, 2, and 3) can be reproduced using `example_bounds_v3.r`
* Numerical experiments on non-Gaussian targets (Figures 4 and 5) can be reproduced using `model_bounds.R`. The MCMC benchmark and ADVI runs are obtained using Stan (https://mc-stan).
* The script `tools.r` contains functions to compute the upper-bounds derived in the paper (including the bounds in the appendix).

### Required packages

The scripts use standard R packages which can be downloaded from Cran. In addition, we use the packages `cmdstanr` (https://mc-stan.org/cmdstanr/) and `bridgestan` (https://github.com/roualdes/bridgestan). Using `bridgestan` in R sometimes causes memory issues -- the issue has been reported and the developers are working on a fix.
