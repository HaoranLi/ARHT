
<!-- README.md is generated from README.Rmd. Please edit that file -->
ARHT
====

Perform the Adaptable Regularized Hotelling's *T*<sup>2</sup> test (ARHT) proposed by Li et al. (2016). Both one- and two- sample mean test are available with various probabilistic alternative prior models. It contains a function to consistently estimate higher order moments of the population covariance spectral distribution using the spectral of the sample covariance matrix. In addition, it contains a function to sample from 3-variate chi-squared random vectors approximately with a given correlation matrix when the degrees of freedom are large.

Installation
------------

You can install ARHT from github with:

``` r
# install.packages("devtools")
devtools::install_github("HaoranLi/ARHT")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(ARHT)
## basic example code
set.seed(10086)
# One-sample test
n1 = 300; p =500
dataX = matrix(rnorm(n1 * p), nrow = n1, ncol = p)
res1 = ARHT(dataX)

# Two-sample test
n2= 400
dataY = matrix(rnorm(n2 * p), nrow = n2, ncol = p )
res2 = ARHT(dataX, dataY, mu_0 = rep(0.01,p))

# Specify probabilistic alternative priors model
res3 = ARHT(dataX, dataY, mu_0 = rep(0.01,p), 
            prob_alt_prior = list(c(1/3, 1/3, 1/3), c(0,1,0)))
# Change Type 1 error calibration method
res4 = ARHT(dataX, dataY, mu_0 = rep(0.01,p),
            Type1error_calib = "sqrt")
RejectOrNot = res4$ARHT_pvalue < 0.05
```

Reference
---------

Li, Haoran, Alexander Aue, Debashis Paul, Jie Peng, and Pei Wang. 2016. “An Adaptable Generalization of Hotelling's *T*<sup>2</sup> Test in High Dimension.” arXiv preprint arXiv:1609.08725.
