
# DMRqr

<!-- badges: start -->
<!-- badges: end -->

The goal of DMRqr is to identify differentially methylated regions using quantile regression.

## Installation

You can install the development version of DMRqr like so:

``` r
#install.packages("pak")
pak::pak("duanww2026/DMRqr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(DMRqr)
set.seed(1234)
m <- 100
k <- 3
y  <- matrix(rnorm(m * k, mean = 0, sd = 1),nrow = m)
x  <- rep(0:1, each = m / 2)
taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
z <- matrix(rbinom(m * 2, 1, 0.2), ncol = 2)
DMRqr(methylation = y, phenotype = x, covariates = z, taus = taus)

set.seed(1234)
m <- 100
k <- 3
n <- m * k
y  <- rnorm(n, mean = 0, sd = 1)
x  <- rep(rep(0:1, each = m / 2), k)
cpg_id <- rep(1:k, each = m)
sample_id <- rep(1:m, k)
taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
z <- matrix(rbinom(m * 2, 1, 0.2), ncol = 2)
z <- z[rep(1: nrow(z), k),]
DMRqr(methylation = y, phenotype = x, covariates = z,
   cpg_id = cpg_id, sample_id = sample_id, taus = taus)
```

