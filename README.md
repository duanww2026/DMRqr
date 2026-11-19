
# DMRqr

DMRqr is an R package for identifying differentially methylated regions based on quantile regression.

## Installation

``` r
#install.packages("devtools")
devtools::install_github("duanww2026/DMRqr")
```

## Example

``` r
library(DMRqr)
## Methylation values in matrix format
set.seed(1234)
m <- 100
k <- 3
y  <- matrix(rnorm(m * k, mean = 0, sd = 1),nrow = m)
x  <- rep(0:1, each = m / 2)
taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
z <- matrix(rbinom(m * 2, 1, 0.2), ncol = 2)
DMRqr(methylation = y, phenotype = x, covariates = z, taus = taus)


## Methylation values in long format
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

