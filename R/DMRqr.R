#' Identification of differentially methylated regions using quantile regression
#'
#' @description
#' A function that tests the association between DNA methylation levels
#'    within a predefined region and a phenotype by obtaining p values at each
#'    quantile using a score test and combining them across quantiles with the
#'    Aggregated Cauchy Association Test (ACAT).
#'
#' @param methylation a matrix or numeric vector of DNA methylation levels
#'    (M-values or β-values) for one region. For a matrix, row names
#'    correspond to sample IDs and column names correspond to CpG IDs.
#'    For a numeric vector, values for all (sample, CpG) pairs are stacked into
#'    a single vector of length \eqn{m \times l}, ordered by sample and then
#'    by CpG site, where \eqn{m} is the number of samples and \eqn{l} is
#'    the number of CpG sites.
#' @param phenotype a numeric, character, or factor vector of phenotype values,
#'    with the same number of rows as `methylation` when it is a matrix, or the
#'    same length as `methylation` when it is a vector.
#' @param covariates a numeric matrix of additional covariates, with the same
#'    number of rows or the same length as `methylation`.
#' @param cpg_id a vector of CpG site identifiers, with the same number of rows
#'    or the same length as `methylation`. `cpg_id` can be omitted when
#'    `methylation` is provided as a matrix.
#' @param sample_id a vector of sample identifiers, with the same number of
#'    rows or the same length as `methylation`. `sample_id` can be omitted when
#'    `methylation` is provided as a matrix.
#' @param taus the quantile levels at which the model is fitted. A single value
#'    or a vector of values between 0 and 1.
#'
#' @return
#' If \code{length(taus) == 1}, returns a numeric p value at that quantile.
#' Otherwise, returns an object of S3 class \code{"DMRqr"} with:
#' \describe{
#'    \item{p}{Named numeric vector of p values at each quantile.}
#'    \item{comb.p}{The ACAT-combined p value with equal weights.}
#'    \item{comb.p.weight}{The ACAT-combined p value using quantile-specific
#'    weights (see \code{acat_weight()}).}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1234)
#' m <- 100
#' k <- 3
#' y  <- matrix(rnorm(m * k, mean = 0, sd = 1),nrow = m)
#' x  <- rep(0:1, each = m / 2)
#' taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' z <- matrix(rbinom(m * 2, 1, 0.2), ncol = 2)
#' DMRqr(methylation = y, phenotype = x, covariates = z, taus = taus)
#'
#' set.seed(1234)
#' m <- 100
#' k <- 3
#' n <- m * k
#' y  <- rnorm(n, mean = 0, sd = 1)
#' x  <- rep(rep(0:1, each = m / 2), k)
#' cpg_id <- rep(1:k, each = m)
#' sample_id <- rep(1:m, k)
#' taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' z <- matrix(rbinom(m * 2, 1, 0.2), ncol = 2)
#' z <- z[rep(1: nrow(z), k),]
#' DMRqr(methylation = y, phenotype = x, covariates = z,
#'    cpg_id = cpg_id, sample_id = sample_id, taus = taus)
#'

DMRqr <- function(methylation, phenotype, covariates=NULL,
                  cpg_id, sample_id, taus=c(0.25, 0.5, 0.75)) {

  if (!is.numeric(taus) || anyNA(taus) ||any(taus <= 0) || any(taus >= 1))
    stop("'taus' must contain values between 0 and 1.")
  if (is.matrix(methylation)) {
    nsample <- nrow(methylation)
    ncpg <- ncol(methylation)
    stopifnot(nsample == length(phenotype))
    methylation <- as.vector(methylation)
    sample_id <- rep(1:nsample, ncpg)
    cpg_id <- rep(1:ncpg, each = nsample)
    if (!is.null(covariates)){
      if (!is.matrix(covariates)) stop("'covariates' must be a matrix")
      covariates <- covariates[rep(1:nsample, ncpg),]
    }
    phenotype <- rep(phenotype, ncpg)
  }
  if (!is.numeric(methylation)) stop("'methylation' must be a numeric vector.")

  y = methylation
  x = phenotype
  z = covariates

  cpg_id <- as.factor(cpg_id)
  if (nlevels(cpg_id) <= 1) stop("at least 2 CpGs in a region")

  x_mat <- stats::model.matrix(~ x + x:cpg_id)[, -1L, drop = FALSE]
  if (is.null(z)) {
    z_mat <- stats::model.matrix(~ cpg_id)
  } else {
    z_mat <- stats::model.matrix(~ z + cpg_id)
  }
  taup <- vapply(taus, function(tau) {
    geescore(y = y, x = x_mat, z = z_mat, k = nlevels(cpg_id), subject = sample_id, tau = tau)$pval
  }, numeric(1L))

  taup[taup==1] <- 0.99999
  names(taup) <- taus
  if (length(taus) > 1){
    result <- list(
      p = taup,
      comb.p = ACAT(taup),
      comb.p.weight = ACAT(taup, Weights=acat_weight(taus))
    )
    class(result) = "DMRqr"
    return(result)
  }
  taup
}

#' @noRd
geescore = function(y, x, z, k, subject, tau=0.5)
{

  if (!is.numeric(tau) || length(tau) != 1L ||
      !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("'tau' must be a single numeric value between 0 and 1.")
  }
  n <- length(y)
  if (nrow(x) != n || nrow(z) != n)
    stop("'methylation', 'phenotype' and 'covariates' must have the same length.")
  if (is.factor(subject)) subject <- stats::model.matrix(~ -1 + subject)
  if (is.vector(subject))  subject <- stats::model.matrix(~ -1 + factor(subject))

  p = qr(x)$rank
  q = qr(z)$rank
  an0 <- suppressWarnings(quantreg::rq.fit.br(z, y, tau = tau)$dual)
  bn0 <- an0 - (1 - tau)
  xstar <- as.matrix(stats::lm(x ~ z - 1)$residuals)
  delta0 <- delta_est(an0, subject, tau = tau)

  J <- matrix(1, k, k)
  I <- diag(k)
  r <- as.numeric(delta0-tau^2)/(tau-tau^2)
  Vi <- (tau-tau^2)*((1-r)*I+r*J)
  invVi <- 1/(tau-tau^2)*(J/(k*(1+(k-1)*r))+1/(1-r)*(I-J/k))
  if (r==1) invVi <- MASS::ginv(Vi)
  U <- rep(0, p)
  Vhat <- matrix(0, p, p)
  for (i in 1:ncol(subject)){
    jind <- which(subject[,i]==1)
    U <- U + crossprod(xstar[jind,], invVi%*%bn0[jind])
    Vhat <- Vhat + crossprod(xstar[jind,], invVi%*%xstar[jind,])
  }
  Ts <- crossprod(U, MASS::ginv(Vhat)%*%U)
  pval <- 1 - stats::pchisq(drop(Ts),df = p)

  out = list( pval = pval, delta0 = delta0, U = U, Vhat = Vhat, Ts = Ts,
              n = n)
  return(out)
}

#' @noRd
delta_est = function(an, subject, tau){

  if (!is.matrix(subject)) {
    subject <- as.matrix(subject)
  }
  K <- colSums(subject)
  one_minus_an <- 1 - an

  num <- as.numeric(
    crossprod(one_minus_an, subject) %*%
      crossprod(subject, one_minus_an) - crossprod(one_minus_an)
  )
  denom <- sum(K * (K - 1))

  delta <- num / denom
  min(max(tau^2, delta), tau)
}


#' @noRd
ACAT<-function(Pvals,Weights=NULL){

  if (anyNA(Pvals)){
    stop("Cannot have NAs in the p-values!")
  }

  if (any(Pvals < 0) || any(Pvals > 1)){
    stop("P-values must be between 0 and 1!")
  }

  is.zero<- any(Pvals == 0)
  is.one<- any(Pvals == 1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  if (is.null(Weights)) {
    Weights <- rep.int(1 / length(Pvals), length(Pvals))
  } else if (length(Weights) != length(Pvals)) {
    stop("The length of weights should be the same as that of the p-values")
  } else if (any(Weights < 0)) {
    stop("All the weights must be positive!")
  } else {
    Weights <- Weights / sum(Weights)
  }

  is_small<-(Pvals<1e-16)
  if (!any(is_small)){
    cct_stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct_stat <- sum((Weights[is_small] / Pvals[is_small]) / pi) +
      sum(Weights[!is_small] * tan((0.5 - Pvals[!is_small]) * pi))
  }

  if (cct_stat>1e+15){
    pval<-(1/cct_stat)/pi
  }else{
    pval<-stats::pcauchy(cct_stat, lower.tail=FALSE)
  }
  return(pval)
}

#' @noRd
acat_weight <-function(mytau){
  lw <- ifelse(mytau <= 0.5, mytau, 1 - mytau)
  return(lw/sum(lw))
}


