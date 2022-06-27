#evaluating truncated normal density
#Based on mvtnorm and tmvtnorm
pmvnorm_1 <- function (lower = -Inf, upper = Inf, mean = rep(0, length(lower)), 
                       corr = NULL, sigma = NULL, algorithm = GenzBretz(), ...) 
{
  carg <- list(lower = lower, upper = upper, mean = mean, 
               corr = corr, sigma = sigma)
  if (!is.null(carg$corr)) {
    corr <- carg$corr
    if (carg$uni) {
      stop(sQuote("sigma"), " not specified: cannot compute pnorm")
    }
    else {
      lower <- carg$lower - carg$mean
      upper <- carg$upper - carg$mean
      mean <- rep(0, length(lower))
      RET <- mvt(lower = lower, upper = upper, df = 0, 
                 corr = corr, delta = mean, algorithm = algorithm, 
                 ...)
    }
  }
  else {
    if (carg$uni) {
      RET <- list(value = pnorm(carg$upper, mean = carg$mean, 
                                sd = sqrt(carg$sigma)) - pnorm(carg$lower, mean = carg$mean, 
                                                               sd = sqrt(carg$sigma)), error = 0, msg = "univariate: using pnorm")
    }
    else {
      lower <- (carg$lower - carg$mean)/sqrt(diag(carg$sigma))
      upper <- (carg$upper - carg$mean)/sqrt(diag(carg$sigma))
      mean <- rep(0, length(lower))
      corr <- cov2cor(carg$sigma)
      RET <- mvt(lower = lower, upper = upper, df = 0, 
                 corr = corr, delta = mean, algorithm = algorithm, 
                 ...)
    }
  }
  structure(RET$value, error = RET$error, msg = RET$msg)
}


dtmvnorm2 <- function (x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                       lower = rep(-Inf, length = length(mean)), upper = rep(Inf, 
                                                                             length = length(mean)), log = FALSE, margin = NULL) 
{
  cargs <- list(mean = mean, sigma = sigma, lower = lower, 
                upper = upper)
  mean <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  if (!is.null(margin)) {
    if (!length(margin) %in% c(1, 2)) 
      stop("Length of margin must be either 1 (one-dimensional marginal density) or 2 (bivariate marginal density).")
    if (any(margin <= 0) || any(margin > length(mean))) {
      stop("All elements in margin must be in 1..length(mean).")
    }
    if (length(margin) == 1) {
      return(dtmvnorm.marginal(xn = x, n = margin, mean = mean, 
                               sigma = sigma, lower = lower, upper = upper, 
                               log = log))
    }
    if (length(margin) == 2) {
      if (margin[1] == margin[2]) 
        stop("Two different margins needed for bivariate marginal density.")
      if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
      }
      if (!is.matrix(x) || ncol(x) != 2) 
        stop("For bivariate marginal density x must be either a (n x 2) matrix or a vector of length 2.")
      return(dtmvnorm.marginal2(xq = x[, 1], xr = x[, 
                                                    2], q = margin[1], r = margin[2], mean = mean, 
                                sigma = sigma, lower = lower, upper = upper, 
                                log = log))
    }
  }
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  T <- nrow(x)
  insidesupportregion <- logical(T)
  for (i in 1:T) {
    insidesupportregion[i] = all(x[i, ] >= lower & x[i, 
    ] <= upper & !any(is.infinite(x)))
  }
  if (log) {
    dvin <- dmvnorm(x, mean = mean, sigma = sigma, log = TRUE) - 
      log(pmvnorm(lower = lower, upper = upper, mean = mean, 
                  sigma = sigma))
    dvout <- -Inf
  }
  else {
    dvin <- dmvnorm(x, mean = mean, sigma = sigma, log = FALSE)/pmvnorm(lower = lower, 
                                                                        upper = upper, mean = mean, sigma = sigma)
    dvout <- 0
  }
  f <- ifelse(insidesupportregion, dvin, dvout)
  return(f)
}
