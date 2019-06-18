dtobit <- function(y, xTbeta, sigma, left = -1, log = FALSE)
{
  d <- y > left

  if (log)
  {
    return(d * (dnorm((y - xTbeta)/sigma, log = TRUE) - log(sigma)) +
             (1 - d) * pnorm((left - xTbeta)/sigma, log.p = TRUE))
  } else
  {
    return((1/sigma * dnorm((y - xTbeta)/sigma))^d * pnorm((xTbeta - left)/sigma, lower.tail = FALSE)^(1-d))
  }
}

objective <- function(theta, y, X, left)
{
  beta <- c(1, theta[1:ncol(X)])
  xTbeta <- cbind(1, X) %*% beta
  gamma <- theta[(ncol(X)+1):length(theta)]

  exp.terms <- sapply(1:length(gamma), function(i){
    gamma[i] * xTbeta^(i - 1)
  })

  sigma <- exp(rowSums(exp.terms))

  return(-sum(dtobit(y = y, xTbeta = xTbeta, sigma = sigma, left = left, log = TRUE)))
}

#' Perform tobit regression assuming polynomial-based heteroscedasticity
#'
#' This function gets maximum-likelihood estimates of beta and gamma assuming the data follow a
#' tobit model with polynomial-based heteroscedasticity. NOTE: This is only for use in the EQ5D5L
#' project at SickKids. There is NO guarantee it will work for a general use-case.
#'
#' @importFrom stats model.matrix optim dnorm pnorm
#' @importFrom compiler cmpfun
#' @param formula a formula specifying the relationship between the tobit response and regressors
#' @param data a data.frame containing the variables specified in formula
#' @param start.beta an initial guess for beta
#' @param start.gamma an initial guess for gamma. This also implicity specifies the degree of the polynomial, i.e. length(gamma) - 1
#' @param left a number specifying where left-censoring occurred
#' @return a list containing the following parameter estimates (from maximum likelihood):\cr
#' \item{beta}{the regression coefficients}
#' \item{gamma}{the polynomial coefficients}
#' @export
tobitpolyvar <- function(formula, data, left = -1, start.beta = NULL, start.gamma = NULL)
{
  X <- model.matrix(formula, data = data)
  y <- data[, all.vars(formula)[1]]

  start.theta <- c(start.beta, start.gamma)

  J <- cmpfun(function(theta){objective(theta = theta, y = y, X = X, left = left)})

  optimum <- optim(par = start.theta, fn = J, method = "L")

  beta <- c(1, optimum$par[1:length(start.beta)])
  gamma <- optimum$par[(length(start.beta) + 1):length(optimum$par)]

  return(list(beta = beta, gamma = gamma))
}

# formula <- tto ~ 0 + mo2 + mo3 + mo4 + mo5 +
#   sc2 + sc3 + sc4 + sc5 +
#   ua2 + ua3 + ua45 +
#   pd2 + pd3 + pd4 + pd5 +
#   ad2 + ad3 + ad4 + ad5
#
# eqdata.tto.pv <- as.data.frame(model.matrix(tto ~ mo + sc + ua + pd + ad, eqdata.tto)[,-1])
# eqdata.tto.pv$tto <- eqdata.tto$tto
# data = eqdata.tto.pv
#
# tobitpolyvar(formula, data, left = -1,
#              start.beta = beta.tto.true[-1] + 0.6,
#              start.gamma = c(0.1, 0.1, 0.1, 0.1, 0.1))
