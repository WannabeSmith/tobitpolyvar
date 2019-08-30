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

llTobit <- function(y, xTbeta, sigma, left)
{
  return(sum(dtobit(y = y, xTbeta = xTbeta, sigma = sigma, left = left, log = TRUE)))
}

objective.tobitpolyvar <- function(params, y, X, left)
{
  beta <- c(1, params[1:ncol(X)])
  xTbeta <- cbind(1, X) %*% beta
  gamma <- params[(ncol(X)+1):length(params)]

  sigma <- exp(X %*% gamma)

  return(llTobit(y = y, xTbeta = xTbeta, sigma = sigma, left = left))
}

#' Perform tobit regression assuming polynomial-based heteroscedasticity
#'
#' This function gets maximum-likelihood estimates of beta and gamma assuming the data follow a
#' tobit model with polynomial-based heteroscedasticity. NOTE: This is only for use in the EQ5D5L
#' project at SickKids. There is NO guarantee it will work for a general use-case.
#'
#' @importFrom stats model.matrix optim dnorm pnorm
#' @importFrom compiler cmpfun
#' @importFrom survival survreg Surv
#' @param formula a formula specifying the relationship between the tobit response and regressors
#' @param data a data.frame containing the variables specified in formula
#' @param start.beta an initial guess for beta
#' @param start.gamma an initial guess for gamma. This also implicity specifies the degree of the polynomial, i.e. length(gamma) - 1
#' @param left a number specifying where left-censoring occurred
#' @return a list containing the following parameter estimates (from maximum likelihood):\cr
#' \item{beta}{the regression coefficients}
#' \item{gamma}{the polynomial coefficients}
#' @export
tobitpolyvar <- function(formula, data, left = -1, start.beta, start.gamma)
{
  X <- model.matrix(formula, data = data)
  y <- data[, all.vars(formula)[1]]

  start.params <- c(start.beta, "gamma" = start.gamma)

  J <- cmpfun(function(params){-objective.tobitpolyvar(params = params, y = y, X = X, left = left)})

  optimum <- optim(par = start.params, fn = J, method = "BFGS")

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

# log-likelihood of logistic regression model
llLogisticReg <- function(y, xTbeta)
{
  sum(-(1 - y) * xTbeta - log(1 + exp(-xTbeta)))
}

# log-likelihood of the hybrid tobit-logit regression model
llHyregPolyVar <- function(y.tobit, y.discrete, sigma,
                            xTbeta.tobit, xTbeta.discrete, left)
{
  llTobit(y = y.tobit, xTbeta = xTbeta.tobit, sigma = sigma, left = left) +
    llLogisticReg(y = y.discrete, xTbeta = xTbeta.discrete)
}

objective.hyregpolyvar <- function(params, y.tobit, y.discrete,
                                    X.tobit, X.discrete, left)
{
  beta.tobit <- c(1, params[1:ncol(X.tobit)])
  gamma <- params[(ncol(X.tobit)+1):(length(params)-1)]
  theta <- params[length(params)]
  beta.discrete <- c(0, theta * beta.tobit[-1])


  xTbeta.tobit <- cbind(1, X.tobit) %*% beta.tobit
  xTbeta.discrete <- cbind(1, X.discrete) %*% (beta.discrete)

  sigma <- exp(X.tobit %*% gamma)

  return(llHyregPolyVar(y.tobit = y.tobit, y.discrete = y.discrete, sigma = sigma,
                         left = left, xTbeta.tobit = xTbeta.tobit,
                         xTbeta.discrete = xTbeta.discrete))
}

#' Perform hybrid tobit-logit regression assuming polynomial-based heteroscedasticity
#'
#' This function gets maximum-likelihood estimates of beta, gamma, and theta assuming
#' the data follow a hybrid tobit-logit model with polynomial-based heteroscedasticity.
#' NOTE: This is only for use in the EQ5D5L project at SickKids.
#' There is NO guarantee it will work for a general use-case.
#'
#' @importFrom stats model.matrix optim dnorm pnorm
#' @importFrom compiler cmpfun
#' @importFrom survival survreg Surv
#' @param formula.tobit a regression formula describing the relationship between the tobit response and the covariates
#' @param formula.discrete a regression formula describing the relationship between the bernoulli response and the covariates
#' @param data.tobit the data.frame containing the tobit responses and covariates
#' @param data.discrete the data.frame containing the bernoulli responses and covariates
#' @param left a number specifying where left-censoring occurred
#' @param start.beta an initial guess for beta
#' @param start.gamma an initial guess for gamma. This also implicity specifies the degree of the polynomial, i.e. length(gamma) - 1
#' @param start.theta an initial guess for theta, the multiplicative factor between betas
#' @return a list containing the following parameter estimates (from maximum likelihood):\cr
#' \item{beta}{the regression coefficients}
#' \item{gamma}{the polynomial coefficients}
#' \item{theta}{the scalar multiplicative factor between betas}
#' @export
hyregpolyvar <- function(formula.tobit, formula.discrete, data.tobit,
                         data.discrete, left = -1, start.beta,
                         start.gamma, start.theta)
{
  X.tobit <- model.matrix(formula.tobit, data = data.tobit)
  y.tobit <- data.tobit[, all.vars(formula.tobit)[1]]

  X.discrete <- model.matrix(formula.discrete, data = data.discrete)
  y.discrete <- data.discrete[, all.vars(formula.discrete)[1]]

  start.params <- c(start.beta, "gamma" = start.gamma, "theta" = start.theta)

  J <- cmpfun(function(params){-objective.hyregpolyvar(params = params, y.tobit = y.tobit,
                                                        y.discrete = y.discrete,
                                                        X.tobit = X.tobit, X.discrete = X.discrete,
                                                        left = left)})

  optimum <- optim(par = start.params, fn = J, method = "BFGS")

  beta <- c(1, optimum$par[1:length(start.beta)])
  gamma <- optimum$par[(length(start.beta) + 1):(length(optimum$par) - 1)]
  theta <- optimum$par[length(optimum$par)]

  return(list(beta = beta,
              gamma = gamma,
              theta = theta))
}

# formula.tobit <- tto ~ 0 + mo2 + mo3 + mo4 + mo5 +
#   sc2 + sc3 + sc4 + sc5 +
#   ua2 + ua3 + ua45 +
#   pd2 + pd3 + pd4 + pd5 +
#   ad2 + ad3 + ad4 + ad5
#
# formula.discrete <- value ~ 0 + mo2 + mo3 + mo4 + mo5 +
#   sc2 + sc3 + sc4 + sc5 +
#   ua2 + ua3 + ua45 +
#   pd2 + pd3 + pd4 + pd5 +
#   ad2 + ad3 + ad4 + ad5
#
# eqdata.tto.pv <- as.data.frame(model.matrix(tto ~ mo + sc + ua + pd + ad, eqdata.tto)[,-1])
# eqdata.tto.pv$tto <- eqdata.tto$tto
# data.tobit = eqdata.tto.pv
#
# data.discrete <- eqdata.dce
#
# start.beta <- beta.tto.true[-1]
# start.gamma <- sim.spec$variance$tto$sigma$gamma
# start.theta <- mean(beta.dce.true[-1] / beta.tto.true[-1])
#
# hybrid.pv.ests <- hyregpolyvar(formula.tobit = formula.tobit, formula.discrete = formula.discrete,
#                                 data.tobit = data.tobit, data.discrete = data.discrete,
#                                 left = -1, start.beta = start.beta, start.gamma = start.gamma,
#                                 start.theta = start.theta)




