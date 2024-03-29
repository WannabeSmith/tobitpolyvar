% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tobitpolyvar.R
\name{hyregpolyvar}
\alias{hyregpolyvar}
\title{Perform hybrid tobit-logit regression assuming polynomial-based heteroscedasticity}
\usage{
hyregpolyvar(formula.tobit, formula.discrete, data.tobit, data.discrete,
  left = -1, start.beta, start.gamma, start.theta)
}
\arguments{
\item{formula.tobit}{a regression formula describing the relationship between the tobit response and the covariates}

\item{formula.discrete}{a regression formula describing the relationship between the bernoulli response and the covariates}

\item{data.tobit}{the data.frame containing the tobit responses and covariates}

\item{data.discrete}{the data.frame containing the bernoulli responses and covariates}

\item{left}{a number specifying where left-censoring occurred}

\item{start.beta}{an initial guess for beta}

\item{start.gamma}{an initial guess for gamma. This also implicity specifies the degree of the polynomial, i.e. length(gamma) - 1}

\item{start.theta}{an initial guess for theta, the multiplicative factor between betas}
}
\value{
a list containing the following parameter estimates (from maximum likelihood):\cr
\item{beta}{the regression coefficients}
\item{gamma}{the polynomial coefficients}
\item{theta}{the scalar multiplicative factor between betas}
}
\description{
This function gets maximum-likelihood estimates of beta, gamma, and theta assuming
the data follow a hybrid tobit-logit model with polynomial-based heteroscedasticity.
NOTE: This is only for use in the EQ5D5L project at SickKids.
There is NO guarantee it will work for a general use-case.
}
