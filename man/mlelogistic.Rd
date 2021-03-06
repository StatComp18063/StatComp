% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlelogistic.R
\name{mlelogistic}
\alias{mlelogistic}
\title{Computing the MLE for logistic distribution}
\usage{
mlelogistic(x, theta0 = mean(x), numstep = 100, eps = 1e-04)
}
\arguments{
\item{x}{vector, the sample}

\item{theta0}{the mean of x}

\item{numstep}{maximum repeat times}

\item{eps}{The iteration error}
}
\value{
theta1 is the MLE

check is the iteration error

realnumsteps is the actual number of iterations
}
\description{
Computing the MLE for logistic distribution using R,
location parameter is theta, scale parameter is 1.
}
\examples{
x=rlogis(100,2,1)
mlelogistic(x)
}
