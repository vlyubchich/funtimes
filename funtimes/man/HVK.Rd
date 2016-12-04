\name{HVK}
\alias{HVK}
\title{HVK Estimator}
\description{
Estimates coefficients in non-parametric autoregression using the difference-based approach by Hall and Van Keilegom (2003).
}

\usage{
HVK(X, m1=NULL, m2=NULL, ar.order=1)
}

\arguments{
  \item{X}{univariate time series. Missing values are not allowed.}
  \item{m1, m2}{subsidiary smoothing parameters. Default \code{m1 = round(length(X)^(0.1))}, \code{m2 = round(length(X)^(0.5))}.}
  \item{ar.order}{order of the non-parametric autoregression (specified by user).}
}

\details{
First, autocovariances are estimated (formula (2.6) by Hall and Van Keilegom, 2003):
\deqn{\hat{\gamma}(0)=\frac{1}{m_2-m_1+1}\sum_{m=m_1}^{m_2}\frac{1}{2(n-m)}\sum_{i=m+1}^{n}\{(D_mX)_i\}^2,}
\deqn{\hat{\gamma}(j)=\hat{\gamma}(0)-\frac{1}{2(n-j)}\sum_{i=j+1}^n\{(D_jX)_i\}^2,}
where \eqn{n}=\code{length(X)} is sample size, \eqn{D_j} is a difference operator such that \eqn{(D_jX)_i=X_i-X_{i-j}}.
Then, Yule-Walker method is used to derive autoregression coefficients. 
}

\value{
Vector of length \code{ar.order} with estimated autoregression coefficients.
}

\references{
Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for inference in nonparametric regression with time series errors. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 65: 443--456. DOI: \href{http://dx.doi.org/10.1111/1467-9868.00395}{10.1111/1467-9868.00395}
}

\author{Yulia R. Gel, Vyacheslav Lyubchich, Xingyu Wang}

\examples{
X <- arima.sim(n=300, list(order=c(1,0,0), ar=c(0.6)))
HVK(as.vector(X), ar.order=1)
}

\keyword{ts}

