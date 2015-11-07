\name{WAVK}
\alias{WAVK}
\title{WAVK statistic}
\description{
Computes statistic for testing the parametric form of a regression function, suggested by Wang, Akritas and Van Keilegom (2008).
}

\usage{
WAVK(z, kn = NULL)
}

\arguments{
  \item{z}{pre-filtered univariate time series (see formula (2.1) by Wang and Van Keilegom, 2007):
  \deqn{Z_i=\left(Y_{i+p}-\sum_{j=1}^p{\hat{\phi}_{j,n}Y_{i+p-j}} \right)-\left( f(\hat{\theta},t_{i+p})-\sum_{j=1}^p{\widehat{\phi}_{j,n}f(\widehat{\theta},t_{i+p-j})} \right)
,} 
where \eqn{Y_i} is observed time series of length \eqn{n}, \eqn{\widehat{\theta}} is an estimator of hypothesized parametric trend \eqn{f(\theta, t)}, and \eqn{\hat{\phi}_p=(\hat{\phi}_{1,n}, \ldots, \hat{\phi}_{p,n})'} are estimated coefficients of an autoregressive filter of order \eqn{p}. Missing values are not allowed.}

  \item{kn}{length of the local window.}
}


\value{
A list with following components:

  \item{Tn}{test statistic based on artificial ANOVA and defined by Wang and Van Keilegom (2007) as a difference of mean square for treatments (MST) and mean square for errors (MSE):
  \deqn{T_n= MST - MSE =\frac{k_{n}}{n-1} \sum_{t=1}^T \biggl(\overline{V}_{t.}-\overline{V}_{..}\biggr)^2 -
\frac{1}{n(k_{n}-1)} \sum_{t=1}^n \sum_{j=1}^{k_{n}}\biggl(V_{tj}-\overline{V}_{t.}\biggr)^2,}
where \eqn{\{V_{t1}, \ldots, V_{tk_n}\}=\{Z_j: j\in W_{t}\}}, \eqn{W_t} is a local window, \eqn{\overline{V}_{t.}} and \eqn{\overline{V}_{..}} are the mean of the \eqn{t}th group and the grand mean, respectively.} 
  \item{Tns}{standardized version of \code{Tn} according to Theorem 3.1 by Wang and Van Keilegom (2007): 
  \deqn{Tns = \left( \frac{n}{kn} \right)^{\frac{1}{2}}Tn \bigg/  \left(\frac{4}{3}\right)^{\frac{1}{2}} \sigma^2,}{Tns = Tn*(n/kn)^0.5 / (sigma^2 * (4/3)^0.5),}
  where \eqn{n} is length and \eqn{\sigma^2}{sigma^2} is variance of the time series. Robust difference-based 
  Rice's estimator (Rice, 1984) is used to estimate \eqn{\sigma^2}{sigma^2}.} 

  \item{p.value}{\eqn{p}-value for \code{Tns} based on its asymptotic \eqn{N(0,1)} distribution.}
  
}

\references{
Rice, J. (1984) Bandwidth choice for nonparametric regression. \emph{The Annals of Statistics} 12, 1215--1230.

Wang, L., Akritas, M. G. and Van Keilegom, I. (2008) An ANOVA-type nonparametric diagnostic test for heteroscedastic regression models. \emph{Journal of Nonparametric Statistics} 20(5), 365--382.
                                                  
Wang, L. and Van Keilegom, I. (2007) Nonparametric test for the form of parametric regression with time series errors. \emph{Statistica Sinica} 17, 369--386.
}

\author{Yulia R. Gel, Vyacheslav Lyubchich}

\seealso{\code{\link{wavk.test}}.
}

\examples{
z <- rnorm(300)
WAVK(z, kn=7)
}

\keyword{ts}
\keyword{trend}