\name{plot_power_curve}
\alias{plot_power_curve}
\title{plot_power_curve}
\usage{
plot_power_curve(result, cexLegend = 1)
}
\arguments{
  \item{result}{the result of \code{\link{sample_size}} or
  \code{\link{est_power_curve}} function. The storeProcess
  parameter should be set as True when performing
  \code{\link{sample_size}} function. If you want to plot
  more than one curves in the same figure, the results from
  \code{\link{sample_size}} function should first be
  combined into a new list. At most five curves were
  allowed in one figure.}

  \item{cexLegend}{the cex for legend.}
}
\description{
A function to plot power curves based on the result of
\code{\link{sample_size}} or \code{\link{est_power_curve}}
function.
}
\examples{
#One power curve
result1<-sample_size(rho=2,phi0=1,lambda0=1,f=0.01,power=0.8,m=20000,m1=500,showMessage=T,storeProcess=T)
plot_power_curve(result1)
#More than one curves
result2<-sample_size(rho=4,phi0=1,lambda0=1,f=0.01,power=0.8,m=20000,m1=500,showMessage=T,storeProcess=T)
result<-list(result1,result2)
plot_power_curve(result)
}

