\name{est_power_curve}
\alias{est_power_curve}
\title{est_power_curve}
\usage{
est_power_curve(n, w = 1, rho = 2, lambda0 = 5, phi0 = 1,
  alpha = 0.05, ...)
}
\arguments{
  \item{n}{Numer of samples.}

  \item{alpha}{alpha level.}

  \item{...}{other parameters for est_power function.}

  \item{w}{Ratio of normalization factors between two
  groups.}

  \item{rho}{minimum fold changes for prognostic genes
  between two groups.}

  \item{lambda0}{Average read counts for prognostic genes.}

  \item{phi0}{Dispersion for prognostic genes.}
}
\description{
A function to estitamete the power curve for differential
expression analysis of RNA-seq data.
}
\examples{
result<-est_power_curve(n=20, alpha=0.05, w=1, rho=1.5, lambda0=5, phi0=0.5)
plot_power_curve(result)
}

