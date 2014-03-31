\name{est_power}
\alias{est_power}
\title{est_power}
\usage{
est_power(n, w = 1, rho = 2, lambda0 = 5, phi0 = 1, alpha = 0.05, ...)
}
\arguments{
  \item{n}{Numer of samples.}

  \item{alpha}{alpha level.}

  \item{...}{other parameters for est_power_root function.}

  \item{w}{Ratio of normalization factors between two
  groups.}

  \item{rho}{minimum fold changes for prognostic genes
  between two groups.}

  \item{lambda0}{Average read counts for prognostic genes.}

  \item{phi0}{Dispersion for prognostic genes.}
}
\description{
A function to estitamete the power for differential
expression analysis of RNA-seq data.
}
\details{
A function to estitamete the power for differential
expression analysis of RNA-seq data.
}
\examples{
n<-20;valpha<-0.05;vw<-1.0;vrho<-1.5;vlambda0<-5;vphi_0<-0.5
est_power(n=n, alpha=valpha, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)
}

