\name{sample_size}
\alias{sample_size}
\title{sample_size}
\usage{
sample_size(m = 10000, m1 = 100, power = 0.8, f = 0.1, w = 1,
  rho = 2, lambda0 = 5, phi0 = 1, nTable = system.file("extdata",
  "nTable.csv", package = "RnaSeqSampleSize"), showMessage = F,
  storeProcess = F)
}
\arguments{
  \item{m}{Total number of genes for testing.}

  \item{m1}{Expected number of prognostic genes.}

  \item{power}{Power to detecte prognostic genes.}

  \item{f}{FDR level}

  \item{w}{Ratio of normalization factors between two
  groups.}

  \item{rho}{minimum fold changes for prognostic genes
  between two groups.}

  \item{lambda0}{Average read counts for prognostic genes.}

  \item{phi0}{Dispersion for prognostic genes.}

  \item{nTable}{A table included the estitameted sample
  numbers (N) with different parameters.}

  \item{showMessage}{Logical. Display the message in the
  estimation process.}

  \item{storeProcess}{Logical. Store the power and n in
  sample size estimation process.}
}
\description{
A function to estitamete the sample size for differential
expression analysis of RNA-seq data.
}
\details{
A function to estitamete the sample size for differential
expression analysis of RNA-seq data.
}
\examples{
#Input initial parameters. These parameters were stored in our N table so you will get results in one second.
vm<-10000;vm1<-100;vpower<-0.8;vf<-0.05;vw<-1.0;vrho<-2.0;vlambda0<-5;vphi_0<-2
sample_size(m=vm, m1=vm1, power=vpower, f=vf, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)
#Estitamete sample size for parameters not stored, may use one minutes.
vm<-30000;vm1<-1000;vpower<-0.8;vf<-0.05;vw<-1.0;vrho<-1.5;vlambda0<-5;vphi_0<-0.5
sample_size(m=vm, m1=vm1, power=vpower, f=vf, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)
}

