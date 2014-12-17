RnaSeqSampleSize
============
* [Introduction](#Introduction)
* [User friendly web interface](#web)
* [Download and install](#download)
* [Example](#example)

<a name="Introduction"/>
# Introduction #
Sample size estimation is the most important issue in the design of RNA sequencing experiments. However, thousands of genes are quantified and tested for differential expression simultaneously in RNA-seq experiments. The false discovery rate for statistic tests should be controlled. At the same time, the thousands of genes have widely distributed read counts and dispersions, which were often estimated by experience or set at the most conservative values in previous sample size estimation methods. As a result, the estimated sample size will be inaccurate or over-estimated.

To solve these issues, we developed a sample size estimation method based on the distributions of gene read counts and dispersions from real data. Datasets from the user's preliminary experiments or the Cancer Genome Atlas (TCGA) can be used as reference. The read counts and their related dispersions will be selected randomly from the reference based on their distributions, and from that, the power and sample size will be estimated and summarized.

<a name="web"/>
# User friendly web interface #
A user friendly web interface for RnaSeqSampleSize package is provided at [CQS website](http://cqs.mc.vanderbilt.edu/shiny/RnaSeqSampleSize/). Most of the features in Example section can be found in this website.


<a name="download"/>
# Download and install #
You can download and install RnaSeqSampleSize package from [github](https://github.com/slzhao/RnaSeqSampleSize) by the following commands.

	#Running the following codes in your R
	library(devtools)
    install_github("slzhao/RnaSeqSampleSizeData")
    install_github("slzhao/RnaSeqSampleSize")

<a name="example"/>
# Example #
After you have installed RnaSeqSampleSize package. You can enter R and use following R codes to view documents and perform examples for it.
	
	#Load package
	library("RnaSeqSampleSize")
	#View vignette
	browseVignettes(package="RnaSeqSampleSize")
	#View help files
	?sample_size
	#Examples for sample size or power estimation by single read count and dispersion
	sample_size(power=0.8, f=0.01,rho=2, lambda0=5, phi0=0.5)
	est_power(n=63, rho=2, lambda0=5, phi0=0.5,f=0.01)
	#Examples for power estimation by prior real data, may use 3 minutes
	est_power_distribution(n=65,f=0.01,rho=2,distributionObject="TCGA_READ",repNumber=5)
	sample_size_distribution(power=0.918,f=0.01,rho=2,distributionObject="TCGA_READ",repNumber=5)
	#Examples for power curve generation
	result1<-est_power_curve(n=63, f=0.01, rho=2, lambda0=5, phi0=0.5)
	result2<-est_power_curve(n=63, f=0.05, rho=2, lambda0=5, phi0=0.5)
	plot_power_curve(list(result1,result2))

