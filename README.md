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
A user friendly web interface for RnaSeqSampleSize package is provided at \url{http://cqs.mc.vanderbilt.edu/shiny/RnaSeqSampleSize/}. Most of the functions in Examples section can be performed in this website.


<a name="download"/>
# Download and install #
You can download and install RnaSeqSampleSize package from [github](https://github.com/slzhao/RnaSeqSampleSize) by the following commands.

	#Running the following codes in your R
	library(devtools)
    install_github("slzhao/RnaSeqSampleSizeData")
    install_github("slzhao/RnaSeqSampleSize")

<a name="example"/>
# Example #
After you have installed RnaSeqSampleSize package. You can enter R and use following R codes to see the examples for it.
	
	#Load package
	library("RnaSeqSampleSize")
	#View vignette
	browseVignettes(package="RnaSeqSampleSize")
	#View help files
	?sample_size
	#Examples for sample size estimation by single read count and dispersion
	example(sample_size)
	#Examples for power estimation by prior real data
	example(est_power_distribution)
