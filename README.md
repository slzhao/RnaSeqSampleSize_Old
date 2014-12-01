RnaSeqSampleSize
============
* [Introduction](#Introduction)
* [User friendly web interface](#web)
* [Download and install](#download)
* [Example](#example)
* [Change log](#Change)

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

<a name="Change"/>
# Change log #
2014-11-18
RnaSeqSampleSize 0.99.2:
 * The vignette was switched to a BiocStyle;
 * Other improvement based on the comments from Bioconductor reviewer;

2014-10-19
RnaSeqSampleSize 0.99.1:
 * Parameter was used to determine the power of genes below minAveCount;
 * Some recommendations from BiocCheck were improved;

2014-10-16
RnaSeqSampleSize 0.99.0:
 * Submit to Bioconductor.

2014-10-15
RnaSeqSampleSize 0.5.0:
 * The datasets from TCGA were moved to RnaSeqSampleSizeData package;
 * The examples were improved;
 * The vignette was generated;

2014-10-10
RnaSeqSampleSize 0.4.0:
 * Rcpp package was used to make some functions faster;
 * The package structure was improved;
 * The web interface was improved;

2014-08-01
RnaSeqSampleSize 0.3.0:
 * Beta approximate method was improved and used in most of the cases;
 * Read count and dispersion distribution estimation function was improved;
 * Estimation power or sample size by read count and dispersion distribution was improved.

2014-06-01
RnaSeqSampleSize 0.2.0:
 * The function of estimation power or sample size by read count and dispersion distribution was imported and improved.

2014-03-30
RnaSeqSampleSize 0.1.7:
 * An approximate estimation method was used to decrease the running time when average read count (lambda0) was larger than 20;
 * A plotting power curve function plot_power_curve was provided;

2014-03-19
RnaSeqSampleSize 0.1.6:
 * The code was improved to fit the web interface;
 * The function est_power was provided for power estimation;
 * A bug was fixed;

2014-03-16
RnaSeqSampleSize 0.1.5:
 * The C code for R function dnbinom was improved to decrease the running time greatly;
 * The powers for different N in estimating sample size were returned to prepare power curve;

2014-03-01
RnaSeqSampleSize 0.1.4:
 * A user friendly web interface was provided at http://cqs.mc.vanderbilt.edu/shiny/RnaSeqSampleSize/;

2014-02-17
RnaSeqSampleSize 0.1.3:
 * The algorithm was improved to decrease the memory usage;

2014-02-04
RnaSeqSampleSize 0.1.2:
 * A more detail parameter to sample size table was provided;
 * The algorithm was improved to decrease the running time and memory usage;

2014-01-27
RnaSeqSampleSize 0.1.1:
 * The algorithm was improved to decrease the running time greatly;

2014-01-24
RnaSeqSampleSize 0.1.0:
 * A more detail parameter to sample size table was provided to make sample size estimation faster;
 * An estimated time to perform sample size estimation will be displayed;
 * A bug fixed when sample size N equals to 1;

2014-01-10
RnaSeqSampleSize 0.0.1:
 * First version.

