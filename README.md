RnaSeqSampleSize
============
* [Introduction](#Introduction)
* [Change log](#Change)
* [Download and install](#download)
* [Example](#example)
* [Usage](#usage)

<a name="Introduction"/>
# Introduction #
Sample size calculation is an important issue in the experimental design of biomedical research. For RNA-seq experiments, the sample size calculation method based on the Poisson model has been proposed; however, when there are biological replicates, RNA-seq data could exhibit variation significantly greater than the mean (i.e. over-dispersion). The Poisson model cannot appropriately model the over-dispersion, and in such cases, the negative
binomial model has been used as a natural extension of the Poisson model. Because the field currently lacks a sample size calculation method based on the negative binomial model for assessing differential expression analysis of RNA-seq data, we propose a method to calculate the sample size.

RnaSeqSampleSize package is based on the paper **Sample size calculation based on exact test for assessing differential expression analysis in RNA-seq data**, which was published in BMC Bioinformatics. The paper can be accessed in [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/?term=24314022).

<a name="Change"/>
# Change log #

2014-02-04
RnaSeqSampleSize 1.1.2:
 * (1) A more detail parameter to sample size table was provided;
 * (2) The algorithm was improved to decrease the running time and memory usage;

2014-01-27
RnaSeqSampleSize 1.1.1:
 * The algorithm was improved to decrease the running time greatly;

2014-01-24
RnaSeqSampleSize 1.1.0:
 * A more detail parameter to sample size table was provided to make sample size estimation faster;
 * An estimated time to perform sample size estimation will be displayed;
 * A bug fixed when sample size N equals to 1;

2014-01-10
RnaSeqSampleSize 1.0.0:
 * First version;

<a name="download"/>
# Download and install #
**You can directly download the zip file (Windows system) or tar.gz file (Linux system) of RnaSeqSampleSize package** from [github release page](https://github.com/slzhao/RnaSeqSampleSize/releases).

Or you can download the source codes of RnaSeqSampleSize package from [github](https://github.com/slzhao/RnaSeqSampleSize) by the following commands (If git has already been installed in your computer). And then you need to built the source codes into a package.

	#This is not necessary if you have already downloaded the zip or tar.gz file from release page
	#The source codes of RnaSeqSampleSize package will be downloaded to your current directory
	git clone https://github.com/slzhao/RnaSeqSampleSize.git
	#Built the source codes into a package
	R CMD build RnaSeqSampleSize
	#The file of RnaSeqSampleSize package will be generated in current folder

RnaSeqSampleSize package requires ssanv packages. You can enter R and use following R codes to install it. 

	#Install required CRAN packages
	install.packages("ssanv")

Then you can enter R and use following R codes to install RnaSeqSampleSize in Windows or Linux.
	
	#Install RnaSeqSampleSize in Windows system, assume RnaSeqSampleSize_1.1.2.zip file was in current directory
	install.packages("RnaSeqSampleSize_1.1.2.zip")
	#Install RnaSeqSampleSize in Linux system, assume RnaSeqSampleSize_1.1.2.tar.gz file was in current directory
	install.packages("RnaSeqSampleSize_1.1.2.tar.gz")

<a name="example"/>
# Example #
After you have installed RnaSeqSampleSize package. You can enter R and use following R codes to see the examples for it.
	
	#Load package
	library("RnaSeqSampleSize")
	#View help files
	?sample_size
	#A simple example
	example(sample_size)

<a name="usage"/>
# Usage #

In RnaSeqSampleSize package, sample_size is the primary function to estimate sample size based on different parameters.
	
	#List of all possible parameters for sample_size function:
	# m: Total number of genes for testing.
	# m1: Expected number of prognostic genes.
	# power: Power to detecte prognostic genes.
	# f: FDR level
	# w: Ratio of normalization factors between two groups.
	# rho: minimum fold changes for prognostic genes between two groups.
	# lambda0: Average read counts for prognostic genes.
	# phi_0: Dispersion for prognostic genes.

	#Set the parameters:
	vm<-10000;vm1<-100;vpower<-0.8;vf<-0.05;vw<-1.0;vrho<-2.0;vlambda0<-5;vphi_0<-0.5
	#Estitamete sample size
	sample_size(m=vm, m1=vm1, power=vpower, f=vf, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)