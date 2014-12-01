# TODO: Add comment
# 
# Author: zhaos
###############################################################################



##' sample_size
##' 
##' A function to estitamete the sample size for differential expression analysis of RNA-seq data.
##' 
##' A function to estitamete the sample size for differential expression analysis of RNA-seq data.
##' 
##' @param m Total number of genes for testing.
##' @param m1 Expected number of prognostic genes.
##' @param power Power to detecte prognostic genes.
##' @param f FDR level
##' @param w Ratio of normalization factors between two groups.
##' @param k Ratio of sample size between two groups.
##' @param rho minimum fold changes for prognostic genes between two groups.
##' @param lambda0 Average read counts for prognostic genes.
##' @param phi0 Dispersion for prognostic genes.
##' @param showMessage Logical. Display the message in the estimation process.
##' @param storeProcess Logical. Store the power and n in sample size or power estimation process.
##' @return Estimate sample size or a list including parameters and sample size in the process.
##' @export
##' @examples power<-0.8;rho<-2;lambda0<-5;phi0<-0.5;f<-0.01
##' sample_size(power=power, f=f,rho=rho, lambda0=lambda0, phi0=phi0)
sample_size<-function(power=0.8,m=20000, m1=200, f=0.1, k=1,w=1, rho=2, lambda0=5, phi0=1,showMessage=FALSE,storeProcess=FALSE){
	lambda0AppCut<-2000

	r1<-m1 * power
	beta<-1-power
	alpha_star<-r1*f/((m-m1)*(1-f))
	z_alpha<-qnorm(1-alpha_star/2, lower.tail=T)
	z_beta<-qnorm(power, lower.tail=TRUE)
	n_w<-( ( z_alpha + z_beta )^2* (1 + rho/w +phi0*lambda0*(1+rho^2)) )/ ( (rho-1)^2*lambda0 )
	
	start.point<-1
	end.point<-round(n_w)+30

	if (lambda0>=lambda0AppCut) {
		p1<-est_power_model(n=start.point,w=w,k=k, rho=rho, lambda0=lambda0, phi0=phi0, beta=beta, alpha=alpha_star)
	} else {
		p1<-est_power_root(n=start.point,w=w,k=k, rho=rho, lambda0=lambda0, phi0=phi0, beta=beta, alpha=alpha_star)
	}
	if (p1>0) {
		return(start.point)
	} else {
		step.power<-6
		if (lambda0>=lambda0AppCut) {
			n_Exact<-uniroot.integer(est_power_model, c(start.point, end.point), w=w,k=k, rho=rho, 
					lambda0=lambda0, phi0=phi0, beta=beta, alpha=alpha_star, pos.side=TRUE,
					step.up=TRUE, step.power=step.power,print.steps=showMessage)
		} else {
			n_Exact<-uniroot.integer(est_power_root, c(start.point, end.point), w=w,k=k, rho=rho, 
					lambda0=lambda0, phi0=phi0, beta=beta, alpha=alpha_star, pos.side=TRUE,
					step.up=TRUE, step.power=step.power,print.steps=showMessage) 
		}
		if (storeProcess) {
			n_Exact$process<-n_Exact$process[order(n_Exact$process[,1]),]
			n_Exact$process[,2]<-n_Exact$process[,2]+power
			n_Exact$parameters<-c(power,m, m1,  f, w, rho, lambda0, phi0,k)
			names(n_Exact$parameters)<-c("power","m", "m1",  "fdr", "w", "rho", "lambda0", "phi0","k")
		} else {
			n_Exact<-n_Exact$root
		}
		return(n_Exact)
	}
}

##' sample_size_distribution
##' 
##' A function to estitamete the sample size based on read counts and dispersion distribution in real data.
##' 
##' A function to estitamete the sample size based on read counts and dispersion distribution in real data.
##' 
##' @inheritParams sample_size
##' @inheritParams est_power_distribution
##' @return Estimate sample size or a list including parameters and sample size in the process.
##' @export
##' @examples \dontrun{
##' #Please note here the parameter repNumber was very small (5) to make the example code faster.
##' #We suggest repNumber should be at least set as 100 in real analysis.
##' sample_size_distribution(power=0.8,f=0.01,distributionObject="TCGA_READ",repNumber=5,showMessage=TRUE)
##' }
sample_size_distribution<-function(power=0.8,m=10000, m1=100, f=0.1, k=1,w=1, rho=2,showMessage=FALSE,storeProcess=FALSE,distributionObject,libSize,minAveCount=5,maxAveCount=2000,repNumber=100,dispersionDigits=1,seed=123,selectedGenes,pathway,species="hsa",countFilterInRawDistribution=TRUE,selectedGeneFilterByCount=FALSE){
	
	r1<-m1 * power
	beta<-1-power
	alpha_star<-r1*f/((m-m1)*(1-f))

	start.point<-1
	end.point<-800
	step.power<-5
	
	#process distribution
	temp<-selectDistribution(distributionObject=distributionObject,libSize=libSize,repNumber=repNumber,dispersionDigits=dispersionDigits,minAveCount=minAveCount,maxAveCount=maxAveCount,seed=seed,selectedGenes=selectedGenes,pathway=pathway,species=species,countFilterInRawDistribution=countFilterInRawDistribution,selectedGeneFilterByCount=selectedGeneFilterByCount)
	dispersionDistribution<-temp$selectedDispersion
	countDistribution<-temp$selectedCount	
	
	n_Exact<-uniroot.integer(est_power_distribution_sampleSize, c(start.point, end.point), w=w,k=k, rho=rho, 
				beta=beta, alpha=alpha_star, dispersionDistribution=dispersionDistribution,countDistribution=countDistribution,
				pos.side=TRUE,	step.up=TRUE, step.power=step.power,print.steps=showMessage)

	if (storeProcess) {
			n_Exact$process<-n_Exact$process[order(n_Exact$process[,1]),]
			n_Exact$process[,2]<-n_Exact$process[,2]+power
			n_Exact$parameters<-c(power,m, m1,  f, w, rho, k)
			names(n_Exact$parameters)<-c("power","m", "m1",  "fdr", "w", "rho","k")
	} else {
			n_Exact<-n_Exact$root
	}
	return(n_Exact)
}

est_power_distribution_sampleSize<-function(beta,...) {
		return(mean(est_power_distribution_root(...))-(1-beta))
}