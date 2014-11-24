##' plot_power_curve
##' 
##' A function to plot power curves based on the result of \code{\link{sample_size}} or \code{\link{est_power_curve}} function.
##' 
##'
##' 
##' @param result the result of \code{\link{sample_size}} or \code{\link{est_power_curve}} function. The storeProcess parameter should be set as True when performing \code{\link{sample_size}} function. If you want to plot more than one curves in the same figure, the results from \code{\link{sample_size}} function should first be combined into a new list. At most five curves were allowed in one figure.
##' @param cexLegend the cex for legend.
##' @param pch Either an integer specifying a symbol or a single character to be used as the default in plotting points.
##' @param lwd The line width.
##' @param las Numeric in {0,1,2,3}; the style of axis labels.
##' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
##' @param col The line color.
##' @inheritParams graphics::plot.default
##' @export
##' @examples result1<-sample_size(rho=2,phi0=1,lambda0=1,f=0.01,power=0.8,m=20000,m1=500,showMessage=TRUE,storeProcess=TRUE)
##' result2<-sample_size(rho=4,phi0=1,lambda0=1,f=0.01,power=0.8,m=20000,m1=500,showMessage=TRUE,storeProcess=TRUE)
##' plot_power_curve(list(result1,result2))
plot_power_curve<-function(result,cexLegend=1,type="b",xlab="Sample Size",ylab="Power",pch=16,lwd=3,las=1,cex=1.5,main="Power Curve",col="red") {
	if (identical(names(result),c("iter","f.root","root","process","parameters")) || identical(names(result),c("process","parameters","power"))) { #Only one result
		plot(result$process[,1],result$process[,2],type=type,xlab=xlab,ylab=ylab,pch=pch,lwd=lwd,las=las,cex=cex,main=main,col=col)
		legend("bottomright",legend=paste(paste(names(result$parameters),result$parameters,sep="="),collapse=";"),bty="n",text.col="red",cex=0.9)
		abline(h= result$parameters["power"],lty=2,col="grey")
	} else { #more than one result
		if (length(result)>5) { #at most 5 curves were allowed
			result<-result[(length(result)-4):length(result)]
			warning("At most 5 curves were allowed in plot_power_curve function, the last 5 in result parameter will be used")
		}
		resultRange<-apply(sapply(result,function(x) (x$process)[nrow(x$process),]),1,max)
		plot(c(0,resultRange[1]),c(0,resultRange[2]),type="n",xlab=xlab,ylab=ylab,las=las,main=main)
		col<-c("brown1","steelblue2","mediumpurple2","seagreen3","lightgoldenrod")
		legendEach<-""
		for (x in 1:length(result)) {
			resultEach<-result[[x]]
			lines(resultEach$process[,1],resultEach$process[,2],type="b",pch=16,lwd=3,cex=1.5,col=col[x])
			legendEach[x]<-paste(paste(names(resultEach$parameters),resultEach$parameters,sep="="),collapse=";")
		}
		abline(h= result[[length(result)]]$parameters["power"],lty=2,col="grey")
		legend("bottomright",legend=legendEach,bty="n",text.col=col[1:length(result)],cex=cexLegend)
	}
}

##' est_power_curve
##' 
##' A function to estitamete the power curve for differential expression analysis of RNA-seq data.
##' 
##'
##'
##' @param ... other parameters for est_power function. 
##' @inheritParams est_power
##' @inheritParams sample_size
##' @return A list including parameters, sample size and power. 
##' @export
##' @examples \dontrun{
##' result1<-est_power_curve(n=63, f=0.01, rho=2, lambda0=5, phi0=0.5)
##' result2<-est_power_curve(n=63, f=0.05, rho=2, lambda0=5, phi0=0.5)
##' plot_power_curve(list(result1,result2))
##' }
est_power_curve<-function(n, w=1, rho=2, lambda0=5, phi0=1,alpha=0.05,f=0.05,...) {
	if (n<=10) {
		sampleSizeList<-1:n
	} else if (! n %% 5) {
		sampleSizeList<-c(1,seq(5,n,by=5))
	} else {
		sampleSizeList<-c(1,seq(5,n,by=5),n)
	}
	powerList<-NULL
	for (i in 1:length(sampleSizeList)) {
		if (!missing(f)) {
			powerList[i]<-est_power(n=sampleSizeList[i], w=w, rho=rho, lambda0=lambda0, phi0=phi0, f=f,...)
		} else {
			powerList[i]<-est_power(n=sampleSizeList[i], w=w, rho=rho, lambda0=lambda0, phi0=phi0, alpha=alpha,...)
		}
	}
	process<-cbind(N=sampleSizeList,Power=powerList)
	if (!missing(f)) {
		parameters<-c(n, f, w, rho, lambda0, phi0)
		names(parameters)<-c("n","fdr", "w", "rho", "lambda0", "phi0")
	} else {
		parameters<-c(n, alpha, w, rho, lambda0, phi0)
		names(parameters)<-c("n","alpha", "w", "rho", "lambda0", "phi0")
	}
	return(list(process=process,parameters=parameters,power=process[nrow(process),2]))
}

