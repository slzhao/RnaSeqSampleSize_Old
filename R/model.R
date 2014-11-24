# TODO: Add comment
# 
# Author: zhaos
###############################################################################

est_power_model<-function(n, w=1,k=1, rho=2, lambda0=5, phi0=1, beta=0.2, alpha=0.05, error=0.001){
	model<-generateModel(lambda0=c(1,5,10),n=n, w=w,k=k, rho=rho, phi0=phi0,showMessage=FALSE,alpha=alpha)
	a<-modelPower(model,n=n,w=w,lambda0=lambda0,rho=rho, phi0=phi0)
	return(a-(1-beta))
}

modelPower<-function(model, n=10,w=1, rho=2, lambda0=5, phi0=1, error=0.001) {
	mu0<-lambda0
	mu1<-mu0*(rho*w)
	phi1<-phi0
	
	modelX1<-model[1]
	modelIntercept1<-model[2]
	modelX2<-model[3]
	modelIntercept2<-model[4]
	a<-0
	q0_u<-qnbinom(1-error, size=n/phi0, mu=n*mu0)
	q0_l<-qnbinom(error, size=n/phi0, mu=n*mu0)
	q1_u<-qnbinom(1-error, size=n/phi1, mu=n*mu1)
	q1_l<-qnbinom(error, size=n/phi1, mu=n*mu1)
#	dnbinomQ1<-dnbinom(q1_l:q1_u, mu=(n*mu1), size=n/phi_1)
	pnbinomQ1<-pnbinom((q1_l-1):q1_u, mu=(n*mu1), size=n/phi1)
	dnbinomQ0<-dnbinom(q0_l:q0_u, mu=(n*mu0), size=n/phi0)
#	pnbinomQ0<-pnbinom((q0_l-1):q0_u, mu=(n*mu0), size=n/phi_0)
	aNRow<-q1_u-q1_l+1
	aNCol<-q0_u-q0_l+1
	
	if (!is.na(modelX1)) {
		y<-round(modelX1*(q0_l:q0_u)+modelIntercept1)
		y[y<q1_l]<-q1_l
		y[y>q1_u]<-aNRow+q1_l-1
		i<-1:aNCol
		a<-a+sum((pnbinomQ1[aNRow]-pnbinomQ1[(y[i]-q1_l+1)])*dnbinomQ0[i])
	}
	
	if (!is.na(modelX2)) {
		y<-round(modelX2*(q0_l:q0_u)+modelIntercept2)
		y[y>q1_u]<-q1_u
		y[y<q1_l]<-aNRow+q1_l-1
		i<-1:aNCol
		a<-a+sum((pnbinomQ1[aNRow]-pnbinomQ1[(y[i]-q1_l+1)])*dnbinomQ0[i])
	}
	return(a)
}

generateModel<-function(lambda0=c(1,5,10,20),n, w=1,k=1, rho=2.0, phi0=1, alpha=0.05, error=0.001,showMessage=FALSE) {
	result<-NULL
	X1<-NULL
	X2<-NULL
	Y1<-NULL
	Y2<-NULL
	for (i in 1:length(lambda0)) {		
		temp<-est_power_root(n=n, w=w,k=k, rho=rho, lambda0=lambda0[i], phi0=phi0, alpha=alpha, error=error,returnDetail=TRUE)
		result[[i]]<-temp[["matrix"]]
		X1<-c(X1,temp[["X1"]])
		X2<-c(X2,temp[["X2"]])
		Y1<-c(Y1,temp[["Y1"]])
		Y2<-c(Y2,temp[["Y2"]])
	}
	
	if (showMessage) {
		maxRange<-max(c(X1,X2,Y1,Y2),na.rm=TRUE)
		maxRange<-maxRange+maxRange/5
		plot(c(0,maxRange),c(0,maxRange),type="l",xlab="q0",ylab="q1",main=paste("n=",n,"; rho=",rho,"; phi=",phi0,sep=""),lty=2)
		
		for (i in 1:length(result)) {
			temp1<-range(as.integer(colnames(result[[i]])))[c(1,2,2,1)]
			temp2<-range(as.integer(row.names(result[[i]])))[c(2,2,1,1)]
			temp3<-which(result[[i]]<=alpha,arr.ind=TRUE)
			
			temp3[,1]<-as.integer(temp3[,2])+as.integer(colnames(result[[i]])[1])
			temp3[,2]<-as.integer(row.names(temp3))
			points((temp3),pch=16,cex=0.2,col="grey")
			polygon(temp1,temp2,border=rainbow(length(result))[i],lwd=2)
		}
		legend("bottomright",legend=c(paste("Counts=",lambda0,sep="")),lwd=2,col=rainbow(length(result)),bty="n")
	}
	
	if (length(X1)>3) {
		model<-lm(Y1~X1)
		modelIntercept1<-model$coefficients[1]
		modelX1<-model$coefficients[2]
		if (showMessage) {
			cat("Model1 summary:\n")
			print(summary(model))
			abline(modelIntercept1,modelX1,col="red")
		}
	} else {
		modelX1<-NA
		modelIntercept1<-NA
	}
	
	if (length(X2)>3) {
		model<-lm(Y2~X2)
		modelIntercept2<-model$coefficients[1]
		modelX2<-model$coefficients[2]
		if (showMessage) {
			cat("Model2 summary:\n")
			print(summary(model))
			abline(modelIntercept2,modelX2,col="green")
		}
	} else {
		modelX2<-NA
		modelIntercept2<-NA
	}
	
	return(c(modelX1=as.numeric(modelX1),modelIntercept1=as.numeric(modelIntercept1),modelX2=as.numeric(modelX2),modelIntercept2=as.numeric(modelIntercept2)))
}

