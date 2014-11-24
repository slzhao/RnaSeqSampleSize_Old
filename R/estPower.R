##' est_power
##' 
##' A function to estitamete the power for differential expression analysis of RNA-seq data.
##' 
##'
##' 
##' @param n Numer of samples.
##' @param alpha alpha level.
##' @return Estimate power
##' @inheritParams sample_size
##' @export
##' @examples n<-63;rho<-2;lambda0<-5;phi0<-0.5;f<-0.01
##' est_power(n=n, rho=rho, lambda0=lambda0, phi0=phi0,f=f)
est_power<-function(n, w=1, rho=2, lambda0=5, phi0=1,alpha=0.05,f,m=20000,m1=200){
	if (!missing(f)) {#FDR power
		power_fdr_100<-est_power_root_fdr(power=1,n=n, w=w, rho=rho, lambda0=lambda0, phi0=phi0,m=m,m1=m1,fdr=f)
		if (power_fdr_100<0.01) {
			return(0)
		} else {
			power_fdr<-uniroot.integer(f=est_power_root_fdr,interval=c(1,100),n=n, w=w, rho=rho, lambda0=lambda0, phi0=phi0,m=m,m1=m1,fdr=f)
			return((power_fdr$root-1)/100)
		}
	} else {#alpha power
		power<-est_power_root(n=n, w=w, rho=rho, lambda0=lambda0, phi0=phi0, alpha=alpha)
		return(power+0.8)
	}
}

est_power_root_fdr<-function(power,n, w, rho, lambda0, phi0,fdr,m,m1,...) {
	alpha_star<-m1 * power/100*fdr/((m-m1)*(1-fdr))
#	cat(paste0(alpha_star,"\n"))
	beta<-1-power/100
	est_power_root(n=n, w=w, rho=rho, lambda0=lambda0, phi0=phi0, alpha=alpha_star,beta=beta,...)
}

est_power_root<-function(n,k=1, w=1, rho=2.0, lambda0=5, phi0=1, beta=0.2, alpha=0.05, bigCount=900,error=0.001,returnDetail=FALSE){
	mu0<-lambda0
	mu1<-mu0*(rho*w)
	phi1<-phi0
	q0_u<-qnbinom(1-error, size=n/phi0, mu=n*mu0)
	q0_l<-qnbinom(error, size=n/phi0, mu=n*mu0)
	q1_u<-qnbinom(1-error, size=k*n/phi1, mu=k*n*mu1)
	q1_l<-qnbinom(error, size=k*n/phi1, mu=k*n*mu1)
	
	if (min(q0_u,q0_l,q1_u,q1_l)>=bigCount) { #beta approx, I am not very sure if we need to *n here to determine if bigCount
		method<-"beta"
	} else { #NB
		method<-"nb"
	}
#	method<-"beta"
#	cat(paste0(method,"\n"))
	temp<-pCutoffMatrix(x1=q1_l:q1_u,x0=q0_l:q0_u, n=n, phi=phi0, w=w,k=k,alpha=alpha,method=method,bigCount=bigCount)
	
	a<-0
	if (returnDetail) {
		b<-matrix(1,nrow=length(q1_l:q1_u),ncol=length(q0_l:q0_u))
		X1<-temp$X1
		Y1<-temp$Y1
		X2<-temp$X2
		Y2<-temp$Y2
	}
	
	temp1<-pnbinom(q1_u,mu=(k*n*mu1), size=k*n/phi1)
	temp2<-pnbinom(q0_u,mu=(n*mu0), size=n/phi0)

	if (!is.na(temp$Y1[1])) {
		a<-sum((temp1-pnbinom(temp$Y1-1,mu=(k*n*mu1), size=k*n/phi1))*dnbinom(temp$X1,mu=(n*mu0), size=n/phi0))
	}
	if (!is.na(temp$Y2[1])) {
		a<-a+sum((temp2-pnbinom(temp$X2-1,mu=(n*mu0), size=n/phi0))*dnbinom(temp$Y2,mu=(k*n*mu1), size=k*n/phi1))
	}
	
	if (returnDetail) {
		#TODO: better b to result matrix
		if (!is.na(Y1[1])) {
			for (i in 1:length(Y1)) {
				b[(Y1[i]-q1_l+1):(q1_u-q1_l+1),(X1[i]-q0_l+1)]<-alpha
			}
		}
		if (!is.na(Y2[1])) {
			for (i in 1:length(Y2)) {
				b[(Y2[i]-q1_l+1),(X2[i]-q0_l+1):(q0_u-q0_l+1)]<-alpha
			}
		}
		colnames(b)<-q0_l:q0_u
		row.names(b)<-q1_l:q1_u
		
#		return(list(matrix=b,X1=X1,X2=X2,Y1=Y1,Y2=Y2,power=a-(1-beta)))
		return(list(matrix=b,X1=X1[X1!=q0_l & X1!=q0_u & Y1!=q1_l & Y1!=q1_u],X2=X2[X2!=q0_l & X2!=q0_u & Y2!=q1_l & Y2!=q1_u],Y1=Y1[X1!=q0_l & X1!=q0_u & Y1!=q1_l & Y1!=q1_u],Y2=Y2[X2!=q0_l & X2!=q0_u & Y2!=q1_l & Y2!=q1_u],power=a-(1-beta)))
	} else {
		return(a-(1-beta))
	}
}

pCutoffMatrix<-function(x1,x0, n, phi, w=1,k=1,alpha=0.05,method=c("nb","beta"),bigCount=900) {
	method<-match.arg(method)
	alphaOneSide<-alpha/2
	
#	y<-unique(as.vector(outer(x1,x0,'+')))
	y<-(min(x1)+min(x0)):(max(x1)+max(x0))
	if (method=="nb") {
		largeCountInd<-y>=(bigCount*4)
		if (any(largeCountInd)) { #Large count*4, beta approciate
			x0Max<-vector( "numeric",length=length(y))
			
			n1<-n
			n2<-k*n
			mu <- y[largeCountInd]/(n1+n2)
			alpha1 <- n1*mu/(1+phi*mu)
			alpha2 <- n2/n1*alpha1
			#d=(x0Cutoff+0.5)/y
			d<-qbeta(alphaOneSide,alpha1,alpha2)
			x0Max[largeCountInd]<-as.integer(d*y[largeCountInd]-0.5)
			
			#Others nb
			if (!all(largeCountInd)) {
				a2Fx<-generateA2FxR(max(y[!largeCountInd]),n,phi)
				yy<-lapply(y[!largeCountInd],function(x) nb_pvalue_store(x, n=n, phi=phi, w=w,k=k,a2Fx=a2Fx))
				x0Max[!largeCountInd]<-sapply(yy,function(x) cumsumBorder(x,alphaOneSide))
			}
		} else {
			a2Fx<-generateA2FxR(max(y),n,phi)
			yy<-lapply(y,function(x) nb_pvalue_store(x, n=n, phi=phi, w=w,k=k,a2Fx=a2Fx))
			x0Max<-sapply(yy,function(x) cumsumBorder(x,alphaOneSide))
		}
	} else { #beta app
		n1<-n
		n2<-k*n
		mu <- y/(n1+n2)
		alpha1 <- n1*mu/(1+phi*mu)
		alpha2 <- n2/n1*alpha1
		#d=(x0Cutoff+0.5)/y
		d<-qbeta(alphaOneSide,alpha1,alpha2)
		x0Max<-as.integer(d*y-0.5)
	}

	temp<-x0Max>=min(x0)
	x0Select<-x0Max[temp]
	x1Select<-(y-x0Max)[temp]
	
	temp<-x1Select<min(x1)
	x0Select[temp]<-x1Select[temp]+x0Select[temp]-min(x1)
	x1Select[temp]<-min(x1)
	temp<-x0Select>max(x0)
	x1Select[temp]<-x1Select[temp]+x0Select[temp]-max(x0)
	x0Select[temp]<-max(x0)
	temp<-x1Select<=max(x1)
	x1Select<-x1Select[temp]
	x0Select<-x0Select[temp]
	
	temp<-c(1,1+which(diff(x0Select)!=0))
	X1<-x0Select[temp]
	Y1<-x1Select[temp]
	
	if (method=="nb") {
		if (any(largeCountInd)) { #Large count*4, beta approciate
			x0Min<-vector( "numeric",length=length(y))
			
			#d=(x0Cutoff-0.5)/y
			d<-qbeta(alphaOneSide,alpha1,alpha2,lower.tail=FALSE)
			x0Min[largeCountInd]<-as.integer(d*y[largeCountInd]+0.5)+1
			
			if (!all(largeCountInd)) {
				x0Min[!largeCountInd]<-y[!largeCountInd]-sapply(yy,function(x) cumsumBorder(rev(x),alphaOneSide))
			}
		} else {
			x0Min<-y-sapply(yy,function(x) cumsumBorder(rev(x),alphaOneSide*sum(x)))
		}
	} else { #beta app
		#d=(x0Cutoff-0.5)/y
		d<-qbeta(alphaOneSide,alpha1,alpha2,lower.tail=FALSE)
		x0Min<-as.integer(d*y+0.5)+1
	}

	temp<-x0Min<=max(x0)
	x0Select<-x0Min[temp]
	x1Select<-(y-x0Min)[temp]
	
	temp<-x0Select<=min(x0)
	x1Select[temp]<-x1Select[temp]+x0Select[temp]-min(x0)
	x0Select[temp]<-min(x0)
	temp<-x1Select>=min(x1)
	x1Select<-x1Select[temp]
	x0Select<-x0Select[temp]
	temp<-x1Select>=max(x1)
	x0Select[temp]<-x1Select[temp]+x0Select[temp]-max(x1)
	x1Select[temp]<-max(x1)
	
	temp<-c(1,1+which(diff(x1Select)!=0))
	X2<-x0Select[temp]
	Y2<-x1Select[temp]
	
	return(list(X1=X1,Y1=Y1,X2=X2,Y2=Y2))
}