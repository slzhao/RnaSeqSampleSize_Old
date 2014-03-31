#' @useDynLib RnaSeqSampleSize generateA2Fx myDnbinom2

generateA2FxR <- function(a,n,phi)
	.Call("generateA2Fx",as.double(a),as.double(n),as.double(phi),PACKAGE="RnaSeqSampleSize")

myDnbinomR <- function(a, mu,size,a2Fx)
	.Call("myDnbinom2",as.double(a),as.double(mu),as.double(size),as.double(a2Fx[a+1]),PACKAGE="RnaSeqSampleSize")

est_pvalue_store<-function(X1X0Sum, n, phi, w=1,a2Fx=NULL) {
	a<-0:X1X0Sum
	if(w==1){
		azanom<-myDnbinomR(a, mu=X1X0Sum/2, size=n/phi,a2Fx=a2Fx)
	}else{
		anom<-dnbinom(a, mu=(w*X1X0Sum/(w+1)), size=n/phi)
		zanom<-dnbinom( X1X0Sum-a, mu=(X1X0Sum/(w+1)), size=n/phi)
		azanom<-anom * zanom
	}
	yy<-azanom/sum(azanom)
	return(yy)
}

est_pvalue_store_get_power3<-function(x1,x0,yy) {
	yy<-yy[[(x1+x0+1)]]
	if (x1>=x0) {
		pvalue<-2* sum(yy[(x1+1):length(yy)])
	} else {
		pvalue<-2*sum(yy[1:(x1+1)])
	}
	return(min(pvalue, 1))
}

est_pvalue<-function(x1,x0,n, phi, w=1,a2Fx=NULL) {
	yy<-est_pvalue_store(X1X0Sum=x1+x0,n=n, phi=phi, w=w,a2Fx=a2Fx)
	if (x1>=x0) {
		pvalue<-2* sum(yy[(x1+1):length(yy)])
	} else {
		pvalue<-2*sum(yy[1:(x1+1)])
	}
	return(min(pvalue, 1))
}

##' est_power
##' 
##' A function to estitamete the power for differential expression analysis of RNA-seq data.
##' 
##' A function to estitamete the power for differential expression analysis of RNA-seq data.
##' 
##' @param n Numer of samples.
##' @param alpha alpha level.
##' @param ... other parameters for est_power_root function. 
##' @inheritParams sample_size
##' @export
##' @examples n<-20;valpha<-0.05;vw<-1.0;vrho<-1.5;vlambda0<-5;vphi_0<-0.5
##' est_power(n=n, alpha=valpha, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)
est_power<-function(n, w=1, rho=2, lambda0=5, phi0=1,alpha=0.05,...){
	power<-est_power_root(n=n, w=w, rho=rho, lambda0=lambda0, phi_0=phi0, alpha=alpha,...)
	return(power+0.8)
}

est_power_root<-function(n, w=1, rho=2.0, lambda0=5, phi_0=1, beta=0.2, alpha=0.05, error=0.001,returnDetail=F){
	mu0<-lambda0
	mu1<-mu0*(rho*w)
	phi_1<-phi_0
	q0_u<-qnbinom(1-error, size=n/phi_0, mu=n*mu0)
	q0_l<-qnbinom(error, size=n/phi_0, mu=n*mu0)
	q1_u<-qnbinom(1-error, size=n/phi_1, mu=n*mu1)
	q1_l<-qnbinom(error, size=n/phi_1, mu=n*mu1)
	a2Fx<-generateA2FxR(max(q0_u+q1_u,q0_l+q1_l),n,phi_0)
	
	a<-0
	if (returnDetail) {
		b<-matrix(1,nrow=length(q1_l:q1_u),ncol=length(q0_l:q0_u))
		X1<-NULL
		Y1<-NULL
		X2<-NULL
		Y2<-NULL
	}
	
	temp1<-dnbinom(q1_l:q1_u, mu=(n*mu1), size=n/phi_1)
	temp2<-dnbinom(q0_l:q0_u, mu=(n*mu0), size=n/phi_0)
	if (max(q0_u,q0_l,q1_u,q1_l)>=10000) { #Method2, doesn't store every pvalue but do it every time
		getPvalue<-function(x1,x0,...) {
			est_pvalue(x1,x0,n=n, phi=phi_0, w=w,a2Fx=a2Fx,...)
		}
	} else { #Method1, store every pvalue
		yMin<-min(q1_l,q1_u)+min(q0_l,q0_u)
		yMax<-max(q1_l,q1_u)+max(q0_l,q0_u)
		yyStore<-list()
		for (y in yMin:yMax) {
			yyStore[[y+1]]<-est_pvalue_store(y,n=n, phi=phi_0, w=w,a2Fx=a2Fx)
		}
		getPvalue<-function(x1,x0,yy=yyStore,...) {
			est_pvalue_store_get_power3(x1,x0,yy=yyStore,...)
		}
	}
	
	aNRow<-length(q1_l:q1_u)
	aNCol<-length(q0_l:q0_u)
	q0_l_loop<-q0_l
	q0_u_loop<-q0_u
	q1_l_loop<-q1_l
	q1_u_loop<-q1_u
	for (x in q0_l_loop:q0_u_loop) {
		for (y in q1_l_loop:q1_u_loop) {
			if (x>y) {next;}
			temp<-getPvalue(x1=y,x0=x)
			if (temp<=alpha) {
				a<-a+sum(temp1[(y-q1_l+1):aNRow]*temp2[(x-q0_l+1)])
				q1_l_loop<-y
				if (returnDetail) {
					b[(y-q1_l+1):aNRow,(x-q0_l+1)]<-temp
					if (y!=q1_l) {
						X1<-c(X1,x)
						Y1<-c(Y1,y)
					}
				}
				break;
			} else if(returnDetail) {
				b[(y-q1_l+1),(x-q0_l+1)]<-temp
			}
		}
	}
	q0_l_loop<-q0_l
	q0_u_loop<-q0_u
	q1_l_loop<-q1_l
	q1_u_loop<-q1_u
	for (y in q1_l_loop:q1_u_loop) {
		for (x in q0_l_loop:q0_u_loop) {
			if (x<y) {next;}
			temp<-getPvalue(x1=y,x0=x)
			if (temp<=alpha) {
				a<-a+sum(temp1[(y-q1_l+1)]*temp2[(x-q0_l+1):aNCol])
				q0_l_loop<-x
				if (returnDetail) {
					b[(y-q1_l+1),(x-q0_l+1):aNCol]<-temp
					if (x!=q0_l) {
						X2<-c(X2,x)
						Y2<-c(Y2,y)
					}
				}
				break;
			} else if(returnDetail) {
				b[(y-q1_l+1),(x-q0_l+1)]<-temp
			}
		}
	}
	if (returnDetail) {
		colnames(b)<-q0_l:q0_u
		row.names(b)<-q1_l:q1_u
	} 
	if (returnDetail) {
		return(list(matrix=b,X1=X1,X2=X2,Y1=Y1,Y2=Y2,power=a-(1-beta)))
	} else {
		return(a-(1-beta))
	}
}

uniroot.integer<-function (f, interval, lower = min(interval), upper = max(interval), 
		step.power = 6, step.up = TRUE, pos.side = FALSE, print.steps = FALSE, 
		maxiter = 1000, ...) 
{
	stored<-NULL
	iter <- 0
	if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
			upper) 
		stop("lower < upper  is not fulfilled")
	if (lower == -Inf && step.up == TRUE) 
		stop("lower cannot be -Inf when step.up=TRUE")
	if (upper == Inf && step.up == FALSE) 
		stop("upper cannot be Inf when step.up=FALSE")
	if (step.up) {
		f.old <- f(lower, ...)
		iter <- iter + 1
		sign <- 1
		xold <- lower
	}
	else {
		f.old <- f(upper, ...)
		iter <- iter + 1
		sign <- -1
		xold <- upper
	}
	if (print.steps) {
		print(paste("x=", xold, " f(x)=", f.old))
	}
	stored<-rbind(stored,c(xold,f.old))
	ever.switched <- FALSE
	tried.extreme <- FALSE
	while (step.power > -1) {
		xnew <- xold + sign * 2^step.power
		if ((step.up & xnew < upper) || (!step.up & xnew > lower)) {
			f.new <- f(xnew, ...)
			iter <- iter + 1
			if (print.steps) {
				print(paste("x=", xnew, " f(x)=", f.new))
			}
			stored<-rbind(stored,c(xnew,f.new))
		}
		else {
			xnew <- xold
			f.new <- f.old
			step.power <- step.power - 1
			if (tried.extreme == FALSE) {
				if (step.up) {
					f.extreme <- f(upper, ...)
					iter <- iter + 1
					x.extreme <- upper
				}
				else {
					f.extreme <- f(lower, ...)
					iter <- iter + 1
					x.extreme <- lower
				}
				tried.extreme <- TRUE
				xswitch <- x.extreme
				f.switch <- f.extreme
				if (print.steps) {
					print(paste("x=", x.extreme, " f(x)=", f.extreme))
				}
				stored<-rbind(stored,c(x.extreme,f.extreme))
				if (f.old * f.extreme >= 0) {
					stop("f() at extremes not of opposite sign")
				}
			}
		}
		if (f.old * f.new < 0) {
			sign <- sign * (-1)
			ever.switched <- TRUE
			xswitch <- xold
			f.switch <- f.old
		}
		if (ever.switched) {
			step.power <- step.power - 1
			if (step.power == -1) {
				(break)()
			}
		}
		xold <- xnew
		f.old <- f.new
	}
	if (pos.side) {
		root <- ifelse(f.new > 0, xnew, xswitch)
		f.root <- ifelse(f.new > 0, f.new, f.switch)
	}
	else {
		root <- ifelse(f.new < 0, xnew, xswitch)
		f.root <- ifelse(f.new < 0, f.new, f.switch)
	}
	colnames(stored)<-c("N","Power")
	return(list(iter = iter, f.root = f.root, root = root,process=stored))
}

find_near<-function(allData,input,method=c("smaller","larger")) {
	if (input %in% allData) {
		return(input)
	}
	method<- if (missing(method)) "smaller" else match.arg(method)
	if (method=="smaller") {
		allData<-allData[which(allData<=input)]
	} else {
		allData<-allData[which(allData>=input)]
	}
	if (length(allData>0)) {
		return(allData[which.min(abs(allData-input))])
	} else {return(NA)}
}

find_near_N<-function(nTable,rho,phi0,lambda0,f,power,m,m1,showMessage=showMessage) {
	rho_near<-find_near(nTable[,1],rho,"larger")
	phi0_near<-find_near(nTable[,2],phi0,"smaller")
	lambda0_near<-find_near(nTable[,3],lambda0,"larger")
	f_near<-find_near(nTable[,4],f,"larger")
	power_near<-find_near(nTable[,5],power,"larger")
	m_near<-find_near(nTable[,6],m,"smaller")
	m1_near<-find_near(nTable[,7],m1,"larger")
	n_near<-nTable[which(nTable[,1]==rho_near & nTable[,2]==phi0_near & nTable[,3]==lambda0_near & nTable[,4]==f_near & nTable[,5]==power_near & nTable[,6]==m_near & nTable[,7]==m1_near),8]
	if (length(n_near)==0) {
		n_near<-1
		if (showMessage) {
			cat("RnaSeqSampleSize can't find nearest sample size in stored N table. Will estimate sample size N from N=1.\n")
		}
	} else if (as.character(rho_near)==as.character(rho) & as.character(phi0_near)==as.character(phi0) & as.character(lambda0_near)==as.character(lambda0) & as.character(f_near)==as.character(f) & as.character(power_near)==as.character(power) & as.character(m_near)==as.character(m) & as.character(m1_near)==as.character(m1)) {
		n_near<-paste("M",n_near,sep="")
	} else {
		estTime<-as.integer(nTable[which(nTable[,1]==rho_near & nTable[,2]==phi0_near & nTable[,3]==lambda0_near & nTable[,4]==f_near & nTable[,5]==power_near & nTable[,6]==m_near & nTable[,7]==m1_near),9])
		if (showMessage) {
			if (length(estTime)==0) {
				cat("RnaSeqSampleSize can't find nearest sample size in stored N table. Will estimate sample size N from N=1.\n")
			} else {
				cat(paste("RnaSeqSampleSize will use about ",estTime," seconds to estimate sample number N.\n",sep=""))
			}
			flush.console()
		}
	}
	return(n_near)
}


est_power_model<-function(n, w=1, rho=2, lambda0=5, phi_0=1, beta=0.2, alpha=0.05, error=0.001,showMessage=F){
	mu0<-lambda0
	mu1<-mu0*(rho*w)
	phi_1<-phi_0
	
	model<-generateModel(lambda0=c(1,5,10),n=n, w=w, rho=rho, phi_0=phi_0,showMessage=F)
	modelX1<-model[1]
	modelIntercept1<-model[2]
	modelX2<-model[3]
	modelIntercept2<-model[4]
	
	a<-0
	q0_u<-qnbinom(1-error, size=n/phi_0, mu=n*mu0)
	q0_l<-qnbinom(error, size=n/phi_0, mu=n*mu0)
	q1_u<-qnbinom(1-error, size=n/phi_1, mu=n*mu1)
	q1_l<-qnbinom(error, size=n/phi_1, mu=n*mu1)
	dnbinomQ1<-dnbinom(q1_l:q1_u, mu=(n*mu1), size=n/phi_1)
	pnbinomQ1<-pnbinom((q1_l-1):q1_u, mu=(n*mu1), size=n/phi_1)
	dnbinomQ0<-dnbinom(q0_l:q0_u, mu=(n*mu0), size=n/phi_0)
	pnbinomQ0<-pnbinom((q0_l-1):q0_u, mu=(n*mu0), size=n/phi_0)
	aNRow<-length(q1_l:q1_u)
	aNCol<-length(q0_l:q0_u)
	
	if (!is.na(modelX1)) {
		y<-round(modelX1*(q0_l:q0_u)+modelIntercept1)
		for (i in 1:aNCol) {
			if (y[i]>q1_u) {
				break;
			} else if (y[i]<q1_l) {
				y[i]<-q1_l
			}
			a<-a+(pnbinomQ1[aNRow]-pnbinomQ1[(y[i]-q1_l+1)])*dnbinomQ0[i]
		}
	}
	
	if (!is.na(modelX2)) {
		y<-round(modelX2*(q0_l:q0_u)+modelIntercept2)
		for (i in 1:aNCol) {
			if (y[i]>q1_u) {
				y[i]<-q1_u
			} else if (y[i]<q1_l) {
				break;
			}
			a<-a+(pnbinomQ1[aNRow]-pnbinomQ1[(y[i]-q1_l+1)])*dnbinomQ0[i]
		}
	}
	
	return(a-(1-beta))
}

generateModel<-function(lambda0=c(1,5,10,20),n, w=1, rho=2.0, phi_0=1, alpha=0.05, error=0.001,showMessage=F) {
	result<-NULL
	X1<-NULL
	X2<-NULL
	Y1<-NULL
	Y2<-NULL
	for (i in 1:length(lambda0)) {		
		temp<-est_power_root(n=n, w=w, rho=rho, lambda0=lambda0[i], phi_0=phi_0, alpha=alpha, error=error,returnDetail=T)
		result[[i]]<-temp[["matrix"]]
		X1<-c(X1,temp[["X1"]])
		X2<-c(X2,temp[["X2"]])
		Y1<-c(Y1,temp[["Y1"]])
		Y2<-c(Y2,temp[["Y2"]])
	}
	
	if (showMessage) {
		maxRange<-max(c(X1,X2,Y1,Y2))+1
		plot(c(0,maxRange),c(0,maxRange),type="l",xlab="q0",ylab="q1",main=paste("n=",n,"; rho=",rho,"; phi=",phi_0,sep=""),lty=2)
		
		for (i in 1:length(result)) {
			temp1<-range(as.integer(colnames(result[[i]])))[c(1,2,2,1)]
			temp2<-range(as.integer(row.names(result[[i]])))[c(2,2,1,1)]
			temp3<-which(result[[i]]<=alpha,arr.ind=T)
			
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
##' @param rho minimum fold changes for prognostic genes between two groups.
##' @param lambda0 Average read counts for prognostic genes.
##' @param phi0 Dispersion for prognostic genes.
##' @param nTable A table included the estitameted sample numbers (N) with different parameters.
##' @param showMessage Logical. Display the message in the estimation process.
##' @param storeProcess Logical. Store the power and n in sample size estimation process.
##' @export
##' @examples #Input initial parameters. These parameters were stored in our N table so you will get results in one second.
##' vm<-10000;vm1<-100;vpower<-0.8;vf<-0.05;vw<-1.0;vrho<-2.0;vlambda0<-5;vphi_0<-2
##' sample_size(m=vm, m1=vm1, power=vpower, f=vf, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)
##' #Estitamete sample size for parameters not stored, may use one minutes.
##' vm<-30000;vm1<-1000;vpower<-0.8;vf<-0.05;vw<-1.0;vrho<-1.5;vlambda0<-5;vphi_0<-0.5
##' sample_size(m=vm, m1=vm1, power=vpower, f=vf, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)
sample_size<-function(m=10000, m1=100, power=0.8, f=0.1, w=1, rho=2, lambda0=5, phi0=1,nTable=system.file("extdata", "nTable.csv", package = "RnaSeqSampleSize"),showMessage=F,storeProcess=F){
	lambda0AppCut<-20
	r1<-m1 * power
	beta<-1-power
	step.power<-5
	
	alpha_star<-r1*f/((m-m1)*(1-f))
	z_alpha<-qnorm(1-alpha_star/2, lower.tail=T)
	z_beta<-qnorm(power, lower.tail=T)
	n_w<-( ( z_alpha + z_beta )^2* (1 + rho/w +phi0*lambda0*(1+rho^2)) )/ ( (rho-1)^2*lambda0 )
	end.point<-round(n_w)+30
	
	if (storeProcess) {
		start.point<-1
	} else {
		if (file.exists(nTable)) {
			nTable<-read.csv(nTable,header=T)
			n_near<-find_near_N(nTable,rho=rho,phi0=phi0,lambda0=lambda0,f=f,m=m, m1=m1, power=power,showMessage=showMessage)
			if (n_near!="") {
				if (length(grep("M",n_near))!=0) {
					n_Exact<-as.numeric(gsub("M","",n_near))
					return(n_Exact)
				} else {
					start.point<-n_near
				}
			} else {
				start.point<-1
			}
		} else {
			start.point<-1
			if (showMessage) {
				cat(paste("nTable can't be found or not specified. RnaSeqSampleSize will estimate sample size N from N=1.\n",sep=""))
			}
		}
	}

	if (lambda0>=lambda0AppCut) {
		p1<-est_power_model(n=start.point,w=w, rho=rho, lambda0=lambda0, phi_0=phi0, beta=beta, alpha=alpha_star)
	} else {
		p1<-est_power_root(n=start.point,w=w, rho=rho, lambda0=lambda0, phi_0=phi0, beta=beta, alpha=alpha_star)
	}
	if (p1>0) {
		return(start.point)
	} else {
		step.power<-6
		if (lambda0>=lambda0AppCut) {
			n_Exact<-uniroot.integer(est_power_model, c(start.point, end.point), w=w, rho=rho, 
					lambda0=lambda0, phi_0=phi0, beta=beta, alpha=alpha_star, pos.side=T,
					step.up=T, step.power=step.power,print.steps=showMessage)
		} else {
			n_Exact<-uniroot.integer(est_power_root, c(start.point, end.point), w=w, rho=rho, 
					lambda0=lambda0, phi_0=phi0, beta=beta, alpha=alpha_star, pos.side=T,
					step.up=T, step.power=step.power,print.steps=showMessage) 
		}
		if (storeProcess) {
			n_Exact$process<-n_Exact$process[order(n_Exact$process[,1]),]
			n_Exact$process[,2]<-n_Exact$process[,2]+power
			n_Exact$parameters<-c(power,m, m1,  f, w, rho, lambda0, phi0)
			names(n_Exact$parameters)<-c("power","m", "m1",  "fdr", "w", "rho", "lambda0", "phi0")
		} else {
			n_Exact<-n_Exact$root
		}
		return(n_Exact)
	}
}

##' plot_power_curve
##' 
##' A function to plot power curves based on the result of \code{\link{sample_size}} or \code{\link{est_power_curve}} function.
##' 
##'
##' 
##' @param result the result of \code{\link{sample_size}} or \code{\link{est_power_curve}} function. The storeProcess parameter should be set as True when performing \code{\link{sample_size}} function. If you want to plot more than one curves in the same figure, the results from \code{\link{sample_size}} function should first be combined into a new list. At most five curves were allowed in one figure.
##' @param cexLegend the cex for legend.
##' @export
##' @examples #One power curve
##' result1<-sample_size(rho=2,phi0=1,lambda0=1,f=0.01,power=0.8,m=20000,m1=500,showMessage=T,storeProcess=T)
##' plot_power_curve(result1)
##' #More than one curves
##' result2<-sample_size(rho=4,phi0=1,lambda0=1,f=0.01,power=0.8,m=20000,m1=500,showMessage=T,storeProcess=T)
##' result<-list(result1,result2)
##' plot_power_curve(result)
plot_power_curve<-function(result,cexLegend=1) {
	if (identical(names(result),c("iter","f.root","root","process","parameters")) || identical(names(result),c("process","parameters","power"))) { #Only one result
		plot(result$process[,1],result$process[,2],type="b",xlab="N",ylab="Power",pch=16,lwd=3,las=1,cex=1.5,main="Power Curve",col="red")
		legend("bottomright",legend=paste(paste(names(result$parameters),result$parameters,sep="="),collapse=";"),bty="n",text.col="red",cex=0.9)
		abline(h= result$parameters["power"],lty=2,col="grey")
	} else { #more than one result
		if (length(result)>5) { #at most 5 curves were allowed
			result<-result[1:5]
		}
		resultRange<-apply(sapply(result,function(x) (x$process)[nrow(x$process),]),1,max)
		plot(c(0,resultRange[1]),c(0,resultRange[2]),type="n",xlab="N",ylab="Power",las=1,main="Power Curve")
		col<-c("brown1","steelblue2","mediumpurple2","seagreen3","lightgoldenrod")
		legendEach<-""
		for (x in 1:length(result)) {
			resultEach<-result[[x]]
			lines(resultEach$process[,1],resultEach$process[,2],type="b",pch=16,lwd=3,cex=1.5,col=col[x])
			legendEach[x]<-paste(paste(names(resultEach$parameters),resultEach$parameters,sep="="),collapse=";")
		}
		legend("bottomright",legend=legendEach,bty="n",text.col=col[1:length(result)],cex=cexLegend)
	}
}

##' est_power_curve
##' 
##' A function to estitamete the power curve for differential expression analysis of RNA-seq data.
##' 
##'
##' 
##' @param n Numer of samples.
##' @param alpha alpha level.
##' @param ... other parameters for est_power function. 
##' @inheritParams est_power
##' @export
##' @examples result<-est_power_curve(n=20, alpha=0.05, w=1, rho=1.5, lambda0=5, phi0=0.5)
##' plot_power_curve(result)
est_power_curve<-function(n, w=1, rho=2, lambda0=5, phi0=1,alpha=0.05,...) {
	if (n<=10) {
		sampleSizeList<-1:n
	} else if (! n %% 5) {
		sampleSizeList<-c(1,seq(5,n,by=5))
	} else {
		sampleSizeList<-c(1,seq(5,n,by=5),n)
	}
	powerList<-NULL
	for (i in 1:length(sampleSizeList)) {
		powerList[i]<-est_power(n=sampleSizeList[i], w=w, rho=rho, lambda0=lambda0, phi0=phi0, alpha=alpha)
	}
	process<-cbind(N=sampleSizeList,Power=powerList)
	parameters<-c(n, alpha, w, rho, lambda0, phi0)
	names(parameters)<-c("n","alpha", "w", "rho", "lambda0", "phi0")
	return(list(process=process,parameters=parameters,power=process[nrow(process),2]))
}

