est_pvalue_store<-function(X1X0Sum, n, phi, w=1) {
	a<-0:X1X0Sum
	if(w==1){
		anom<-dnbinom(a, mu=X1X0Sum/2, size=n/phi)
		zanom<-rev(anom)
	}else{
		anom<-dnbinom(a, mu=(w*X1X0Sum/(w+1)), size=n/phi)
		zanom<-dnbinom( X1X0Sum-a, mu=(X1X0Sum/(w+1)), size=n/phi)
	}
	azanom<-anom * zanom
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

est_pvalue<-function(x1,x0,n, phi, w=1) {
	yy<-est_pvalue_store(X1X0Sum=x1+x0,n=n, phi=phi, w=w)
	if (x1>=x0) {
		pvalue<-2* sum(yy[(x1+1):length(yy)])
	} else {
		pvalue<-2*sum(yy[1:(x1+1)])
	}
	return(min(pvalue, 1))
}

est_power3<-function(n, w=1, rho=2.0, lambda0=5, phi_0=1, beta=0.2, alpha=0.05, error=0.001){
	mu0<-lambda0
	mu1<-mu0*(rho*w)
	phi_1<-phi_0
	q0_u<-qnbinom(1-error, size=n/phi_0, mu=n*mu0)
	q0_l<-qnbinom(error, size=n/phi_0, mu=n*mu0)
	q1_u<-qnbinom(1-error, size=n/phi_1, mu=n*mu1)
	q1_l<-qnbinom(error, size=n/phi_1, mu=n*mu1)
	
	a<-matrix(nrow=length(q1_l:q1_u),ncol=length(q0_l:q0_u),F)
	if (max(q0_u,q0_l,q1_u,q1_l)>=10000) { #Method2, doesn't store every pvalue but do it every time
		getPvalue<-function(n=n, phi=phi_0, w=w,...) {
			est_pvalue(n=n, phi=phi_0, w=w,...)
		}
	} else { #Method1, store every pvalue
		yMin<-min(q1_l,q1_u)+min(q0_l,q0_u)
		yMax<-max(q1_l,q1_u)+max(q0_l,q0_u)
		yyStore<-list()
		for (y in yMin:yMax) {
			yyStore[[y+1]]<-est_pvalue_store(y,n=n, phi=phi_0, w=w)
		}
		getPvalue<-function(yy=yyStore,...) {
			est_pvalue_store_get_power3(yy=yyStore,...)
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
			if (x>y) {break;}
			temp<-getPvalue(x1=y,x0=x)
			if (temp<=alpha) {
				a[(y-q1_l+1):aNRow,(x-q0_l+1)]<-T
				q1_l_loop<-y
				break;
			}
		}
	}
	q0_l_loop<-q0_l
	q0_u_loop<-q0_u
	q1_l_loop<-q1_l
	q1_u_loop<-q1_u
	for (y in q1_l_loop:q1_u_loop) {
		for (x in q0_l_loop:q0_u_loop) {
			if (x<y) {break;}
			temp<-getPvalue(x1=y,x0=x)
			if (temp<=alpha) {
				a[(y-q1_l+1),(x-q0_l+1):aNCol]<-T
				q0_l_loop<-x
				break;
			}
		}
	}
	temp<-dnbinom(q1_l:q1_u, mu=(n*mu1), size=n/phi_1)
	b<-apply(a, 2, function(x) x*temp )
	temp<-dnbinom(q0_l:q0_u, mu=(n*mu0), size=n/phi_0)
	power<-sum( apply(b, 1, function(x) x*temp ) )
	return(power-(1-beta))
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
##' @param phi_0 Dispersion for prognostic genes.
##' @param nTable A table included the estitameted sample numbers (N) with different parameters.
##' @param showMessage Logical. Display the message in the estimation process.
##' @export
##' @examples #Input initial parameters. These parameters were stored in our N table so you will get results in one second.
##' vm<-10000;vm1<-100;vpower<-0.8;vf<-0.05;vw<-1.0;vrho<-2.0;vlambda0<-5;vphi_0<-2
##' sample_size(m=vm, m1=vm1, power=vpower, f=vf, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)
##' #Estitamete sample size for parameters not stored, may use two minutes.
##' vm<-30000;vm1<-1000;vpower<-0.8;vf<-0.05;vw<-1.0;vrho<-1.5;vlambda0<-5;vphi_0<-0.5
##' sample_size(m=vm, m1=vm1, power=vpower, f=vf, w=vw, rho=vrho, lambda0=vlambda0, phi_0=vphi_0)
sample_size<-function(m=10000, m1=100, power=0.8, f=0.1, w=1, rho=2, lambda0=5, phi_0=1,nTable=system.file("extdata", "nTable.csv", package = "RnaSeqSampleSize"),showMessage=T){
	r1<-m1 * power
	beta<-1-power
	
	alpha_star<-r1*f/((m-m1)*(1-f))
	z_alpha<-qnorm(1-alpha_star/2, lower.tail=T)
	z_beta<-qnorm(power, lower.tail=T)
	n_w<-( ( z_alpha + z_beta )^2* (1 + rho/w +phi_0*lambda0*(1+rho^2)) )/ ( (rho-1)^2*lambda0 )
	end.point<-round(n_w)+30
#	require("ssanv")
	
	if (file.exists(nTable)) {
		nTable<-read.csv(nTable,header=T)
		n_near<-find_near_N(nTable,rho=rho,phi0=phi_0,lambda0=lambda0,f=f,m=m, m1=m1, power=power,showMessage=showMessage)
		if (length(n_near)!=0) {
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
	p1<-est_power3(n=start.point,w=w, rho=rho, lambda0=lambda0, phi_0=phi_0, beta=beta, alpha=alpha_star)
	if (p1>0) {
		return(start.point)
	} else {
		n_Exact<-uniroot.integer(est_power3, c(start.point, end.point), w=w, rho=rho, 
				lambda0=lambda0, phi_0=phi_0, beta=beta, alpha=alpha_star, pos.side=T,
				step.up=T, step.power=5,print.steps=showMessage)$root  
		return(n_Exact)
	}
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
		n_near<-""
		if (showMessage) {
			cat("RnaSeqSampleSize can't find nearest sample size in stored N table. Will estimate sample size N from N=1.\n")
		}
	}	else if (as.character(rho_near)==as.character(rho) & as.character(phi0_near)==as.character(phi0) & as.character(lambda0_near)==as.character(lambda0) & as.character(f_near)==as.character(f) & as.character(power_near)==as.character(power) & as.character(m_near)==as.character(m) & as.character(m1_near)==as.character(m1)) {
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
