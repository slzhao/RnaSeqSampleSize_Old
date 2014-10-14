#' @useDynLib RnaSeqSampleSize generateA2Fx myDnbinom2
#' @importFrom Rcpp sourceCpp
#' @importFrom KEGGREST keggLink


generateA2FxR <- function(a,n,phi)
	.Call("generateA2Fx",as.double(a),as.double(n),as.double(phi),PACKAGE="RnaSeqSampleSize")

myDnbinomR <- function(a, mu,size,a2Fx)
	.Call("myDnbinom2",as.double(a),as.double(mu),as.double(size),as.double(a2Fx[a+1]),PACKAGE="RnaSeqSampleSize")

nb_pvalue_store<-function(X1X0Sum, n, phi, w=1,k=1,a2Fx=NULL) {
	a<-0:X1X0Sum
	weight<-k*w/(k*w+1)
	if(w==1 && k==1){
		azanom<-myDnbinomR(a, mu=X1X0Sum/2, size=n/phi,a2Fx=a2Fx)
#		anom<-dnbinom(a, mu=(X1X0Sum/2), size=n/phi)
#		zanom<-rev(anom)
#		azanom<-anom * zanom
	}else{
		anom<-dnbinom(a, mu=X1X0Sum*(weight), size=n*k/phi)
		zanom<-dnbinom(X1X0Sum-a, mu=X1X0Sum*(1-weight),size=n/phi)
		azanom<-anom * zanom
	}
	yy<-azanom/sum(azanom)
	return(yy)
}





