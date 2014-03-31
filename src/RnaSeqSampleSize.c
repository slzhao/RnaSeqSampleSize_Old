#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "nmath.h"
#include "dpq.h"
/* #include "bd0.c" */
/* #include "stirlerr.c" */

#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

SEXP generateA2Fx(SEXP a,SEXP n,SEXP phi){
	int i, aV;
	aV = INTEGER_VALUE(a);
	double nDividePhi=REAL(n)[0]/REAL(phi)[0];
	SEXP result = PROTECT(allocVector(REALSXP, (aV+1)));
	double *xresult=REAL(result);
	xresult[0]=0;
	for (i = 1; i <= aV; i++) xresult[i]=stirlerr(i+(nDividePhi)) - stirlerr(nDividePhi) - stirlerr(i)-0.5*(M_LN_2PI + log(nDividePhi) + log1p(- (nDividePhi)/(i+(nDividePhi))))+log(nDividePhi/(i+nDividePhi));
	UNPROTECT(1);
	return result;
}

double attribute_hidden
dbinom_raw1(double x, double n, double p, double q, double aFxValue)
{
    double lc;
    if (q == 0) return(0);
    /* n*p or n*q can underflow to zero if n and p or q are small.  This
       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
    lc = aFxValue - bd0(x,n*p) - bd0(n-x,n*q);
    return lc;
}


double attribute_hidden
dbinom_raw2(double x, double n1,double n2, double p, double q, double aFxValue)
{
    double lc;
    if (q == 0) return(0);
    /* n*p or n*q can underflow to zero if n and p or q are small.  This
       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
    lc = aFxValue - bd0(x,n1*p) - bd0(n1-x,n1*q)- bd0(x,n2*p) - bd0(n2-x,n2*q);
    return exp(lc);
}


double dnbinom_mu2(double x1,double x2, double size, double mu,double aFxValue)
{
    /* originally, just set  prob :=  size / (size + mu)  and called dbinom_raw(),
     * but that suffers from cancellation when   mu << size  */
    double ans;

 /*   x = R_forceint(x);*/
    if(x1 == 0 || x2 ==0) {
		double ans1,ans2;
		if (x1 == 0) {/* be accurate, both for n << mu, and n >> mu :*/
			ans1=size * (size < mu ? log(size/(size+mu)) : log1p(- mu/(size+mu)));
		} 
		else {/*We can use the aFxValue for x1+x2 to caculate only x1 or x2 because when one of them ==0, its aFxValue==0*/
			ans1=dbinom_raw1(size, x1+size, size/(size+mu), mu/(size+mu), aFxValue);
		}
		if (x2 == 0) {
			ans2=size * (size < mu ? log(size/(size+mu)) : log1p(- mu/(size+mu)));
		}
		else {
			ans2=dbinom_raw1(size, x2+size, size/(size+mu), mu/(size+mu), aFxValue);
		}
		return exp(ans1+ans2);
	}
	
    /* else: no unnecessary cancellation inside dbinom_raw, when
     * x_ = size and n_ = x+size are so close that n_ - x_ loses accuracy
     */
    ans = dbinom_raw2(size, x1+size, x2+size, size/(size+mu), mu/(size+mu), aFxValue);
	return(ans);
/*    p = ((double)size)/(size+x);
 *   return(p * ans);*/
}


SEXP myDnbinom2(SEXP a, SEXP mu, SEXP size, SEXP aFxValue){
	int i, na;
	na = length(a);
	SEXP result = PROTECT(allocVector(REALSXP, na));
	double *xa=REAL(a);
	double *xmu=REAL(mu),*xsize=REAL(size), *xaFxValue=REAL(aFxValue),*xresult=REAL(result);
	for (i = 0; i < na/2+1; i++) {
		xresult[i] = dnbinom_mu2(xa[i],xa[(na-i-1)], xsize[0], xmu[0],(xaFxValue[i]+xaFxValue[na-i-1])); /*Please note the first parameter we used i insteat of xi[i], which is the really value*/
		xresult[(na-i-1)] = xresult[i];
	}
	UNPROTECT(1);
	return result;
}

