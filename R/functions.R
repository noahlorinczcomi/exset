#' Monte Carlo Integration
#'
#' This internal function allows you to integrate over a half plane in the joint space of independent weighted chi-squares
#' @param lam eigenvalues of positive-definite LD matrix
#' @param niter number of Monte Carlo replicates used to approximate integral
#' @param alpha type I error rate
#' @keywords integration
#' @export
#' @examples
#' mc_integrate()
mc_integrate=function(lam,niter=10000,alpha=0.05) {
  lam=sort(lam)
  p=length(lam)
  X=matrix(rchisq(niter*p,1),niter,p)
  X=X*matrix(lam,niter,p,byrow=TRUE)
  r=rowSums(X)
  tau=quantile(r,prob=1-alpha)
  return(tau)
}

#' Exact Gene-Based Testing
#'
#' This function allows you to perform an exact hypothesis test of the gene-based null hypothesis
#' @param ld_eigenvalues eigenvalues of positive-definite LD matrix
#' @param iterations number of Monte Carlo replicates used to approximate integral
#' @param type1_error corrected type I error rate
#' @keywords integration
#' @export
#' @examples
#' exset()
exset=function(statistic,ld_eigenvalues,iterations=1e5,type1_error=0.05/12727) {
  crit=mc_integrate(ld_eigenvalues,niter=iterations,alpha=type1_error)
  boo=statistic>crit
  dec=boo
  dec[boo]='reject H0'
  dec[!boo]='fail to reject H0'
  return(dec)
}
