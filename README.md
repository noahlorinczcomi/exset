# Overview
This repository contains a brief introduction to an **ex**act **set**-based gene association test statistic (**exset**). 

Existing methods to test a gene-based association null hypothesis using the sum of SNP chi-squares in a gene-specific set rely on permutation, simulation, or approximation since the null distribution of the test statistic isn't known. This can lead to uncontrolled Type I errors and enormous computational burden. We overcome these computational and theoretical limitations by providing a theoretically exact hypothesis test which does not require stating the null distribution of the testing statistic. 

The principle of our method is to find the null region in the joint space of a random vector whose 2-norm is distributionally equivalent to the standard gene-based association test statistic, then to compare the standard test statistic to the critical value defining this null region. Please read our [Technical Note](https://github.com/noahlorinczcomi/exset/blob/main/Technical%20Note.pdf) for more details, and see our [simulation code](https://github.com/noahlorinczcomi/exset/blob/main/plot_and_simulation_code.R) to replicate all of our results.

# Implementation
We assume that researchers have gene-based association test statistics and intend to test their null hypotheses of no association with a phenotype. The [GenT Shiny application](https://nlorinczcomi.shinyapps.io/gent/) provides these statistics for over 50 phenotypes, along with quantities related to the linkage disequilibrium (LD) structure of the SNPs which are used in each gene-based association test.

```R
remotes::install_github('noahlorinczcomi/exset')
```

We simulate the following data which is necessary to perform a gene-based association test, then use the ```exset``` function to test the null hypothesis of no association between any SNPs in the gene-specific set and the phenotype.
```R
library(exset)
# a function to generate simulated data under H0 of the gene-based null hypothesis
simdata=function(niter,R) {
  m=nrow(R)
  Z=mvnfast::rmvn(niter,rep(0,m),R)
  S=rowSums(Z^2)
  S
}
m=50 # number of SNPs tested for this gene
R=0.5^toeplitz(0:(m-1)) # first-order autoregressive structure of the LD matrix with correlation parameter 0.5
S=simdata(niter=10,R=R) # generate 10 simulated gene-based test statistics under the gene-based null hypothesis
lam=eigen(R)$values # eigenvalues of the LD matrix
exset(S,lam,type1_error=0.05) # test results for each of the 10 genes

 [1] "fail to reject H0" "fail to reject H0"
 [3] "fail to reject H0" "fail to reject H0"
 [5] "fail to reject H0" "fail to reject H0"
 [7] "reject H0"         "fail to reject H0"
 [9] "fail to reject H0" "fail to reject H0"
```
The ```exset``` function circumvents having to specify the null distribution of the gene-based association test statistic and so only returns a decision as to whether the null hypothesis should be rejected at level ```type1_error``` or not, i.e. it does not return a P-value because it is not possible to exactly calculate.
