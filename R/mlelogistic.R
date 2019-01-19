#' @title Computing the MLE for logistic distribution
#' @description Computing the MLE for logistic distribution using R,
#' location parameter is theta, scale parameter is 1.
#' @param x  vector, the sample
#' @param theta0 the mean of x
#' @param numstep maximum repeat times
#' @param eps The iteration error
#' @return theta1 is the MLE
#' @return check is the iteration error
#' @return realnumsteps is the actual number of iterations
#' @examples
#' x=rlogis(100,2,1)
#' mlelogistic(x)
#' @export
mlelogistic=function(x,theta0=mean(x),numstep=100,eps=.0001){
  n=length(x)
  numfin=numstep
  small=1.0*10^(-8)
  ic=0
  istop=0
  while(istop==0){
    ic=ic+1
    expx=exp(-(x-theta0))
    lprime=n-2*sum(expx/(1+expx))
    ldprime=-2*sum(expx/(1+expx)^2)
    theta1=theta0-(lprime/ldprime)
    check=abs(theta0-theta1)/abs(theta0+small)
    if(check<eps){istop=1}
    theta0=theta1
  }
  list(theta1=theta1,check=check,realnumsteps=ic)
}
