#' @title Simulation Poisson Process
#' @description A Poisson Process sampler using R
#' @importFrom stats runif
#' @param n the number of samples
#' @param lambda the mean of the Poisson distribution
#' @return a random sample of size \code{n} for Poisson Process
#' @examples

#' temp=poisrand(1000,5)
#' head(temp)
#' mean(temp)
#' @export
poisrand=function(n,lambda){
  #
  # n is the number of simulations
  # lambda is the mean of the Poisson distribution.
  #
  poisrand=rep(0,n)
  for(i in 1:n){
    xt=0
    t=0
    while(t<1){
      x=xt
      y=-(1/lambda)*log(1-runif(1))
      t=t+y
      xt=xt+1
    }
    poisrand[i]=x
  }
  poisrand
}
