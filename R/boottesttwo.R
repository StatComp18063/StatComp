#' @title Two-sample bootstrap test
#' @description This function is used for the two-sample bootstrap test,
#' where the test statistic is the difference between the two sample means
#' @param x vector containing first sample.
#' @param y vector containing second sample.
#' @param b number of bootstrap replications
#' @return origtest: value of test statistic on original samples
#' @return pvalue: bootstrap p-value
#' @return teststatall: vector of bootstrap test statistics
#' @examples
#' x=rnorm(15,2,4)
#' y=rnorm(10,4,4)
#' boottesttwo(x,y,50)
#' @export
boottesttwo=function(x,y,b){
  #
  # x vector containing first sample.
  # y vector containing second sample.
  # b number of bootstrap replications
  # origtest: value of test statistic on original samples
  # pvalue: bootstrap p-value
  # teststatall: vector of bootstrap test statistics

  n1=length(x)
  n2=length(y)
  v=mean(y)-mean(x)
  z=c(x,y)
  counter=0
  teststatall=rep(0,b)
  for(i in 1:b){
    xstar=sample(z,n1,replace=T)
    ystar=sample(z,n2,replace=T)
    vstar=mean(ystar)-mean(xstar)
    if(vstar>=v) {counter=counter+1}
    teststatall[i]=vstar
  }
  pvalue=counter/b
  list(origtest=v,pvalue=pvalue,teststatall=teststatall)

}
