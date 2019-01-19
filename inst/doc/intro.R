## ------------------------------------------------------------------------
library(lattice)

## ------------------------------------------------------------------------
n=seq(5,45,5)
x=rnorm(sum(n))
y=factor(rep(n,n),labels = paste('n=',n))
densityplot(~x|y,
            panel = function(x, ...) {
              panel.densityplot(x, col="DarkOliveGreen", ...)
              panel.mathdensity(dmath=dnorm,
                                args=list(mean=mean(x), sd=sd(x)),
                                col="darkblue")
              })



## ------------------------------------------------------------------------

data(quakes)
mini = min(quakes$depth)
maxi = max(quakes$depth)
int = ceiling((maxi - mini)/9)
inf = seq(mini, maxi, int)
quakes$depth.cat = factor(floor(((quakes$depth - mini) / int)),
labels=paste(inf, inf + int, sep="-"))
xyplot(lat ~ long | depth.cat, data = quakes)

## ------------------------------------------------------------------------
data(InsectSprays)
aov.spray = aov(sqrt(count) ~ spray, data = InsectSprays)
aov.spray

## ------------------------------------------------------------------------

termplot(aov.spray, se=TRUE, partial.resid=TRUE, rug=TRUE)

## ------------------------------------------------------------------------
search()

## ------------------------------------------------------------------------
library(grid)

## ------------------------------------------------------------------------
library(help = grid)

## ------------------------------------------------------------------------
set.seed(977)   ##  set a random seed
x=c(0:4)        ## generate the x
p=c(0.1,0.2,0.2,0.2,0.3)  ## generate the probability
cp=cumsum(p)   ## give the cumulate probability
n=1000      ## set the sample of size 1000 
r=numeric(n)
r=x[findInterval(runif(n),cp)+1] ## the inverse transform method
table(r) ## construct a relative frequency table 
ct=as.vector(table(r))
ct/sum(ct)  ## compare the empirical probability
ct/sum(ct)/p  ## the ratio of emparical probability and theorical probability


## ----fig.height=4--------------------------------------------------------
Y=p*n
counts = rbind(table(r),Y) ## combine r and Y into a data.frame
rownames(counts)=c('r','Y')  ## gives the row name of the data.frame
counts
barplot(counts/1000, main="multinomial distribution generation",
        xlab="x",ylab='prob', col=c("darkblue","red"),
        beside=TRUE,legend.text = c('empirical','theoretical'),
        args.legend = list(x = "topleft")) 

## ------------------------------------------------------------------------
set.seed(53) ## set a random seed
n=1000
k=0  ## counter for accepted
j=0  ## iterations
y=numeric(n)

 while (k<n) {
  u=runif(1)
  j=j+1
  x=runif(1)  ## #random variate from g
  if(x^2*(1-x)>u) {
    k=k+1
    y[k]=x
  }
}
j
hist(y,breaks=15,freq=F, main='empirical vs theoretical') ## plot the histogram
x1=seq(0,1,0.02)
y1=12*x1^2*(1-x1)
lines(x1,y1)   ## add the theoretical curve


## ------------------------------------------------------------------------
#generate a Exponential-Gamma mixture
set.seed(473)
n=1000
r=4
beta=2
lambda = rgamma(n, r, beta) ## generate the $\lambda$ from the Gamma(4,2) distribution

#now supply the sample of lambda??s as the parametric of Exponential distribution

y_lambda = rexp(n, lambda) ##  generate the mixture

hist(y_lambda,col='gray',ylim=c(0,800),,xlim=c(0,12),main='Exponential-Gamma mixture')


## ------------------------------------------------------------------------
set.seed(432)
m=10000
u=runif(m)
##  the function to generate cdf of Beta(3,3)
CDF=function(x){
  theta_hat=mean(x^3*u^2*(1-x*u)^2)
  CDF_hat=30*theta_hat
  #print(CDF_hat)
}

q=seq(0.1,0.9,0.1) 
Phi=pbeta(q,3,3)   ## q-quantile of Beta(3,3) by function 'pbeta'


##  q-quantile of Beta(3,3) by function CDF
cdf=rep(0,length(q))  
for (i in 1 : length(q)) { 
  cdf[i] = CDF(q[i]) 
  }
## Compare the estimator 
print(round(rbind(q,cdf,Phi),3))

## ------------------------------------------------------------------------
## the function of  generate samples from a Rayleigh distribution

  mc=function(x,R=10000,antithetic = TRUE){
    u=runif(R/2) 
    if (!antithetic)  v=runif(R/2) else
      v=1-u
    u=c(u,v)
    me=numeric(length(x))
    for (i in 1:length(x)){
      g=x[i]*(sqrt(2*log(1/(1-u))))
      me[i]=mean(g)
    }
    me
  }
sigma=c(1,2,3,4,5)  ## given the $\sigma$
mc1=mc(sigma)  ## antithetic variables method
mc2=mc(sigma,anti=F) ## inverse method.
print(rbind(sigma,mc1,mc2))  

### compare  the percent reduction in variance
  m=1000  ## repeat 1000 times 
  mc1=mc2=numeric(m)
  sigma=1  ## let $\sigma=1$
  for(i in 1:m){
    mc1[i]=mc(sigma,R=1000)
    mc2[i]=mc(sigma,R=1000,anti=F)
  }
print(sd(mc1))  ## the variance of antithetic variables method
print(sd(mc2))  ## the variance of inverse method
print((var(mc2)-var(mc1))/var(mc2)) ## the percent reduction in variance



## ------------------------------------------------------------------------
set.seed(345)
m=10000
u=runif(m)
x=sqrt(-2*log(1-(1-exp(-0.5))*u)) ## inverse method generte sample for inportance function
theta_hat=1/2-(1-exp(-1/2))/sqrt(2*pi)*mean(x)
theta_hat


## ------------------------------------------------------------------------
set.seed(436)
### geerate the random numbers, divide to 1000 groups, each group contains 5000 random numbers.

data=matrix(rlnorm(1000*2000,mean=0,sd=1),nrow=2000,ncol=1000)

logn_mean=rep(0,ncol(data)) ## definate the estimator of mean for 1000 times 
logn_median=rep(0,ncol(data)) ## defnate the estimator of median for 1000 times 
logn_G=rep(0,ncol(data)) ## definate the estimator of G for 1000 times

## repaet 1000 times

for (i in 1:ncol(data)){
  logn_mean[i]=mean(data[,i]) ## compute the sample mean
  data[,i]=sort(data[,i],decreasing=F) ## sort the random number from minimum to maximum.
  logn_median[i]=(data[,i][ncol(data)/2]+data[,i][ncol(data)/2+1])/2 ## compute the sample median
  
  ##compute the Gini ratio via order statistics
  n=nrow(data)
  d=rep(0,n)
  for(j in 1:n) {
    d[j]=(2*j-n-1)*data[,i][j]
  }
logn_G[i]=mean(d)/(n^2*logn_mean[i])
}
hist(logn_G)  
hist(logn_mean)
hist(logn_median)


## ------------------------------------------------------------------------
set.seed(1436)
### geerate the random numbers, divide to 1000 groups, each group contains 5000 random numbers.

data=matrix(runif(1000*2000),nrow=2000,ncol=1000)

unif_mean=rep(0,ncol(data)) ## definate the estimator of mean for 1000 times 
unif_median=rep(0,ncol(data)) ## defnate the estimator of median for 1000 times 
unif_G=rep(0,ncol(data)) ## definate the estimator of G for 1000 times

## repaet 1000 times

for (i in 1:ncol(data)){
  unif_mean[i]=mean(data[,i]) ## compute the sample mean
  data[,i]=sort(data[,i],decreasing=F) ## sort the random number from minimum to maximum.
  unif_median[i]=(data[,i][ncol(data)/2]+data[,i][ncol(data)/2+1])/2 ## compute the sample median
  
  ##compute the Gini ratio via order statistics
  n=nrow(data)
  d=rep(0,n)
  for(j in 1:n) {
    d[j]=(2*j-n-1)*data[,i][j]
  }
unif_G[i]=mean(d)/(n^2*unif_mean[i])
}
hist(unif_G)  
hist(unif_mean)
hist(unif_median)


## ------------------------------------------------------------------------
## write a function  n denote sample size, m denote repaet times, mu and sigma denote the parameter  in lognorm-distribution

mc_fun=function(n,m,mu,sigma){
  set.seed(528)
  data=matrix(rlnorm(m*n,mu,sigma),nrow=n,ncol=m)
  logn_G=rep(0,m)
  logn_mean=rep(0,m)
  for(i in 1:m){
    logn_mean[i]=mean(data[,i])
     data[,i]=sort(data[,i],decreasing=F) ##sort the random number from minimum to maximum.(obtain the order statistics)
     
    ##compute the Gini ratio via order statistics]
     
       d=rep(0,n)
       for(j in 1:n)  d[j]=(2*j-n-1)*data[,i][j]
      
       logn_G[i]=mean(d)/(n^2*logn_mean[i])  
  }
  ##print the sample mean.sample sd,  and 95% CI??
print(c(mean(logn_G),sd(logn_G),mean(logn_G)-2*sd(logn_G),
        mean(logn_G)+2*sd(logn_G)))
}
Gini=mc_fun(1000,1000,1,1)

## ------------------------------------------------------------------------


##Simulate from a Multivariate Normal Distribution 
library(MASS) 
set.seed(4356)
Sigma = matrix(c(1,0.5,0.5,1),2,2) 
Sigma 
x=mvrnorm(n=1000, rep(0, 2), Sigma)
cor.test(x[,1],x[,2],method='pearson')  ##  person correlation coefficient
cor.test(x[,1],x[,2],method='spearman') ## spearman rank correlation coefficient
cor.test(x[,1],x[,2],method='kendall')  ## kendall correlation coefficient



## ------------------------------------------------------------------------
set.seed(1)
library(bootstrap)
data=law  # data set
head(data)
n=nrow(data)  # sample size
theta.hat=cor(data[,1],data[,2]) # the sample correlation.
#### jackknife
theta.jack=numeric(n) 
for (i in 1:n) {   ##leave one out
  data.jack=data[-i,]
  theta.jack[i]=cor(data.jack[,1],data.jack[,2])
}
jack.bias=(n-1)*mean(theta.jack-theta.hat) # jackknife bias
jack.sd=sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2)) # jackknife sd
round(c(theta.hat=theta.hat,theta.jack=mean(theta.jack),
        jack.bias=jack.bias,jack.sd=jack.sd),5)

## ------------------------------------------------------------------------
library(boot)

set.seed(1)
#  define a function to compute the boostrap estimator for sample mean.
air.fun <- function(data,i) {
  m <- mean(data$hours[i]) # sample mean
  n <- length(i)
  v <- (n-1)*var(data$hours[i])/n^2 # sample variance
  c(m, v)
}

# use  boot function, repeat 1000 times
air.boot <- boot(aircondit, air.fun, R = 1000) 
air.boot #t1 denote sample mean bootstrap, t2 denote sample variance bootstrap

# use boot.ci function to obtain the CI
ci <- boot.ci(air.boot,type=c("norm","basic","perc","bca")) 
ci

## ------------------------------------------------------------------------
library(bootstrap) #loading packages
data=as.matrix(scor) # let data.frame become matrix
head(data)
n=nrow(data) # sample size
lambda.hat=eigen(t(data)%*%data/n)$values ## sample eigenvalues
theta.hat=lambda.hat[1]/sum(lambda.hat)  ## sample estimator 
#### jackknife
theta.jack=numeric(n)
for(i in 1:n){          ##leave one out
  data.jack=data[-i,]
  lambda.jack=eigen(t(data.jack)%*%data.jack/(n-1))$values
  theta.jack[i]=lambda.jack[1]/sum(lambda.jack)
}
jack.bias=(n-1)*mean(theta.jack-theta.hat) #jackknife  bias
jack.sd=sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2)) #jackknife se
round(c(theta.jack=mean(theta.jack),theta.hat=theta.hat,
        jack.bias=jack.bias,jack.sd=jack.sd),7)


## ------------------------------------------------------------------------
library(DAAG) # loading packages
data=ironslag #data set
n=nrow(data) # sample size
# creat a matrix to storage the index of test for leave-two-out CV.
index=matrix(nrow=n*(n-1)/2,ncol=2) 

#we obtain all the index via random sampling 2 observations, 
# there is 53*52/2=1378 outcomes.

a=c(1:n) # index of sample size, is 1:n
l=1
for(i in 1: (n-1)){
  for(j in (i+1):n){
index[l,]=(a[c(i,j)]) 
l=l+1
  }    
}
##
e1= e2 = e3 = e4 = matrix(nrow=nrow(index),ncol=2)

# fit models on leave-two-out samples
for(l in 1:nrow(index)){
  k=index[l,] #the index of testing data 
  y=data[-k,2] # the y of training data
  x=data[-k,1] # the x of training data
  
## linear model
  J1 = lm(y ~ x) # linear model
  yhat1 = J1$coef[1] + J1$coef[2] * data[k,1] # yhat on testing data
  e1[l,] = data[k,2] - yhat1
  
##quadratic model
  J2 <- lm(y ~ x + I(x^2)) # quadratic model 
yhat2 <- J2$coef[1] + J2$coef[2] * data[k,1] + # yhat on testing data
J2$coef[3] * data[k,1]^2
e2[l,] = data[k,2] - yhat2
  
##Exponential model
J3 = lm(log(y) ~ x) #exp model
logyhat3 =J3$coef[1]+J3$coef[2] * data[k,1] # yhat on testing data
yhat3 <- exp(logyhat3)
e3[l,] = data[k,2] - yhat3

## Log-Log model
  J4 <- lm(log(y) ~ log(x)) # log-log model
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(data[k,1]) # yhat on testing data
  yhat4 <- exp(logyhat4)
  e4[l,] = data[k,2] - yhat4

}
print(c(mse1=mean(e1^2), mse2=mean(e2^2), mse3=mean(e3^2), mse4=mean(e4^2)))

## ------------------------------------------------------------------------
cvm.d<-function(x,y){ #compute the Cramer-von Mises statistic
  n<-length(x);m<-length(y)
  Fn<-ecdf(x);Gm<-ecdf(y)
  W2<-((m*n)/(m+n)^2)*
      (sum((Fn(x)-Gm(x))^2)+sum((Fn(y)-Gm(y))^2))
  W2
}

## ------------------------------------------------------------------------
library(boot)
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

set.seed(611)
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
n<-length(x)
m<-length(y)
N <- 1:(m+n)
reps <- numeric(R) #storage for replicates
W2 <- cvm.d(x, y)
W2
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(N, size = n, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- cvm.d(x1, y1)
}
p <- mean(c(W2, reps) >= W2)
p


## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
f <- function(x, theta,eta) {
  stopifnot(theta > 0)
  return(1/(theta*pi*(1+((x-eta)/theta)^2))) 
} 
xt <- x[i-1] 
y <-  xt+runif(1,-1,1)
m <- 10000
theta <- 1
eta <- 0
x <- numeric(m) 
x[1] <- runif(1,-1,1) 
k <- 0
u <- runif(m)

for (i in 2:m) {
  xt <- x[i-1]
  y <- xt+runif(1,min=-1,max=1)
  num <- f(y, theta,eta)*dnorm(xt,mean=y,sd=1) 
  den <- f(xt, theta,eta)*dnorm(y,mean=xt,sd=1)
  if (u[i] <= num/den) x[i] <- y else {
    x[i] <- xt
    k <- k+1         #y is rejected 
  } 
}
print(k) 

## ---- echo=TRUE----------------------------------------------------------
index <- 1500:2000
y1 <- x[index] 
plot(index, y1, type="l", main="", ylab="x") 

## ---- echo=TRUE----------------------------------------------------------
gene<- quantile(x[-(1:1000)],seq(0.1,0.9,0.1))
empi<-qt(seq(0.1,0.9,0.1),1)
rbind(gene,empi)

## ---- echo=TRUE----------------------------------------------------------
b <- 1001     #discard the burnin sample 
y <- x[b:m] 
probs <- seq(0.1,0.9,0.01)
QR <- qt(probs,df=1)
Q <- quantile(x, probs)

qqplot(QR, Q, main="", 
       xlab="Cauchy Quantiles", ylab="Sample Quantiles")

hist(y, breaks="scott", main="", xlab="", freq=FALSE) 
lines(QR, f(QR, 1,0))

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)

w <- .25   #width of the uniform support set 
m <- 5000  #length of the chain 
burn <- 1000  #burn-in time 
total <- 197 
x <- numeric(m)  #the chain

## ---- echo=TRUE----------------------------------------------------------
group <- c(125,20,18,34)

prob <- function(y, group) { 
  # computes (without the constant) the target density 
  if (y < 0 || y >= 1) 
    return (0) 
  return(((2+y)/4)^group[1] * 
           ((1-y)/4)^group[2] * ((1-y)/4)^group[3] *
           ((y/4)^group[4]))
} 

u <- runif(m)          #for accept/reject step 
v <- runif(m, -w, w)   #proposal distribution
x[1] <- 0.5
for (i in 2:m) {
  y <- x[i-1] + v[i]
  if (u[i] <= prob(y, group) / prob(x[i-1], group)) 
    x[i] <- y else 
      x[i] <- x[i-1]
}
print(group) 
print(round(group/total, 3))

## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
 xb <- x[(burn+1):m] 
 print(mean(xb))    

## ------------------------------------------------------------------------
lden=function(theta)125*log(2+theta)+38*log(1-theta)+34*log(theta)
MCMC=function(x,burn_in=0,N=burn_in+5000,print_acc=F){
  y=1:N#for MCMC generation
  ld=lden#use log-likelihood, things could be a little bit different.
  acc=0
  for(i in 1:N){
    p=runif(1)
    y[i]=x=if(runif(1)<exp(ld(p)-ld(x))){acc=acc+1;p}else x
  }
  if(print_acc)print(acc/N)
  y[(burn_in+1):N]
}
plot(MCMC(0.5,print_acc=T))#accept rate~0.15

## ------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
 # psi[i,j] is the statistic psi(X[i,1:j])
 # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)

   psi.means <- rowMeans(psi)     #row means
   B <- n * var(psi.means)        #between variance est.
   psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + B/n+(B/(n*k))     #upper variance est.
   r.hat <- v.hat / W             #G-R statistic
return(r.hat)
  }

## ------------------------------------------------------------------------
M1=cumsum((MC1=MCMC(0.2)))/1:5000#mean_cumsum
M2=cumsum((MC2=MCMC(0.4)))/1:5000#mean_cumsum
M3=cumsum((MC3=MCMC(0.6)))/1:5000#mean_cumsum
M4=cumsum((MC4=MCMC(0.8)))/1:5000#mean_cumsum
psi=rbind(M1,M2,M3,M4)
plot((R=sapply(1:100,function(i)Gelman.Rubin(psi[,1:i]))),main="R value of Gelman-Rubin method",type='l',ylab = "")

res=(1:100)[R[2:100]>1.2][sum(R[2:100]>1.2)]+1
c(mean(c(MC1[res:5000],MC2[res:5000],MC3[res:5000],MC4[res:5000])),
  var(c(MC1[res:5000],MC2[res:5000],MC3[res:5000],MC4[res:5000])))


## ------------------------------------------------------------------------
n=25
A=rep(0,length(c(4:n)))
for(k in 4:n){
f=function(x){
pt(sqrt(x^2*k/(k+1-x^2)),k)-pt(sqrt(x^2*(k-1)/(k-x^2)),k-1)} 
A[k-3]=uniroot(f,c(0.00001,sqrt(k)-0.00001))$root #solve the intersection points
}
plot(A,type='l')

## ------------------------------------------------------------------------
n=100
A=rep(0,length(c(4:n)))
for(k in 4:n){
f=function(x){
pt(sqrt(x^2*k/(k+1-x^2)),k)-pt(sqrt(x^2*(k-1)/(k-x^2)),k-1)}
A[k-3]=uniroot(f,c(0.00001,sqrt(k)-0.00001))$root #solve the intersection points
}
plot(A,type='l')

## ------------------------------------------------------------------------
n=500
A=rep(0,length(c(4:n)))
for(k in 4:n){
f=function(x){
pt(sqrt(x^2*k/(k+1-x^2)),k)-pt(sqrt(x^2*(k-1)/(k-x^2)),k-1)}
A[k-3]=uniroot(f,c(0.00001,sqrt(k)-0.00001))$root #solve the intersection points
}
plot(A,type='l')

## ------------------------------------------------------------------------
n=1000
A=rep(0,length(c(4:n)))
for(k in 4:n){
f=function(x){
pt(sqrt(x^2*k/(k+1-x^2)),k)-pt(sqrt(x^2*(k-1)/(k-x^2)),k-1)}
A[k-3]=uniroot(f,c(0.00001,sqrt(k)-0.00001))$root #solve the intersection points
}
plot(A,type='l') 

## ------------------------------------------------------------------------
## define the pdf of cauchy distribution
pdf_cauchy=function(t,theta,eta){
  1/(theta*pi*(1+((t-eta)/theta)^2))
}
x=seq(0,10,0.1)
num_int=rep(0,length(x))
p=rep(0,length(x))
# let theta=1,eta=0.
for(i in 1:length(x)){
num_int[i]=integrate(pdf_cauchy, -Inf, x[i],theta=1,eta=0)$value
p[i]=pcauchy(x[i],location = 0, scale = 1) # by pcauchy
}
plot(x,num_int) # plot the num_int
lines(x,p,col='green') # plot the pcauchy

# let theta=1,eta=1.
for(i in 1:length(x)){
num_int[i]=integrate(pdf_cauchy, -Inf, x[i],theta=1,eta=1)$value
p[i]=pcauchy(x[i],location = 1, scale = 1) # by pcauchy
}
plot(x,num_int) # plot the num_int
lines(x,p,col='green') # plot the pcauchy

## ------------------------------------------------------------------------
## use loops to fit linear model.
formulas = list(
mpg ~ disp,
mpg ~ I(1 / disp), mpg~disp+wt, mpg~I(1/disp)+wt
)
out1=list()
for(i in 1:4){
  out1[[i]]=lm(formulas[[i]],data=mtcars)
}
out1

### use lapply() function to fit linear model
out2=lapply(formulas,lm,data=mtcars)
out2

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
         rows <- sample(1:nrow(mtcars), rep = TRUE)
         mtcars[rows, ]
})
### loops
out3=list()
for(i in 1:10){
  out3[[i]]=lm(mpg ~ disp,data=bootstraps[[i]])
}
out3

### lapply
newdata=list()
for(i in 1:10){
  
  newdata[[i]]=bootstraps[[i]][,c(1,3)]
}
out4=lapply(newdata,lm)
out4


## ------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
## compute R2 for out1 
R2_1=lapply(out1,rsq)
R2_1
## compute R2 for out2
R2_2=lapply(out2,rsq)
R2_2
## compute R2 for out3
R2_3=lapply(out3,rsq)
R2_3
## compute R2 for out4
R2_4=lapply(out4,rsq)
R2_4

## ------------------------------------------------------------------------
trials <- replicate(
         100,
         t.test(rpois(10, 10), rpois(7, 10)),
         simplify = FALSE
       )
## define a function to extract p.value
p_fun=function(f)(f$p.value)
sapply(trials,p_fun)

## ------------------------------------------------------------------------
chisq.test1=function (x, y = NULL, correct = TRUE, p = rep(1/length(x), length(x)), 
    rescale.p = FALSE, simulate.p.value = FALSE, B = 2000) 
{
    DNAME <- deparse(substitute(x))
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        if (min(dim(x)) == 1L) 
            x <- as.vector(x)
    }
    if (!is.matrix(x) && !is.null(y)) {
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        DNAME2 <- deparse(substitute(y))
        xname <- if (length(DNAME) > 1L || nchar(DNAME, "w") > 
            30) 
            ""
        else DNAME
        yname <- if (length(DNAME2) > 1L || nchar(DNAME2, "w") > 
            30) 
            ""
        else DNAME2
        OK <- complete.cases(x, y)
        x <- factor(x[OK])
        y <- factor(y[OK])
        if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) 
            stop("'x' and 'y' must have at least 2 levels")
        x <- table(x, y)
        names(dimnames(x)) <- c(xname, yname)
        DNAME <- paste(paste(DNAME, collapse = "\n"), "and", 
            paste(DNAME2, collapse = "\n"))
    }
    if (any(x < 0) || anyNA(x)) 
        stop("all entries of 'x' must be nonnegative and finite")
    if ((n <- sum(x)) == 0) 
        stop("at least one entry of 'x' must be positive")
    if (simulate.p.value) {
        setMETH <- function() METHOD <<- paste(METHOD, "with simulated p-value\n\t (based on", 
            B, "replicates)")
        almost.1 <- 1 - 64 * .Machine$double.eps
    }
    if (is.matrix(x)) {
        METHOD <- "Pearson's Chi-squared test"
        nr <- as.integer(nrow(x))
        nc <- as.integer(ncol(x))
        if (is.na(nr) || is.na(nc) || is.na(nr * nc)) 
            stop("invalid nrow(x) or ncol(x)", domain = NA)
        sr <- rowSums(x)
        sc <- colSums(x)
        E <- outer(sr, sc, "*")/n
        v <- function(r, c, n) c * r * (n - r) * (n - c)/n^3
        V <- outer(sr, sc, v, n)
        dimnames(E) <- dimnames(x)
        if (simulate.p.value && all(sr > 0) && all(sc > 0)) {
            setMETH()
            tmp <- .Call(C_chisq_sim, sr, sc, B, E)
            STATISTIC <- sum(sort((x - E)^2/E, decreasing = TRUE))
            PARAMETER <- NA
            PVAL <- (1 + sum(tmp >= almost.1 * STATISTIC))/(B + 
                1)
        }
        else {
            if (simulate.p.value) 
                warning("cannot compute simulated p-value with zero marginals")
            if (correct && nrow(x) == 2L && ncol(x) == 2L) {
                YATES <- min(0.5, abs(x - E))
                if (YATES > 0) 
                  METHOD <- paste(METHOD, "with Yates' continuity correction")
            }
            else YATES <- 0
            STATISTIC <- sum((abs(x - E) - YATES)^2/E)
            PARAMETER <- (nr - 1L) * (nc - 1L)
            PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
        }
    }
    else {
        if (length(dim(x)) > 2L) 
            stop("invalid 'x'")
        if (length(x) == 1L) 
            stop("'x' must at least have 2 elements")
        if (length(x) != length(p)) 
            stop("'x' and 'p' must have the same number of elements")
        if (any(p < 0)) 
            stop("probabilities must be non-negative.")
        if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
            if (rescale.p) 
                p <- p/sum(p)
            else stop("probabilities must sum to 1.")
        }
        METHOD <- "Chi-squared test for given probabilities"
        E <- n * p
        V <- n * p * (1 - p)
        STATISTIC <- sum((x - E)^2/E)
        names(E) <- names(x)
        if (simulate.p.value) {
            setMETH()
            nx <- length(x)
            sm <- matrix(sample.int(nx, B * n, TRUE, prob = p), 
                nrow = n)
            ss <- apply(sm, 2L, function(x, E, k) {
                sum((table(factor(x, levels = 1L:k)) - E)^2/E)
            }, E = E, k = nx)
            PARAMETER <- NA
            PVAL <- (1 + sum(ss >= almost.1 * STATISTIC))/(B + 
                1)
        }
        else {
            PARAMETER <- length(x) - 1
            PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
        }
    }
    names(STATISTIC) <- "X-squared"
    names(PARAMETER) <- "df"
    if (any(E < 5) && is.finite(PARAMETER)) 
        warning("Chi-squared approximation may be incorrect")
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name = DNAME, observed = x, 
        expected = E, residuals = (x - E)/sqrt(E), stdres = (x - 
            E)/sqrt(V)), class = "htest")
}

## ------------------------------------------------------------------------
chisq.test2=function (x, y)
{
  if (length(x) != length(y)) 
  stop("'x' and 'y' must have the same length")
  OK <- complete.cases(x, y)
        x <- factor(x[OK])
        y <- factor(y[OK])
        if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) 
            stop("'x' and 'y' must have at least 2 levels")
        x <- table(x, y)
        n <- sum(x)
        p = rep(1/length(x), length(x))
        p <- p/sum(p)
  
        
        sr <- rowSums(x)
        sc <- colSums(x)
        E <- outer(sr, sc, "*")/n
        v <- function(r, c, n) c * r * (n - r) * (n - c)/n^3
        V <- outer(sr, sc, v, n)
        dimnames(E) <- dimnames(x)
        STATISTIC <- sum(sort((x - E)^2/E, decreasing = TRUE))
         STATISTIC 
      
}

options(warn=-1)
x <- c(89,37,30,28,2,14,25)
y <- c(29,47,10,8,14,32,52)

chisq.test1(x,y)$statistic
chisq.test2(x,y)


library(microbenchmark)
microbenchmark(
  chisq.test1(x,y)$statistic,
  chisq.test2(x,y)
)


