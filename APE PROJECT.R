install.packages("AdequacyModel")
library(AdequacyModel)

DATA=c(1,4,4,7,11,13,15,15,17,18,19,19,20,20,22,23,28,29,31,32,36,37,47,48,49,50,54,54,55,59,59,61,61
       ,66,72,72,75,78,78,81,93,96,99,108,113,114,120,120,120,123,124,129,131,137,145,151,156,171,
       176,182,188,189,195,203,208,215,217,217,217,224,228,233,255,271,275,275,275,286,291,312,
       312,312,315,326,326,329,330,336,338,345,348,354,361,364,369,378,390,457,467,498,517,566,
       644,745,871,1312,1357,1613,1630)
y=DATA

#FIGURE 3
hist(DATA,main="Relative Histogram of the Coal Mine Dataset.",xlab="Days",ylab="Relative Frequency",xlim = c(0,1800),prob=TRUE)

#FIGURE 4( Scaled TTT plot of the coal mine dataset.) 
TTT(DATA,lwd=2,lty=2,col="black",grid="True")

# Minimum, maximum, first quartile, median, and third quartile of the dataset
summary(DATA)

#FOR APE DISTRIBUTION

#FOR THE MLE AND LOG_LIKELIHOOD OF APE

#Set the MLE function for alpha
n=109
LOG_ALPHA_HAT <- function(y, lambda) {
  n   <- length(y);
  NUM <- ((1/n)*sum(y)) - 1/lambda;
  DEN <- (1/n)*sum(y*exp(-lambda*y));
  NUM/DEN; }

#Set the log-likelihood function
LOGLIKE <- function(y, lambda) {
  la <- LOG_ALPHA_HAT(y, lambda);
  if (la == 0) {
    LL <- n*log(lambda) - lambda*sum(y); } else {
      LL <- n*la + n*log(la/expm1(la)) + n*log(lambda) - 
        lambda*sum(y) - la*sum(exp(-lambda*y)); }
  LL; }


#Maximise the log-likelihood function
OBJECTIVE  <- function(lambda) { - LOGLIKE(y , lambda) }
START      <- c(1/mean(DATA))
NLM        <- nlm(OBJECTIVE, p = START);
LLMAX      <- - NLM$minimum;
MLE_LAMBDA <- NLM$estimate;
MLE_ALPHA  <- exp(LOG_ALPHA_HAT(y, MLE_LAMBDA));
MLE        <- data.frame(alpha = MLE_ALPHA, lambda = MLE_LAMBDA, loglike = LLMAX);
rownames(MLE) <- 'MLE'
MLE


#FOR THE K-S STATISTIC AND P-VALUE OF APE

#CDF OF APE

alpha_ape=0.00366583
lambda_ape=0.0009550325
cdf <- function(y, alpha,lambda){ 
  
  if(alpha!=1){
    
    apecdf<-((alpha^(1-exp(-lambda*y)))-1)/ (alpha-1)}else if(alpha==1){
      
      apecdf<- 1-(exp(-lambda*y))}
  
  return(apecdf)
  
}

t1 <- ks.test(x = y,y = "cdf", alpha_ape, lambda_ape)
t1

#For the use in table 2
ks_ape = 0.061673
p_ape = 0.8014


#FOR WEIBULL DISTRIBUTION
#MLE and LOG_LIKELIHOOD for Weibull
library(weibullness)
weibull.mle(DATA, threshold=0 )
loglike<- (n*log(0.8848074)-n*0.8848074*log( 218.6772) +sum((0.8848074 -1)*log(y))-sum((y/ 218.6772)^(0.8848074)))
loglike

#k-s test statistic and p-Value for Weibull
alpha_weibull=0.8848074
lambda_weibull=218.6772
weibull <- ks.test(DATA, "pweibull", shape =alpha_weibull, scale = lambda_weibull)
weibull

#For the use in table 2
ks_weibull = 0.078444
p_weibull =   0.5136 


#FOR GAMMA DISTRIBUTION
#MLE and LOG_LIKELIHOOD for GAMMA
install.packages("EnvStats")
library(EnvStats)
est.par=egamma(DATA, method = "mle", ci = FALSE, 
       ci.type = "two-sided", ci.method = "normal.approx", 
       normal.approx.transform = "kulkarni.powar", conf.level = 0.95)
est.par
LLab <- sum(dgamma(y,shape=0.8559939,scale=272.5733252,log=TRUE))
LLab

#k-s test statistic and p-Value for GAMMA
alpha_gamma=0.8559939
lambda_gamma=272.5733252
gamma <- ks.test(DATA, "pgamma", shape = alpha_gamma , scale =lambda_gamma )
gamma


#For the use in table 2
ks_gamma = 0.082138
p_gamma = 0.4539


#TABLE 2

table = data.frame(
  The_Model = c("Gamma","Weibull","APE"),
  data.frame(
    alpha_cap = c(alpha_gamma,alpha_weibull,alpha_ape),
    lambda_cap = c(lambda_gamma,lambda_weibull,lambda_ape)),
  Log_Likelihood = c(LLab,loglike,LLMAX),
  K_S_Statistic = c(ks_gamma,ks_weibull,ks_ape),
  p_value = c(p_gamma,p_weibull,p_ape)
  )

table


#Figure 5(a) and 5(b) using the results obtained above:


alpha=0.00366583
lambda=0.0009550325
#PDF OF APE DISTRIBUTION
pdf=if(alpha != 1){
  (log(alpha)/(alpha-1))*(lambda*exp(-lambda*y))*(alpha^(1-exp(-lambda*y)))}else{
    lambda*exp(-lambda*y)}
pdf

#Histogram of coal mine dataset
hist(DATA,main=" The Histogram and the fitted APE distribution",xlab="x",ylab="Density",xlim = c(0,1800),prob=TRUE)
#The fitted APE distribution
lines(DATA,pdf,col="blue","l")


#FIGURE 5(b)
#SURVIVAL FUNCTION OF APE DISTRIBUTION

# functions: 
cdf <- function(y, alpha,lambda) {
  if (alpha != 1) {
    apecdf <-((alpha^(1-exp(-lambda * y))) - 1)/ (alpha - 1)
  } else {
    apecdf<- 1-(exp(-lambda * y))
  }
  return(apecdf)
} 

S <- function(y, alpha, lambda) {
  if (alpha != 1) {
    (alpha / (alpha - 1)) * (1 - alpha^(-exp(-lambda * y)))
  } else {
exp(-lambda * y) 
  }
}

# theoretical cdf and survival function for given parameters: 
S_th <- S(y = y, alpha = alpha, lambda = lambda)
cdf_th <- cdf(y = y, alpha = alpha, lambda = lambda)

# implementation check, S(x) = 1 - F(x) where F(.) is the cdf: 
all.equal(S_th, 1 - cdf_th)
#> [1] TRUE

# empirical cdf and survival function: 
## argument x of is specified as a numeric vector of observations; ie. y in your case: 

cdf_emp_fun <- ecdf(x = y) #  is a function
cdf_emp <- cdf_emp_fun(y) # values

# S as the counterpart of the cdf: 
S_emp <- 1 - cdf_emp

# Plotting the survival functions: 
plot(stepfun(x = y, y = c(1, S_emp)), main = "The fitted APE survival function and empirical
     survival function for coal-mining data.", xlab = "x", ylab = "Survival Function") # empirical S
lines(x = y, y = S_th, col = "red", type = "l") # theoretical S
legend("topright",legend=(c(p1="Empirical Survival Function",p2="Fitted Survival Function")),col=c("black","red",15),lty=1)

