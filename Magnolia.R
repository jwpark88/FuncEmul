########## ERGM exmpales  ##########
## change statistics of gwdegree: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2031865/#FD24 ##
## ergm basics: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2743438/

rm(list=ls())
library(FuncEmul)
library(coda)
library(Bergm)
library(Rcpp)
library(RcppArmadillo)
library(DiceKriging)
library(DiceDesign)
library(MASS)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


########## Call functions ##########
source("http://www.stat.psu.edu/~mharan/batchmeans.R")


########## Run step ##########
data(package='ergm')

###  faus.magnolia network (1461 by 1461)  ###

data(faux.magnolia.high)
faux.magnolia.high
plot(faux.magnolia.high)

# summary statistics for data X
X = as.matrix(faux.magnolia.high)

# Parameter estimation via MLE or MPLE #
formula <- faux.magnolia.high ~ edges + gwesp(0.25,fixed=TRUE)
stat = summary(formula)
stat
Summary(X)

m <-ergm(formula,estimate="MPLE")
hat <- m$coef
summary(m)
COV <- solve(-m$hessian)

# system.time(  fit <- ergm(faux.magnolia.high ~ edges+gwesp(0.25,fixed=T), control = control.ergm(seed=1),verbose=T) )
# mcmc.diagnostics(fit, center=F)
# fit
# summary(fit)


### 1. Make Design points ###

# (1) from MPLE information 
p <- 2                                          # dimension of parameter 
summary(m)
m$coef[1] - 12*0.03813
m$coef[1] + 12*0.03813
m$coef[2] - 12*0.02861 
m$coef[2] + 12*0.02861 


Domain <- matrix(c(m$coef[1] - 12*0.03813, m$coef[1] + 12*0.03813,
                   m$coef[2] - 12*0.02861, m$coef[2] + 12*0.02861),2,p) # corresponding domain
num.point <- 3000                                                        # suppose we have 'num.point' particles
th <- matrix(0,num.point,p)

# Design for model 
ptm <- proc.time()
Design <- lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
proc.time()- ptm
th <- Design$design
th[,1] <- (Domain[2,1]-Domain[1,1])*th[,1]+Domain[1,1]
th[,2] <- (Domain[2,2]-Domain[1,2])*th[,2]+Domain[1,2]


### (2) Conduct ABC approach for generating particles

N <- 1      # number of importance sampling estimate
cycle <- 1  # numberof inner sampler
ptm <- proc.time()
summat = pAuxSamp(X, cycle, th, N, 20)  
ABCtime = proc.time()-ptm 
ABCtime 

dist = sqrt(  apply(  ( matrix( rep(stat,num.point), num.point, byrow = T) - summat[1:num.point,,1] )^2, 1, sum ) )
hist(dist)
eps = 50

m =  apply(  th[which(dist < eps),], 2, mean) 
S = cov( th[which(dist < eps),] )
thABC = mvrnorm( 400, m, S)


Domain <- rbind(apply(thABC,2,min),apply(thABC,2,max))                                        # Domain from ABC
# num.point <- 400                                                        # suppose we have 'num.point' particles
num.point <- 400
th <- matrix(0,num.point,p)

# Design for model 
ptm <- proc.time()
Design <- lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
proc.time()- ptm
th <- Design$design
th[,1] <- (Domain[2,1]-Domain[1,1])*th[,1]+Domain[1,1]
th[,2] <- (Domain[2,2]-Domain[1,2])*th[,2]+Domain[1,2]



### 2. Generate set of auxiliary variables for each particle point ###


N <- 1000      # number of importance sampling estimate
cycle <- 1     # numberof inner sampler
hat <- apply(thABC,2,mean)


ptm <- proc.time()
Sample = pResponseErgm(X, cycle, hat, N, 20)
IStime = proc.time()-ptm 
IStime 

y <- c()
for(i in 1:num.point){  
    cal = Sample%*%(th[i,]-hat) 
    mx = max(cal)
    y[i] = mx + log( mean( exp(cal-mx) ) )
    }

# unnormalized likelihood 
lhX <- c()
for(i in 1:num.point){ lhX[i] = stat%*%th[i,] }

# full likelihood 
y <- lhX - y

### 3. Final run of MCMC by using krigging estimate ###
Niter <- 25000
theta <- matrix(hat,1)
l.h.X <- as.numeric(hat%*%stat)

# calculating initial value of normalize const #
ptm <- proc.time()
m <- km(~ ., design = th[1:num.point,], response = matrix(y[1:num.point],num.point,1),covtype = "matern3_2")
GPtime = proc.time()-ptm
GPtime

x.point <- data.frame( t(theta[1,]) )

# Kriging
pred.m <- predict(m,x.point, "UK")
lhXZ <- ypred <- pred.m$mean

# starting cov matrix for proposal
COV <- cov(thABC)

# GLS bets coefficients 
beta.hat.gls <- c(coef(m)$trend1,coef(m)$trend2,coef(m)$trend3)

# MCMC run 
ptm <- proc.time()
parameter = GPmcmcErgm(Niter,theta,COV, lhXZ,beta.hat.gls,c(coef(m)$range,coef(m)$sd2),th,y,stat)
proc.time()-ptm



save.image("Magnolia.RData")




