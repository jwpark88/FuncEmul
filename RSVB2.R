########## Conduct GP Interpolation (Attraction repulsion) ##########

rm(list=ls())
library(FuncEmul)
library(sp)
library(gstat)
library(fields)
library(classInt)
library(spatstat)
library(MASS)
library(DiceKriging)
library(DiceDesign)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
library(Rcpp)
library(RcppArmadillo)

###   Call C++ function   ###
source("http://www.stat.psu.edu/~mharan/batchmeans.R")


# BDmcmc(vec center, double range, mat initial, mat initialdist, vec parameter, int n)
# DMH(vec center, double range, mat initial, mat initialdist, double lhX, int Niter, mat theta, mat COV, int n)
# pResponse(vec center, double range, mat initial, mat initialdist, vec hatparameter, int inner, mat Designmat, int m, int num)
# GPmcmc(int Niter, mat theta, mat COV, double lhX, double logconst, vec betahat, vec phihat, mat Designmat, vec y, mat Distmat)



###                   run step                  ###
###      1. Load RSV-B2 datasets                ###
###    A-2A_16h_A1 is used                      ###

# Followings also can be used
# A-2A_16h_A1 (3077)
# A-2A_16h_A2 (3710)
# A-2A_16h_A3 (3286)
# A-2A_3h_A1  (4412)
# A-2A_3h_A2  (4551)
# A-2A_3h_A3  (4241)
# A-1B2A_16h_A1 (3949)


# Interaction function 
R <- 0 
rho <- function(x){
       result <- 0 + (x>100)*1 + ((x>r1) & (x<=100))*(1 + 1/( t3*(x-r2) )^2) + 
       ((x>R) & (x<= r1)) * (t1 - ( sqrt(t1)*(x-t2)/(t2-R) )^2)
       return(result)}  
data(RSVB2)
X <- as.matrix( RSVB2 )
r.X <- rdist(X)
dim(X)
plot(X)

# draw pair correlation function as EDA
center <- c(2000,2000); radius <- (1350)
W <- disc(radius, center)
pppp <- pcf(ppp(X[,1],X[,2],window=W))
plot(pppp,xlim=c(1,70))


# calculating initial likelihood 
hat = c(3e-4,1.35,12,0.25)  
lambda = hat[1]
t1 <- hat[2]; t2 <- hat[3]; t3 <- hat[4]   # theta parameter (for interaction function)

# solution for r1 r2 solving two equations #
a <- (t2-R)^2 / (t1*t3^2)
b <- t3^2 * (t1 - 1)
d <- 27*a*b^2
c <- ( d + sqrt((d+2)^2-4) + 2 )^(1/3)
deltar <- 1/(3*b) * ( c/2^(1/3) + 2^(1/3)/c + 1 )

r1 <- a / deltar^1.5 + t2
r2 <- r1 - sqrt(deltar)

rho.X <- rho(r.X)                            
diag(rho.X) <- NA
rho.X.sum <- apply(log(rho.X),2,sum,na.rm=T)
l.h.X <- sum( pmin(rho.X.sum,1.2) ) + (dim(X)[1])*log(lambda)   # likelihood for the initial data


### 2. Make Design points ###
hat                                            # initial guess  
l.h.X                                            # likelihood for the initial data
n <- 30000                       # inner sampler
COV <- diag(c(3e-10,0.01,0.1,0.01))

#short run of dmh
ptm <- proc.time() 
dmh <- DMH(center, radius, X, r.X, l.h.X, 3000, matrix(hat,1), COV, n)
proc.time() - ptm


# calculating initial likelihood 
Liangshorttime = 68352*(3000/40000)
parameter = dmh[1000:3000,]
hat = apply(parameter[1:2000,],2,mean)

### 2. Make Design points ###

p <- 4                                          # dimension of parameter 
num.point <- 400                                # suppose we have 'num.point data'
th <- matrix(0,num.point,p)
for(i in 1:4){ th[,i] = unique(parameter[,i])[1:num.point] }
length(unique(th[,1]))

### 3. Calculate importance sampling estimate  (sampling and calculating) ###
N <- 2000 # number of importance sampling estimate
n <- 30000 # number of Birth Death MCMC run 
ptm <- proc.time()
Sample = pResponse(center,(radius),X,r.X,hat,n,th,N,20)
IStime = proc.time()-ptm 
IStime

Max = apply(Sample,2,max) 
maxx = matrix(rep(Max,N),N,byrow=T)
y = Max + log( apply( exp((Sample)-maxx), 2, mean ) )

# unnormalized likelihood 
lhX <- c()
for(i in 1:400){ lhX[i] <- evalunNormX(th[i,], r.X) }

# full likelihood
y = lhX - y 

### 5. Final run of MCMC by using krigging estimate ###
Niter <- 40000
theta <- matrix(hat,1)

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
COV <- diag(c(3e-10,0.01,0.1,0.001))                                 

# GLS bets coefficients 
beta.hat.gls <- c(coef(m)$trend1,coef(m)$trend2,coef(m)$trend3,coef(m)$trend4,coef(m)$trend5)

# MCMC run 
set.seed(1)
ptm <- proc.time()
parameter = GPmcmc(Niter,theta,COV,lhXZ,beta.hat.gls,c(coef(m)$range,coef(m)$sd2),th,y,r.X)
MCMCtime = proc.time()-ptm
MCMCtime



save.image(file = "RSVB2.RData") 




















