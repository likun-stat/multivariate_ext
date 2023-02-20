setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/multivariate_ext/")
source("./utils.R")
library(mvtnorm)
library(rmutil)

# --------------------------- Generating data based on Theorem 2.1 ---------------------------
d <- 4
N <- 100
sigma <- matrix(c(1,0.3,0.4,0.2,
                  0.3,1,0.5,0.3,
                  0.4,0.5,1,0.6,
                  0.2,0.3,0.6,1), ncol=d)
W <- matrix(c(1,0,0,0,
              1,1,0,0,
              0,1,1,0,
              1,0,0,1), nrow=d)
phi <- c(0.52,0.56,0.6,0.65)


g <- function(Z){
  x <- 1/(1-pnorm(Z))
  return(x)
}

Z <- array(NA,dim = c(N,d))
X <- array(NA,dim = c(N,d))
A <- array(NA,dim = c(N,d))
S <- array(NA,dim = c(N,d))

set.seed(123)
for (i in 1:N) {
  Z[i,] <- rmvnorm(1, mean=c(0,0,0,0), sigma=sigma)
  S[i,] <- rlevy(d, m=0, s=1)
  R <- (W)%*%S[i,]
  for (j in 1:d) {
    X[i,j] <- (R[j]^phi[j])*g(Z[i,j])
    A[i,j] <- g(Z[i,j])^(phi[j]/2)
  }
}



# --------------------------- Check full conditionals ---------------------------
index = 1
phi.start1 <- seq(phi[index]-0.06, phi[index]+0.03, length.out = N)
new.vec1 <- array(NA,N)
for (i in 1:N) {
  new.vec1[i] <- temp(phi.start1[i], index=index)
}

index = 2
phi.start2 <- seq(phi[index]-0.003, phi[index]+0.003, length.out = N)
new.vec2 <- array(NA,N)
for (i in 1:N) {
  new.vec2[i] <- temp(phi.start2[i], index=index)
}

index = 3
phi.start3 <- seq(phi[index]-0.003, phi[index]+0.02, length.out = N)
new.vec3 <- array(NA,N)
for (i in 1:N) {
  new.vec3[i] <- temp(phi.start3[i], index=index)
}

index = 4
phi.start4 <- seq(phi[index]-0.003, phi[index]+0.02, length.out = N)
new.vec4 <- array(NA,N)
for (i in 1:N) {
  new.vec4[i] <- temp(phi.start4[i], index=index)
}

par(mfrow = c(2,2))
plot(phi.start1,new.vec1, type = "l", main=expression(phi[1]))
plot(phi.start2,new.vec2, type = "l", main=expression(phi[2]))
plot(phi.start3,new.vec3, type = "l", main=expression(phi[3]))
plot(phi.start4,new.vec4, type = "l", main=expression(phi[4]))
par(mfrow = c(1,1))


index = 2; row=100
s.start1 <- seq(S[row,index]-5, S[row,index]+18, length.out = N)
new.vec1 <- array(NA,N)
for (i in 1:N) {
  new.vec1[i] <- tmp.s(s.start1[i], index=index, row=row)
}
plot(s.start1, new.vec1, type = "l", main=expression(phi[1]))
abline(v=S[row,index], col='red')

s.start1 <- seq(1, S[row,index]+5, length.out = N)
new.vec1 <- array(NA,N)
for (i in 1:N) {
  new.vec1[i] <- tmp.s(s.start1[i], index=index, row=row)
}
plot(s.start1,new.vec1, type = "l", main=expression(phi[1]))



# --------------------------- MH Algorithm ----------------------------
# -------------------------- Adaptive  MCMC ---------------------------

cov.s.fun <- function(dat){
  cov.s.mat <- apply(dat, 2, cov)
  cov.s.list <- list()
  for( iter in 1:ncol(cov.s.mat)){
    cov.s.list[[iter]] <- matrix(cov.s.mat[,iter],4,4)
  }
  return(cov.s.list)
}

cov.s <- list()
for( iter in 1:N){
  cov.s[[iter]] <- diag(4)
}

cov.phi <- diag(4)
#cov.phi <- cov(phi.mc)


phi <-  c(0.52,0.56,0.6,0.65)
M <- 1e6
k <- 100
s.mc  <- array(NA,dim = c(M,N,4))
phi.mc <- array(NA,dim = c(M,4))
acc.s <- array(0, dim = c(M, N))
acc.phi  <- array(0,M)
r.s  <- array(NA,dim = c(M,N))
r.phi  <- array(NA,M)

d <- 4
c.0 <- 1
c.1 <- 0.8
r.opt <- 0.234
sig.phi <- 0.00001
sig.s <- array(2.4^2/d,N)
t <- 1


for (i in 1:M) {
  if(i %% 10000 == 0) cat("Done with", i, "iterations\n")
  for (j in 1:N){
    s.t <- as.vector(S[j,] + rmvnorm(1,c(0,0,0,0),sig.s[j]*cov.s[[j]]))
    r.s[i,j] <- exp(log.post.s(X[j,], phi, s.t) - log.post.s(X[j,], phi, S[j,]))
    if (is.na(r.s[i,j])){
      s.mc[i,j,] <- S[j,]
    }
    if (!is.na(r.s[i,j]) & runif(1) < r.s[i,j]){
      S[j,] <- s.t
      acc.s[i,j] <- 1
    }
    s.mc[i,j,] <- S[j,]
  }
  phi.t <- as.vector( phi + rmvnorm(1,c(0,0,0,0),sig.phi*cov.phi))
  r.phi[i] <- exp(log.post.phi.total(X, phi.t, S) - log.post.phi.total(X, phi, S))
  if (is.na(r.phi[i])){
    phi.mc[i,] <- phi
  }
  if (!is.na(r.phi[i]) & runif(1) < r.phi[i]){
    phi <- phi.t
    acc.phi[i] <- 1
  }
  phi.mc[i,] <- phi
  
  #adaptive part
  if(i %% k == 0){
    gamma.1 <- 1/(t+3)^c.1
    gamma.2 <- c.0*gamma.1
    
    r.phi.est <- sum(acc.phi[(i-k+1):i], na.rm = TRUE)/k
    cov.phi.est <- cov(phi.mc[(i-k+1):i, ])
    cov.phi <- cov.phi + gamma.1*(cov.phi.est - cov.phi)
    log.sig.sq.new.phi <- log(sig.phi^2) + gamma.2*(r.phi.est-r.opt)
    sig.phi <- sqrt(exp(log.sig.sq.new.phi))
    
    
    r.s.est <- apply(acc.s[(i-k+1):i,],2, function(x) sum(x,na.rm = TRUE)/k)
    cov.s.est <- cov.s.fun(s.mc[(i-k+1):i,,])
    cov.s.list <- list()
    for (iter in 1:length(cov.s.est)) {
      cov.s.list[[iter]] <- cov.s[[iter]] + gamma.2*(cov.s.est[[iter]] - cov.s[[iter]])
    }
    cov.s <- cov.s.list
    log.sig.sq.new.s <- log(sig.s^2) + gamma.2*(r.s.est-r.opt)
    sig.s <- sqrt(exp(log.sig.sq.new.s))
    
    t <- t + 1
  }
}

pdf(paste0("phi with M = ",M," block size = ",k," N = ",N,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=TRUE)
plot(phi.mc[,1], type = "l")
plot(phi.mc[,2], type = "l")
plot(phi.mc[,3], type = "l")
plot(phi.mc[,4], type = "l")
mtext(paste0("phi with M = ",M," block size = ",k," N = ",N), side = 3,line = - 2,outer = TRUE)
dev.off()

pdf(paste0("phi with 5000 burnin M = ",M," block size = ",k," N = ",N,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=TRUE)
plot(phi.mc[5000:M,1], type = "l")
plot(phi.mc[5000:M,2], type = "l")
plot(phi.mc[5000:M,3], type = "l")
plot(phi.mc[5000:M,4], type = "l")
mtext(paste0("phi with M = ",M," block size = ",k," N = ",N), side = 3,line = - 2,outer = TRUE)
dev.off()

pdf(paste0("s with M = ",M," block size = ",k," N = ",N,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=TRUE)
plot(s.mc[,2,1], type = "l")
plot(s.mc[,2,2], type = "l")
plot(s.mc[,2,3], type = "l")
plot(s.mc[,2,4], type = "l")
mtext(paste0("s with M = ",M," block size = ",k," N = ",N), side = 3,line = - 2,outer = TRUE)
dev.off()

phi.mc[1:10,1]
