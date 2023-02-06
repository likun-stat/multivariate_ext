library(mvtnorm)
library(rmutil)

# Generating data based on Theorem 2.1 ------------------------------------
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

for (i in 1:N) {
  Z[i,] <- rmvnorm(1, mean=c(0,0,0,0), sigma=sigma)
  S[i,] <- rlevy(d, m=0, s=1)
  R <- (W)%*%S[i,]
  for (j in 1:d) {
    X[i,j] <- (R[j]^phi[j])*g(Z[i,j])
    A[i,j] <- g(Z[i,j])^(phi[j]/2)
  }
}


# Full Conditionals -------------------------------------------------------

g.inv <- function(x,R,phi){
  out1 <- 1-(R^phi)/x
  suppressWarnings(out2 <- sapply(c(out1),qnorm))
  return(out2)
}

log.post.s <- function(x,phi,s){
  R <- (W)%*%s
  out2 <- sapply(c(g.inv(x,R,phi)),dnorm)
  if(any(is.na(out2))){
    return(-Inf)
  }
  out3 <- ((R^phi)/x^2)/out2
  if(any(s<0)){
    out4 = 0
  }else{
    out4 <- sapply(c(s),plevy)
  }
  if(is.null(dim(g.inv(x,R,phi)))){
    out1 <- 0
    out <- out1 + sum(log(out3) + log(out4))
    return(out)
  }
  out1 <- dmvnorm(g.inv(x,R,phi), mean=c(0,0,0,0), sigma=sigma, log = T)
  out <- out1 + sum(log(out3) + log(out4))
  return(out)
}

log.post.phi <- function(x,phi,s){
  R <- (W)%*%s
  out2 <- sapply(c(g.inv(x,R,phi)),dnorm)
  if(any(is.na(out2))){
    return(-Inf)
  }
  out3 <- ((R^phi)/x^2)/out2
  if(is.null(dim(g.inv(x,R,phi)))){
    out1 <- 0
    out <- out1 + sum(log(out3))
    return(out)
  }
  out1 <- dmvnorm(g.inv(x,R,phi), mean=c(0,0,0,0), sigma=sigma, log = T)
  out <- out1 + sum(log(out3))
  return(out)
}

log.post.phi.total <- function(X,phi,S){
  n <- nrow(X)
  out1 <- array(0,n)
  for (i in 1:n) {
    out1[i] <- log.post.phi(X[i,],phi,S[i,])
  }
  return(sum(out1))
}


temp <- function(x){
  current.phi <- phi
  current.phi[1] <- x
  n <- nrow(X)
  out1 <- array(0,n)
  for (i in 1:n) {
    out1[i] <- log.post.phi(X[i,],current.phi,S[i,])
  }
  return(sum(out1))
}

phi.start <- seq(0.519,0.521,length.out = 100)
new.vec <- array(NA,100)
for (i in 1:100) {
  new.vec[i] <- temp(phi.start[i])
}


plot(phi.start,new.vec, type = "l")

# MH Algorithm ------------------------------------------------------------

#Arbitrary starting values

phi <-  c(0.52,0.56,0.6,0.65)
lambda1 <- 0.0001
N <- 100
M <- 10000
s.mc  <- array(NA,dim = c(M,N,4))
phi.mc <- array(NA,dim = c(M,4))
acc.s <- array(0, dim = c(M, N))

r.s  <- array(NA,dim = c(M,N))
r.phi  <- array(NA,M)
acc.phi  <- array(NA,M)


oldlambda <- c(0.0001,0.00000001)

for (i in 1:M) {
  for (j in 1:N){
    s.t <- as.vector(S[j,] + rmvnorm(1,c(0,0,0,0),diag(4))*0.001)
    r.s[i,j] <- exp(log.post.s(X[j,], phi, s.t) - log.post.s(X[j,], phi, S[j,]))
    if (is.na(r.s[i,j])){
      s.mc[i,j,] <- S[j,]
      next
    }
    if (runif(1) < r.s[i,j]){
      S[j,] <- s.t
      acc.s[i,j] <- 1
    }
    s.mc[i,j,] <- S[j,]
  }
  phi.t <- as.vector( phi + abs(rmvnorm(1,c(0,0,0,0),diag(4)*0.000001) ))
  r.phi[i] <- exp(log.post.phi.total(X, phi.t, S)/N - log.post.phi.total(X, phi, S)/N)
  if (is.na(r.phi[i])){
    phi.mc[i,] <- phi
    next
  }
  if (runif(1) < r.phi[i]){
    phi <- phi.t
    acc.phi[i] <- 1
  }
  phi.mc[i,] <- phi
  
}

pdf("s scatter plots lamda new.pdf",height=6,width=7.5)
par(mfrow = c(2,3),mar=c(4,4,5,3),xpd=TRUE)
plot(s.mc[,1,1],s.mc[,1,2])
plot(s.mc[,1,1],s.mc[,1,3])
plot(s.mc[,1,1],s.mc[,1,4])
plot(s.mc[,1,2],s.mc[,1,3])
plot(s.mc[,1,3],s.mc[,1,4])
dev.off()

pdf("s for M=10000, N=100 lamda new.pdf",height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=TRUE)
plot(s.mc[,1,1], type = "l")
plot(s.mc[,1,2], type = "l")
plot(s.mc[,1,3], type = "l")
plot(s.mc[,1,4], type = "l")
mtext("s (M=10000, N=100)", side = 3,line = - 2,outer = TRUE)
dev.off()

pdf("s[,10,] for M=10000, N=100 lamda new.pdf",height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=TRUE)
plot(s.mc[,2,1], type = "l")
plot(s.mc[,2,2], type = "l")
plot(s.mc[,2,3], type = "l")
plot(s.mc[,2,4], type = "l")
mtext("s[,10,] (M=10000, N=100)", side = 3,line = - 2,outer = TRUE)
dev.off()

pdf("phi for M=10000, N=100 lamda new.pdf",height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=TRUE)
plot(phi.mc[,1], type = "l")
plot(phi.mc[,2], type = "l")
plot(phi.mc[,3], type = "l")
plot(phi.mc[,4], type = "l")
mtext("phi (M=10000, N=100)", side = 3,line = - 2,outer = TRUE)
dev.off()




cov.test <- apply(s.mc,2, cov)
Sliced <- aperm(`dim<-`(t(cov.test), list(4, 4, 100)))
