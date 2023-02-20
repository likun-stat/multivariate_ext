# --------- Arbitrary starting values ---------

phi <-  c(0.52,0.56,0.6,0.65)
lambda1 <- 0.0001
M <- 10000
s.mc  <- array(NA,dim = c(M,N,4))
phi.mc <- array(NA,dim = c(M,4))
acc.s <- array(0, dim = c(M, N))

r.s  <- array(NA,dim = c(M,N))
r.phi  <- array(NA,M)
acc.phi  <- array(NA,M)


oldlambda <- c(0.001,0.00001)

for (i in 1:M) {
  if(i %% 500 == 0) cat("Done with", i, "iterations\n")
  for (j in 1:N){
    s.t <- as.vector(S[j,] + rmvnorm(1,c(0,0,0,0),diag(4))*0.01)
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
  phi.t <- as.vector( phi + rmvnorm(1,c(0,0,0,0),diag(4)*0.00005))
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

pdf("phi for M=10000, N=100 lamda new.pdf",height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=TRUE)
plot(phi.mc[,1], type = "l")
plot(phi.mc[,2], type = "l")
plot(phi.mc[,3], type = "l")
plot(phi.mc[,4], type = "l")
mtext("phi (M=10000, N=100)", side = 3,line = - 2,outer = TRUE)
dev.off()

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




# --------- Proposal covariance ---------
choose <- tail(1:M, 5000)
cov.phi <- cov(phi.mc[choose, ])
cov.s.mat <- apply(s.mc, 2, function(x) cov(x[choose, ]))
cov.s <- list()
for( iter in 1:ncol(cov.s.mat)){
  cov.s[[iter]] <- matrix(cov.s.mat[,iter],4,4)
}

lambda.s <- 0.5
lambda.phi <- 0.5

for (i in 1:100000) {
  if(i %% 500 == 0) cat("Done with", i, "iterations\n")
  for (j in 1:N){
    s.t <- as.vector(S[j,] + rmvnorm(1,c(0,0,0,0),cov.s[[j]])*lambda.s)
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
  phi.t <- as.vector( phi + rmvnorm(1,c(0,0,0,0),cov.phi)*lambda.phi)
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

pdf("phi for M=10000, N=100 lamda new.pdf",height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=TRUE)
plot(phi.mc[,1], type = "l")
plot(phi.mc[,2], type = "l")
plot(phi.mc[,3], type = "l")
plot(phi.mc[,4], type = "l")
mtext("phi (M=10000, N=100)", side = 3,line = - 2,outer = TRUE)
dev.off()

