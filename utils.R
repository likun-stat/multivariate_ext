# Full Conditionals -------------------------------------------------------
g.inv <- function(x,R,phi){
  out1 <- 1-(R^phi)/x
  suppressWarnings(out2 <- sapply(c(out1),qnorm))
  return(out2)
}

library(mvtnorm)
log.post.s <- function(x,phi,s){
  R <- (W)%*%s
  Z_tmp <- g.inv(x,R,phi)
  out2 <- sapply(c(Z_tmp),dnorm)
  if(any(is.na(out2))){
    return(-Inf)
  }
  out3 <- ((R^phi)/x^2)/out2
  if(any(s<0)){
    out4 = 0
  }else{
    out4 <- sapply(c(s),dlevy)
  }
  # if(is.null(dim(Z_tmp))){
  #   out1 <- 0
  #   out <- out1 + sum(log(out3) + log(out4))
  #   return(out)
  # }
  out1 <- dmvnorm(Z_tmp, mean=c(0,0,0,0), sigma=sigma, log = TRUE)
  out <- out1 + sum(log(out3) + log(out4))
  return(out)
}

log.post.phi <- function(x,phi,s){
  if(any(phi<0.45 | phi >0.7)){
    return(-Inf)
  }
  R <- (W)%*%s
  Z_tmp <- g.inv(x,R,phi)
  out2 <- sapply(c(Z_tmp),dnorm)
  if(any(is.na(out2))){
    return(-Inf)
  }
  out3 <- ((R^phi)/x^2)/out2
  # if(is.null(dim(Z_tmp))){
  #   out1 <- 0
  #   out <- out1 + sum(log(out3))
  #   return(out)
  # }
  out1 <- dmvnorm(Z_tmp, mean=c(0,0,0,0), sigma=sigma, log = TRUE)
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


temp <- function(x, index){
  current.phi <- phi
  current.phi[index] <- x
  n <- nrow(X)
  out1 <- array(0,n)
  for (i in 1:n) {
    out1[i] <- log.post.phi(X[i,],current.phi,S[i,])
  }
  return(sum(out1))
}

tmp.s <- function(s_x, index, row){
  x <- X[row,]
  s <- S[row,]
  s[index] <- s_x
  R <- (W)%*%s
  Z_tmp <- g.inv(x,R,phi)
  out2 <- sapply(c(Z_tmp),dnorm)
  if(any(is.na(out2))){
    return(-Inf)
  }
  out3 <- ((R^phi)/x^2)/out2
  if(any(s<0)){
    out4 = 0
  }else{
    out4 <- sapply(c(s),dlevy)
  }
  # if(is.null(dim(Z_tmp))){
  #   out1 <- 0
  #   out <- out1 + sum(log(out3) + log(out4))
  #   return(out)
  # }
  out1 <- dmvnorm(Z_tmp, mean=c(0,0,0,0), sigma=sigma, log = TRUE)
  out <- out1 + sum(log(out3) )
  return(out)
}
