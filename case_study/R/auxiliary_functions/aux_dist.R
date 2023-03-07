######## GEV ########

dgev <- function (x, xi = 0, alpha = 1, k = 0, log = FALSE) {
  
  x <- (x - xi)/alpha
  
  if (k == 0) {
    d <- log(1/alpha) - x - exp(-x)
  } else {
    n <- length(x)
    x <- 1 + k * x
    i <- x > 0
    alpha <- rep(alpha, length.out = n)[i]
    d[i] <- log(1/alpha) - x[i]^(-1/k) - (1/k + 1)*log(x[i])
    d[!i] <- -Inf
  }
  
  if (!log) {
    d <- exp(d)
  } 
    
  d
}

pgev <- function(q, xi = 0, alpha = 1, k = 0, lower.tail = TRUE) {
  
  if (k == 0) {
    p <- (q - xi)/alpha
  } else {
    p <- -1/k*log(pmax(0, 1 - k*(q - xi)/alpha))
  }
  
  if (!lower.tail) {
    p <- 1 - p
  }
  
  return(exp(-exp(-p)))
}

qgev <- function(p, xi = 0, alpha = 1, k = 0, lower.tail = TRUE) {
  
  if (!lower.tail) {
    p <- 1 - p
  }
    
  if (k == 0) {
    out <- xi - alpha*log(-log(p))
  } else {
    out <- xi + alpha/k*(1 - (-log(p))^k)
  }
  
  return(out)
}

rgev <- function(n, xi = 0, alpha = 1, k = 0) {
  
  qgev(runif(n), xi, alpha, k)
}

######## GPA ########

dgpa <- function (x, xi = 0, alpha = 1, k = 0, log = FALSE) {

  d <- (x - xi)/alpha
  
  n <- length(d)
  alpha <- rep(alpha, length.out = n)
  i <- (d > 0 & ((1 + k*d) > 0))
  
  if (k == 0) {
    d[i] <- log(1/alpha[i]) - d[i]
    d[!i] <- -Inf
  } else {
    d[i] <- log(1/alpha[i]) - (1/k + 1)*log(1 + k*d[i])
    d[!i] <- -Inf
  }
  
  if (!log) {
    d <- exp(d)
  }
  
  d
}

pgpa <- function(q, xi = 0, alpha = 1, k = 0, lower.tail = TRUE) {
  
  if (k == 0) {
    p <- (q - xi)/alpha
  } else {
    p <- -1/k*log(pmax(0, 1 - k*(q - xi)/alpha))
  }
  
  if (!lower.tail) {
    p <- 1 - p
  }
  
  1 - exp(-pmax(p, 0))
}

qgpa <- function(p, xi = 0, alpha = 1, k = 0, lower.tail = TRUE) {
  
  if (!lower.tail) {
    p <- 1 - p
  }
  
  if (k == 0) {
    out <- xi + alpha*(-log(1 - p))
  } else {
    out <- xi + alpha/k*(1 - (1 - p)^(k))
  }
  
  return(out)
}

rgpa <- function(n, xi = 0, alpha = 1, k = 0) {
  
  qgpa(runif(n), xi, alpha, k)
}

######## multivariate normal ########

rmvnorm <- function (n, sigma) {
  
  r <- chol(sigma, pivot = TRUE)
  r[, order(attr(r, 'pivot'))]
  
  out <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% r
  
  out
}

######## L-mom ########

lmanual <- function (x) {
  
  x <- sort(x)
  
  n <- length(x)
  nn <- rep(n - 1, n)
  
  pp <- seq(0, n - 1)
  p1 <- pp/nn
  p2 <- p1*(pp - 1)/(nn - 1)
  p3 <- p2*(pp - 2)/(nn - 2)
  
  b0 <- sum(x)/n
  b1 <- sum(p1*x)/n
  b2 <- sum(p2*x)/n
  b3 <- sum(p3*x)/n
  
  l1 <- b0
  l2 <- 2*b1 - b0
  l3 <- 2*(3*b2 - b0)/(2*b1 - b0) - 3
  l4 <- 5*(2*(2*b3 - 3*b2) + b0)/(2*b1 - b0) + 6
  
  unlist(mget(paste0('l', 1:4)))
}

gpa.para <- function(l) {
  
  k <- (1 - 3*l[3])/(1 + l[3])
  alpha <- (1 + k)*(2 + k)*l[2]
  xi <- l[1] - (2 + k)*l[2]
  
  unlist(setNames(c(xi, alpha, k), c('xi', 'alpha', 'k')))
}

gev.para <- function(l) {
  
  c <- 2/(3 + l[3]) - log(2)/log(3)
  
  k <- 7.859*c + 2.9554*c^2
  alpha <- (l[2]*k)/((1 - 2^(-k))*gamma(1 + k))
  xi <- l[1] - alpha*(1 - gamma(1 + k))/k
  
  unlist(setNames(c(xi, alpha, k), c('xi', 'alpha', 'k')))
}
