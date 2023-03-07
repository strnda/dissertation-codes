N <- function(par, val, dist, norm, n.points, pp) {

  edf <- ECDF(x = val,
              pp = pp) ## get ecdf

  aux <- dim(x = edf)[1] ## get number of values

  if(aux != n.points) { ## select number of values for norms

    if(aux < n.points) {

      n.points <- aux
    }

    edf <- edf[seq(from = 1,
                   to = aux,
                   length.out = n.points),]
  }

  a <- getDistArg(dist = dist) ## get distribution arguments

  arg <- as.list(setNames(par, ## ready dist arg for do.call
                          a))
  if((norm == "N1") | (norm == "N2")) {

    FN <- edf$value ## empirical values
    F <- do.call(what = paste0("q", ## theoretical values
                               dist),
                 args = c(list(edf$p),
                          arg))


    if(norm == "N1") { ## N1

      out <- rMSE(x = F,
                  y = FN)
    }

    if(norm == "N2") { ## N2

      out <- MSE(x = F,
                 y = FN)
    }
  }

  if((norm == "N3") | (norm == "N4")) {

    Xi <- edf$p ## empirical probs
    Xu <- do.call(what = paste0("p", ## theoretical probs
                                dist),
                  args = c(list(edf$value),
                           arg))

    if(norm == "N3") { ## N3

      out <- rMSE(x = Xu,
                  y = Xi)
    }

    if(norm == "N4") { ## N4

      out <- MSE(x = Xu,
                 y = Xi)
    }
  }

  if(!is.finite(x = out)) {

    out <- 10E9 ## if the norm result is not finite return a big number
  }

  return(out)
}

fitDist <- function(data,
                    dist,
                    n.points,
                    norm,
                    constrain,
                    pp,
                    opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                "xtol_rel" = 1.0e-8,
                                "maxeval" = 1.0e4)) {

  a <- getDistArg(dist = dist)

  lwr <- NULL

  if (constrain) {

    if (any(a == "shape2")) {

      start <- c(rep(x = 1,
                     times = length(x = a) - 1),
                 .1)
      upr <- c(rep(x = Inf,
                   times = length(x = a) - 1),
               .48)
    } else {

      start <- rep(x = 1,
                   times = length(x = a))
      upr <- rep(x = Inf,
                 times = length(x = a))
    }

  } else {

    if (any(a == "shape2")) {

      start <- c(rep(x = 1,
                     times = length(x = a) - 1),
                 .1)
      upr <- rep(x = Inf,
                 times = length(x = a))
    } else {

      start <- rep(x = 1,
                   times = length(x = a))
      upr <- rep(x = Inf,
                 times = length(x = a))
    }
  }

  # fit <-  neldermead(x0 = start, ## parameter starting position
  #                    fn = N, ## function to minimize
  #                    val = data, ## input data
  #                    dist = dist, ## distribution to be fitted
  #                    n.points = n.points, ## number of points for fitting
  #                    norm = norm, ## norm used to fit
  #                    lower = lwr, ## lower bound
  #                    upper = upr) ## upper bound

  fit <- nloptr(x0 = start,
                eval_f = N,
                lb = lwr,
                ub = upr,
                val = data, ## input data
                dist = dist, ## distribution to be fitted
                pp = pp, ## plotting position
                n.points = n.points,
                norm = norm,
                opts = opts)


  out <- as.list(setNames(object = fit$solution, ## return list of dist arguments
                          nm = a))

  structure(.Data = out,
            dist = dist,
            edf = ECDF(data),
            nfo = fit,
            class = "fitDist")
}

getDistArg <- function(dist) {

  if (exists(paste0('d', dist))) { ## see whether CDF is defined

    x <- unlist(strsplit(deparse(args(paste0('q', dist)))[1], ',')) ## get string of all arguments
    out <- gsub(' |= [0-9]+', ## strip the arguments of the balast
                '',
                x[grep('function|lower.tail|log.p = FALSE|/',
                       x,
                       invert = TRUE)])

    return(out)
  } else {

    message('Distribution in not defined')
  }
}

ECDF <- function(x, pp = "weibull") {

  ## sort the values: x1 < x2 < x3 < ... < xN
  st <- sort(x = x)

  ## weibull plotting position

  if (pp == "weibull") {

    p <- rank(x = st,
              ties.method = "first") / (length(x = st) + 1)
  }

  ## weibull plotting position
  if (pp == "chegodaev") {

    p <- (rank(x = st,
               ties.method = "first") - .3) / (length(x = st) + .4)
  }

  ## weibull plotting position
  aux <- data.table(
    p = p,
    value = st,
    key = "value"
  )

  J <- value <- NULL

  out <- aux[J(unique(x = value)),
             mult = "last"]

  structure(.Data = out)
}

rMSE <- function(x, y) {

  sum((x/y - 1)^2)/length(y) ## ratio MSE
}

MSE <- function(x, y) {

  sum((x - y)^2)/length(y) ## ordinary MSE
}

AD_test <- function(data, dist, ...) {
  
  val <- unlist(x = data)
  val <- val[!is.na(x = val)]
  
  pars <- list(...)
  
  p <- do.call(what = paste0("p", dist),
               args = c(list(q = val),
                        pars))
  
  u <- sort(x = p[(p != 0) & (p != 1)])
  nr <- length(x = u)
  -nr - 1/nr * sum((2 * 1:nr - 1) * log(x = u) +
                     (2 * nr - 2 * 1:nr + 1) * log(x = 1 - u))
  
}

pggamma <- function(q, scale, shape1, shape2, lower.tail = TRUE, log.p = FALSE) { ## cdf
  
  if((shape1 <= 0) | (shape2 <= 0)) {
    
    return(NaN)
  } else {
    
    p <- pgamma((q/scale)^shape2, scale = 1, shape = shape1/shape2)
    
    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}
    
    return(p)
  }
}

qggamma <- function(p, scale, shape1, shape2, lower.tail = TRUE, log.p = FALSE) { ## cdf
  
  if((shape1 <= 0) | (shape2 <= 0)) {
    
    return(NaN)
  } else {
    
    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}
    
    q <- scale*qgamma(p, scale = 1, shape = shape1/shape2)^(1/shape2)
    
    return(q)
  }
}

smp <- function(x, n = 1000, method = 1) {
  
  x <- as.data.table(x = x)
  
  if (method == 1) {
    
    out <- lapply(X = 1:n, 
                  FUN = function(i) {
                    x[base::sample(x = 1:nrow(x = x), 
                                   size = nrow(x = x),
                                   replace = TRUE)]
                  }
    )
  } else {
    
    out <- lapply(X = 1:n, 
                  FUN = function(i) {
                    
                    p0 <- length(x = which(x = is.na(x = x))) / length(x = unlist(x = x[, -1]))
                    
                    
                    e <- try(expr = gs <- rmvnorm(n = nrow(x = x),
                                                  mean = rep(x = 0,
                                                             times = (ncol(x = x) - 1)),
                                                  sigma = cor(x = x[, -1],
                                                              use = "pairwise.complete.obs")), 
                             silent = TRUE)
                    
                    if (inherits(x = e, what = "try-error")) {
                      
                      gs <- matrix(data = rnorm(n = nrow(x = x) * (ncol(x = x) - 1)),
                                   ncol = ncol(x = x) - 1,
                                   nrow = nrow(x = x))
                    }
                    
                    
                    gs <- as.data.table(x = cbind(seq(to = year(x = Sys.Date()), 
                                                      by = 1,
                                                      length.out = nrow(x = gs)),
                                                  gs))
                    names(x = gs) <- names(x = x)
                    r_fit <- r_fit_aux <- regsamlmu(x = x[, -1])
                    
                    r_fit_aux$t <- r_fit$l_1 * r_fit$t
                    
                    para <- apply(X = r_fit_aux[,-1:-2], 
                                  MARGIN = 1, 
                                  FUN = pelgev)
                    para <- as.data.table(x = t(x = para))
                    para[, site := r_fit$name]
                    
                    p_m <- melt(data = gs, 
                                id.vars = 1, 
                                variable.name = "site", 
                                value.name = "g")
                    
                    p_m <- merge(x = p_m, 
                                 y = para,
                                 by = "site")
                    
                    p_m[, p := pnorm(q = g)]
                    
                    if (method == 3) {
                      
                      p_m[, p := (p - p0)/(1 - p0)]
                      p_m[p < 0, p := NA]
                    } 
                    
                    p_m[, value := extraDistr::qgev(p = p,
                                                    mu = xi,
                                                    sigma = alpha,
                                                    xi = -k)]
                    
                    mx <- dcast(data = p_m,
                                formula = YR ~ site,
                                value.var = "value")
                    
                    if (method == 2) {
                      
                      mx[is.na(x = x)] <- NA
                    } 
                    
                    mx
                  }
    )
    
  }
}

fitsmp <- function(s) {
  
  out <- lapply(
    X = s, 
    FUN = function(x) {
      
      e <- try(expr = {
        
        r_mom <- regsamlmu(x = x[, -1])
        
        r_para <- t(x = apply(X = r_mom[,-1:-2], 
                              MARGIN = 1, 
                              FUN = pelgev))
        
        at_site_q <- t(
          x = apply(X = r_para, 
                    MARGIN = 1, 
                    FUN = function(x) {
                      quagev(f = 1 - (1 / c(5, 10, 20, 50, 100)),
                             para = x)
                    }
          )
        )
        
        r_fit <- regfit(regdata = r_mom, 
                        dist = "gev")
        
        at_site_q_reg <- sitequant(f = 1 - (1 / c(5, 10, 20, 50, 100)), 
                                   rfd = r_fit)
        
        dif <- apply(X = at_site_q, 
                     MARGIN = 2, 
                     FUN = IQR) / 
          apply(X = at_site_q_reg, 
                MARGIN = 2, 
                FUN = IQR)
        
      })
      
      if (inherits(x = e, 
                   what = "try-error")) {
        
        NA
      } else {
        
        c(r_fit$para, dif)
      }
    }
  )
  
  out <- do.call(what = rbind,
                 args = out)
  
  as.data.table(x = out)
}
