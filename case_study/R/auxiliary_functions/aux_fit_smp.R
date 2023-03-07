sim <- function(extremes, dist = 'gpa', trim = c(0, 0)) { # stationary index flood method ;)
  
  l.atsite <- t(apply(extremes, 2, function(x, ...) samlmu(x[is.finite(x)], trim = trim)))
  para <- data.table(.id = names(extremes), t(apply(l.atsite, 1, function(x) {do.call(paste0('pel', dist), list(lmom = x))})))

  if(dim(l.atsite)[1] > 1)  {
    
    l.atsite[,2] <- l.atsite[,2]/l.atsite[,1]
    
    # w <- apply(extremes, 2, function(x) length(x[!is.na(x)]))
    # w.l <- apply(l.atsite[,-1], 2, weighted.mean, w = w)
    w.l <- apply(l.atsite[,-1], 2, mean)
    
    structure(.Data = list(data = as.data.frame(extremes), 
                           scaling_factor = l.atsite[,1], 
                           REG = do.call(paste0(dist, '.para'), args = list(c(1, w.l))), 
                           dist = dist,
                           para = para),
              sim.call = match.call(),
              class = 'sim')
  } else {
    
    structure(.Data = list(data = as.data.frame(extremes), 
                           scaling_factor = 1, 
                           REG = do.call(paste0(dist, '.para'), args = list(l.atsite)), 
                           dist = dist,
                           para = para),
              sim.call = match.call(),
              class = 'sim')
  }
}

# zero - r..., nezavis.
# one - NA jsou nepar. nasamplovany, radky rmvnorm, hodnoty jsou po sloupcich r..., orankovani
# two - cela matice rmvnorm, NA jsou vlozeny do matice, zbytek jako metoda one
# three - cela matice rmvnorm, NA jsou vlozeny do matice, hodnoty jsou pres residua prevedeny na gev/gpa

sample <- function(...) {UseMethod('sample')}

sample.default <- base::sample

sample.sim <- function(model_object, length = 1, type = 'nonpar', na = T) { # testovaci nonpar sample (pro test fitovani)
  
  dta <- as.data.table(model_object$data)
  para <- model_object$para
  dist <- attr(model_object, 'sim.call')$dist
  
  i <- 1:length
  
  if(type == 'nonpar') {
    
    out <- mapply(function(i) {
      
      m <- resid.sim(model_object)
      m <- m[sample(1:dim(m)[1], dim(m)[1], replace = T),]
      m <- backtodata(m, model_object)
      return(m)
    }, i, SIMPLIFY = FALSE)
  }
  
  if(type == 'zero') {
    
    out <- mapply(function(i) {
      
      m <- suppressWarnings(melt(dta, variable.name = '.id'))
      m <- data.table(para)[m, on = c('.id')]
      
      m[!is.na(value), value := do.call(paste0('r',dist), list(.N, xi = unique(xi), alpha = unique(alpha), k = unique(k))), by = .id]
      m[, aux := 1:.N, by = .id]
      
      return(dcast(m[, .(.id, value, aux)], aux ~ .id, value.var = 'value')[,-1])
      
    }, i, SIMPLIFY = FALSE)
  }
  
  if(type == 'para_cor') {
    
    out <- mapply(function(i) {
      cvr <- cov(resid.sim(model_object), use = 'pairwise.complete.obs')
      
      crr <- cov2cor(cvr)
      crr[] <- mean(crr[upper.tri(crr)], na.rm = TRUE)
      diag(crr) <- 1
      sdMat <- diag(sqrt(diag(cvr)))
      cvr <- sdMat %*% crr %*% t(sdMat)
      
      norm.res <- data.table(rmvnorm(dim(dta)[1], sigma = cvr))
      names(norm.res) <- names(dta)
      
      if(na) {
        
        norm.res <- norm.res*(dta/dta) 
      }
      
      m <- suppressWarnings(melt(norm.res, variable.name = '.id'))
      m <- data.table(para)[m, on = c('.id')]
      
      m[!is.na(value), value := do.call(paste0('r',dist), list(.N, xi = unique(xi), alpha = unique(alpha), k = unique(k))), by = .id]
      m[, aux := 1:.N, by = .id]
      
      return(dcast(m[, .(.id, value, aux)], aux ~ .id, value.var = 'value')[,-1])
      
    }, i, SIMPLIFY = FALSE)
  }
  
  structure(.Data = out,
            names = paste0('b_sample_',i),
            class = 'simsample',
            type = type)

}

fit <- function(...) {UseMethod('fit')}

fit.simsample <- function(smp, model_object) {
  
  cl <- attr(model_object, 'sim.call')

  structure(.Data = lapply(1:length(smp), function(x) {cl$extremes <- smp[[x]]; eval(cl)}),
            class = 'simsample',
            type = attributes(smp)$type)
}
