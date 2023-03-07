# library(nim)
# 
# data('precip_max')
# 
# dta <- precip_max[,c(1,12)]
# extremes(dta) = -1
# 
# n <- nim( ~1, data = dta)
# 
# smp <- sample(n, length = 50)
# f <- fit(n, smp)

quantile.nim <- function(nim, p = NULL, Tm = c(2, 5, 10, 50), at_site = TRUE){
  
  cl = match.call()
  qO = Vectorize(function(p){
    r = data.frame(XI*(1 - G/K*(1 -(-log(p))^(-K))))
    if (model_info(nim)$cvrt == 'I') {
      r = data.frame(I = 1, r[1, ])
      if (is.null(names(XI))) {
        names(r) = c(names(r)[1], 'value')
      } else {
        names(r) = c(names(r)[1], names(XI))#if (!at_site) {names(r)[2] = 'regional'}
      }
    } else {
      r = data.frame(nim$REG[[1]], r, check.names = FALSE)
      names(r)[1] = model_info(nim)$cvrt
      names(r)[2:ncol(r)] = if (at_site) (names(nim$XI)) else ('regional')
    }
    return(r)
  }, SIMPLIFY = FALSE)
  
  if (is.null(p)) p = 1 - 1 / Tm
  
  if(at_site) {
    XI <- outer(regional(nim)$XI, atsite(nim))
    G <- matrix(data = rep(x = exp(nim$REG$G), times = dim(XI)[2]),
                nrow = dim(XI)[1])
    K <- matrix(data = rep(x = nim$REG$K, times = dim(XI)[2]),
                nrow = dim(XI)[1])
  } else {
    XI <- regional(nim)$XI
    G <- exp(nim$REG$G)
    K <- nim$REG$K
  }
  
  res = qO(p)
  names(res) = p
  res = rbindlist(res, idcol = 'p')
  res
}

AD.test <- function(val, location = 0, scale = 1, shape = 0, dist = 'gev') {
  
  val <- unlist(val)
  val <- val[!is.na(val)]
  
  p <- do.call(what = paste0('p', dist), args = list(q = val, xi = location, alpha = scale, k = shape))
  
  u <- sort(p[(p != 0) & (p != 1)])
  nr <- length(u)
  -nr - 1/nr*sum((2*1:nr - 1)*log(u) + (2*nr - 2*1:nr + 1)*log(1 - u))
  
}

AD.regional <- function(dta, a = .1,  N, type = 'para_cor', save.samples = NULL) {
  
  s <- sample(dta, length = N,  type = 'para_cor', na = T)
  
  if(!is.null(save.samples)) assign(as.character(save.samples), s, envir = .GlobalEnv)
  
  ss <- lapply(s, function(x) {x <- as.data.frame(x); names(x) <- names(MX); x})
  
  para.s <- rbindlist(lapply(ss, function(f) {as.data.table(merge(suppressMessages(melt(f)),
                                                                 data.table(SP_ID = names(f), t(apply(f, 2, function(x) pelgpa(samlmu(x))))), 
                                                                 by.x = 'variable',
                                                                 by.y = 'SP_ID'))}), idcol = T)
  ad.s <- dcast(para.s[, .(ad.val = AD.test(val = value,
                                            location = unique(xi), 
                                            scale = unique(alpha), 
                                            shape = unique(k), 
                                            dist = 'gpa')),
                       by = .(.id, variable)], variable ~ .id, value.var = 'ad.val')
  
  ad.res <- data.table(t(ad.s[,-1]))
  setnames(ad.res, ad.s[, as.character(variable)])

  ad.res[, I := 1:nrow(ad.res)]
  
  mad <- melt(ad.res, id.vars = 'I')
  mad[, RNK := rank(value), by = variable]
  mad[, MX := max(RNK), by = I]
  
  x <- 1:N
  aux <- mad[, mapply(function(x) length(which(dta$data >= x))/.N, x = x)]
  
  gp <- mad[, min(value[RNK == which(abs(aux - a) == min(abs(aux - a)))]), by = variable][, mean(V1)]
  
  return(gp)    
}

n.max <- function(x, n = 2){
  l <- length(x)
  if (n > l){
    warning('vole !!!')
    n <- length(x)
  }
  if (n == 0) {
    x[NULL]
  } else {
    sort(x)[(l - n + 1):l]
  }
}

n.min <- function(x, n = 2){
  l <- length(x)
  if(n > l){
    warning('vole !!!')
    n <- length(x)
  }  
  if (n == 0) {
    x[NULL]
  } else {
    sort(x)[1:n]
  }
}

# pelgpa(samlmu(x))
# gpa.para(moje.lm(x))

################################### resid

resid.sim <- function(model_object){
  
  dta <- model_object$data
  para <- model_object$para
  sf <- model_object$scaling_factor
  
  resi <- suppressMessages(melt(dta))
  resi <- data.table(variable = names(dta), sf = sf, para)[resi, on = c('variable')]
  resi <- resi[, `:=`(resi = 1/-k*log(1 + -k*((value/sf)/alpha)),
                      aux = seq_len(.N)), by = 'variable'] # coles - sigma = alpha, xi = -kappa #######
  
  dcast(resi, aux ~ variable, value.var = 'resi')[,-1]
}

backtodata <- function(res, model_object){
  
  para <- model_object$para
  sf <- model_object$scaling_factor
  
  res <- suppressWarnings(melt(res))
  res <- data.table(variable = names(model_object$data), sf = sf, para)[res, on = c('variable')]
  dta <- res[, `:=`(val = ((alpha*(exp(value*-k) - 1))/-k)*sf,
                    aux = seq_len(.N)), by = 'variable']
  
  dcast(dta, aux ~ variable, value.var = 'val')[,-1]
}

synth.data <- function(model_object, dependent = T, use.NA = T, sample_para = T, dist = NULL) {
  
  dta <- as.data.table(model_object$data)
  if(is.null(dist)) {dist <- attr(model_object, 'sim.call')$dist}
  
  
  if(sample_para) {
    
    n <- dim(dta)[2]
    
    reg.mom <- do.call(paste0('lmr',dist), list(model_object$REG))
    
    at.site.mom <- data.table(l_1 = sample(seq(min(model_object$scaling_factor), max(model_object$scaling_factor), .1), n, replace = T), l_2 = reg.mom[2], t_3 = reg.mom[3])
    at.site.mom[, l_2 := l_1*l_2] # !!!! Filipe !!!!
    
    para <- as.data.table(t(apply(at.site.mom, 1, noquote(paste0(dist,'.para')))))
    para[, .id := dimnames(dta)[2]]
  
  } else {
    
    para <- model_object$para
  }
  
  if(dependent) {
    
    cvr <- cov(resid.sim(model_object), use = 'pairwise.complete.obs')
    
    crr <- cov2cor(cvr)
    crr[] <- mean(crr[upper.tri(crr)], na.rm = TRUE)
    diag(crr) <- 1
    sdMat <- diag(sqrt(diag(cvr)))
    cvr <- sdMat %*% crr %*% t(sdMat)
    
    norm.res <- data.table(rmvnorm(dim(dta)[1], sigma = cvr))
    names(norm.res) <- names(dta)
    
    if(use.NA) {
      
      norm.res <- norm.res*(dta/dta) 
    }
    
    m <- suppressWarnings(melt(norm.res, variable.name = '.id'))
    m <- data.table(para)[m, on = c('.id')]
    
    m[!is.na(value), value := do.call(paste0('r',dist), list(.N, xi = unique(xi), alpha = unique(alpha), k = unique(k))), by = .id]
    m[, aux := 1:.N, by = .id]
    
  } else {
    
    m <- suppressWarnings(melt(dta, variable.name = '.id'))
    m <- data.table(para)[m, on = c('.id')]
    
    if(use.NA) {
      
      m[!is.na(value), value := do.call(paste0('r',dist), list(.N, xi = unique(xi), alpha = unique(alpha), k = unique(k))), by = .id]
    } else {
      
      m[, value := do.call(paste0('r',dist), list(.N, xi = unique(xi), alpha = unique(alpha), k = unique(k))), by = .id]
    }
    
    m[, aux := 1:.N, by = .id]
    
  }
  
  return(dcast(m[, .(.id, value, aux)], aux ~ .id, value.var = 'value')[,-1])
}

lseq <- function(from = 1E-5, to = 1, length.out = 6) {
  exp(seq(log(from), log(to), length.out = length.out))
}
