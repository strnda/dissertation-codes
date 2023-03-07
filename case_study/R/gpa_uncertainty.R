library(data.table); library(ggplot2); library(lmom)

source('R/auxiliary_functions/imports.R')

mx <- data.table(readRDS('data/DV.rds'))

mx <- split(mx, mx$cluster)

MX <- lapply(mx, function(x) data.table(dcast(x[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV')))
MX <- lapply(MX, function(x) as.data.frame(x[, year := NULL]))

u <- sapply(as.data.frame(t(sapply(MX, function(y) {
  
  lmom.atsite <- as.data.frame(t(apply(y, 2, function(x) samlmu(x))))
  at.site.para <- t(apply(lmom.atsite, 1, pelgpa))
  
  site.para <- rbindlist(lapply(1:dim(y)[2], function(j) {
    
    l <- sapply(lapply(1:500, function(i) rgpa(length(which(!is.na(y[,j]))),
                                               at.site.para[j, 'xi'], 
                                               at.site.para[j, 'alpha'],
                                               at.site.para[j, 'k'])), samlmu)
    p <- apply(l, 2, gpa.para)
    p <- melt(as.data.table(rbind(p, l), keep.rownames = T), id.vars = 1)
    
    p
    
  }), idcol = 'id')
  
  site.para <- dcast(site.para, id + variable ~ rn, value.var = 'value')
  site.para[, gamma := alpha/xi]
  
  site.q <- site.para[, lapply(.SD, function(x) quantile(x, c(.25, .75))), by = id, .SDcols = c('gamma', 'k')]
  site.q[, id := rep(c('25%', '75%'), length.out = .N)]
  site.mq <- t(site.q[, lapply(.SD, mean), by = id, .SDcols = c('gamma', 'k')])[-1,]
  
  site.iqr <- as.numeric(site.mq[, 2]) - as.numeric(site.mq[, 1])
  
  dta.fit <- sim(y, dist = 'gpa', trim = c(0,0))
  
  f <- fit(sample.sim(dta.fit, length = 500,  type = 'nonpar'), dta.fit)
  
  reg.para <- as.data.table(t(sapply(1:500, function(i) f[[i]]$REG)))
  reg.para[, gamma := alpha/xi]
  
  reg.iqr <- reg.para[, lapply(.SD, IQR), .SDcols = c('alpha', 'k')]
  
  # site.iqr
  # reg.iqr
  
  uncer.para <- reg.iqr/site.iqr
  
  ## Nlet 2, 50
  
  site.N <- site.para[, `:=`(N2 = qgpa(1 - 1/2, xi, alpha, k)/l_1,
                             N50 = qgpa(1 - 1/50, xi, alpha, k)/l_1), by = .(id, variable)]
  site.Nq <- site.N[, lapply(.SD, function(x) quantile(x, c(.25, .75))), by = id, .SDcols = c('N2', 'N50')]
  site.Nq[, id := rep(c('25%', '75%'), length.out = .N)]
  site.mNq <- t(site.Nq[, lapply(.SD, mean), by = id, .SDcols = c('N2', 'N50')])[-1,]
  
  site.Niqr <- as.numeric(site.mNq[, 2]) - as.numeric(site.mNq[, 1])
  
  reg.N <- reg.para[, `:=`(N2 = qgpa(1 - 1/2, xi, alpha, k),
                           N50 = qgpa(1 - 1/50, xi, alpha, k))] #########################################
  reg.Niqr <- reg.N[, lapply(.SD, IQR), .SDcols = c('N2', 'N50')]  
  
  uncer.N <- reg.Niqr/site.Niqr
  
  c(uncer.para, uncer.N)
  
}))), unlist)

rownames(u) <- paste('Cluster', rownames(u))

u <- round(100*(1 - u), 2)

# saveRDS(u, 'paper/data/uncer.rds')
