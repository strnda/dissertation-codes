library(data.table); library(ggplot2); library(lmom)

source('R/auxiliary_functions/imports.R')

cluster.number <- 1
a <- .10
N <- 3000

{
  
  mx <- data.table(readRDS('data/DV.rds'))
  
  mx <- mx[which(cluster %in% cluster.number),]
  
  MX <- data.table(dcast(mx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
  MX <- MX[, year := NULL]
  # MX <- MX[, '324000' := NULL]
  
  dta.fit <- sim(MX, dist = 'gpa', trim = c(0,0)) 
  
  y <- dta.fit$data
  
  para.mx <- as.data.table(merge(melt(y), 
                                 data.table(SP_ID = names(MX), t(apply(y, 2, function(x) pelgpa(samlmu(x))))), 
                                 by.x = 'variable',
                                 by.y = 'SP_ID'))
  
  ad.base <- para.mx[, .(base_ad = AD.test(val = value,
                                           location = unique(xi), 
                                           scale = unique(alpha), 
                                           shape = unique(k), 
                                           dist = 'gpa')),
                     by = variable]
  
  gp <- AD.regional(dta.fit, save.samples = 's', N = N)
  
  ad.base[base_ad > gp]
  
  f <- fit(s, dta.fit)
  rp.res <- sapply(f, function(x) AD.regional(x, N = N))
  
}

fff <- function(x,y) {x[base_ad > y]}
fff(ad.base, rp.res[4])
