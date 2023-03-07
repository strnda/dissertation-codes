library(data.table); library(ggplot2); library(lmomRFA)

source('R/auxiliary_functions/imports.R')

cluster.number <- 1:3

mx <- data.table(readRDS('data/DV.rds'))

for (i in cluster.number) {
  
  mxx <- mx[which(cluster %in% cluster.number[i]),]
  
  MX <- data.table(dcast(mxx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
  MX <- MX[, year := NULL]
  
  print(mean(apply(MX, 2, function(x) {length(which(!is.na(x)))})))
  
  xxx <- regtst(regsamlmu(MX)) 
  print(length(which(xxx$D > xxx$Dcrit)))
}
