# site.para[, .(Q25 = quantile(k, .25), Q75 = quantile(k, .75)), by = id][, mean(Q75) - mean(Q25)]/reg.para[, IQR(k)]

library(data.table); library(ggplot2); library(lmom); library(bilan)

source('R/auxiliary_functions/imports.R')

mx <- data.table(readRDS('data/DV.rds'))

mx <- split(mx, mx$cluster)

MX <- lapply(mx, function(x) data.table(dcast(x[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV')))
MX <- lapply(MX, function(x) as.data.frame(x[, year := NULL]))

range(sapply(as.data.frame(t(sapply(MX, function(y) {lmom.atsite <- as.data.frame(t(apply(y, 2, function(x) samlmu(x))))}))), unlist)[, 'l_1'])

sapply(MX, dim)


cl <- readRDS('~/ownCloud/Active Docs/nlet/data/cluster_numbers_final.rds')
dv <- readRDS('~/ownCloud/Active Docs/nlet/data/de_vol_p3_0.2.rds')
dta <- merge(dv, cl, by.x = 'SP_ID', by.y = 'DBCN')

xtable::xtable(round(cbind(dta[, .(P = mean(P, na.rm = TRUE)*12 , DV = mean(DV, na.rm = TRUE)), by = cluster],
                           data.table(p = sapply(MX, function(x) length(which(is.na(x)))/length((is.na(x)))))), 2))
