library(data.table); library(ggplot2); library(lmom)

source('R/auxiliary_functions/imports.R')

## select cluster number
cluster.number <- 1

{

## data import
mx <- data.table(readRDS('data/DV.rds'))

mx <- mx[which(cluster %in% cluster.number),]

MX <- data.table(dcast(mx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
MX <- MX[, year := NULL]

# lmom.atsite <- as.data.frame(t(apply(MX, 2, function(x) samlmu(x, trim = trim))))
# at.site.para <- t(apply(lmom.atsite, 1, pelgpa))

}

# data('precip_max', package = 'nim')
# 
# MX <- precip_max[,-1]

# x <- regtst(regsamlmu(MX))
# 
# MX <- MX[which(x$D < min(x$Dcrit))]

## fit stationary model
dta.fit <- sim(MX, dist = 'gpa', trim = c(0,0)) ################ trim

MX <- synth.data(dta.fit, use.NA = F, dist = 'gpa', dependent = T, sample_para = F)

dta.fit <- sim(MX, dist = 'gpa', trim = c(0,0)) ################ trim

gumbelplot(dta.fit, maxima_with_NA = F)
# qq(dta.fit, maxima_with_NA = F)

s <- sample(dta.fit, length = 100,  type = 'nonpar')
f <- fit(s, dta.fit)

growthcurve(dta.fit, f)
para.boxplot(dta.fit, f)

s <- sample(dta.fit, length = 100,  type = 'zero')
f <- fit(s, dta.fit)

growthcurve(dta.fit, f)
para.boxplot(dta.fit, f)

s <- sample(dta.fit, length = 100,  type = 'para_cor', na = T)
f <- fit(s, dta.fit)

growthcurve(dta.fit, f)
para.boxplot(dta.fit, f)

