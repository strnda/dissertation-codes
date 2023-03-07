library(readxl); library(xtable)

## tab. 1

DV_1 <- read_excel('paper/data/DV_1.xlsx')
names(DV_1)[1] <- 'Period'

xtable(DV_1, caption = 'X')

## tab. 2

DV_val <- read_excel('paper/data/DV_val.xlsx')

xtable(DV_val, caption = 'X')

## tab. 3

`Cluster 1` <- c(xi = -0.01, alpha = 0.86, kappa = -0.15, ADcrit = 2.42) 
`Cluster 2` <- c(xi = -0.02, alpha = 0.83, kappa = -0.19, ADcrit = 2.64) 
`Cluster 3` <- c(xi = -0.04, alpha = 0.71, kappa = -0.32, ADcrit = 2.79) 

xtable(rbind(`Cluster 1`, `Cluster 2`, `Cluster 3`), caption = 'X')

## tab. 4

u <- readRDS('paper/data/uncer.rds')

xtable(u)
