library(stats4); library(parallel); library(nloptr); library(extraDistr); library(lmom)
library(data.table); library(ggplot2); library(fst)

source(file = "./R/aux_fun.R")

dist_smp <- "gpd"
dist <- "gpd"
# 
# dta <- read_fst(path = paste0("./data/fitting/", toupper(dist_smp),
#                               "_fitted_to_", toupper(dist),
#                               ".fst"), 
#                 as.data.table = TRUE)

ifelse(
  test = dist_smp == "gpd",
  yes = {
    loc <- 1
    scale_p <- .5
    shape_p <- .2}, 
  no = {
    loc <- 10
    scale_p <- 1.5
    shape_p <- .1
  }
)

n <- 100

mx <- data.table(mx = do.call(what = paste0("r", dist), 
                              args = list(n = n,
                                          mu = loc,
                                          sigma = scale_p,
                                          xi = shape_p)))

par <- fitDist(data = mx[,mx],
               dist = dist,
               n.points = 200,
               norm = "N2",
               constrain = FALSE,
               pp = "chegodaev")
par <- unlist(par)

base_ad <- AD_test(data = mx[, mx],
                   dist = "gev",
                   mu = par["mu"],
                   sigma = par["sigma"],
                   xi = par["xi"])

base_ad

n_smp <- 200

pars <- list(mu = par["mu"],
             sigma = par["sigma"],
             xi = par["xi"])

smp <- replicate(n = n_smp,
                 expr = do.call(what = paste0("r", dist),
                                args = c(list(n = dim(x = mx)[1]),
                                         pars)))

smp_ad <- apply(X = smp,
                MARGIN = 2,
                FUN = function(x) {
                  
                  par_aux <- fitDist(data = x,
                                     dist = dist,
                                     n.points = 200,
                                     norm = "N2",
                                     constrain = FALSE,
                                     pp = "chegodaev")
                  
                  par_aux <- unlist(x = par_aux)

                  AD_test(data = mx[, mx],
                          dist = "gev",
                          mu = par_aux["mu"],
                          sigma = par_aux["sigma"],
                          xi = par_aux["xi"])
                }
)

p_val <- length(x = which(x = smp_ad >= base_ad)) / n_smp

p_val
