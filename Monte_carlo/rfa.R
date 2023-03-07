library(lmomRFA); library(mvtnorm); library(nim); library(extraDistr); library(fst)

dir.create(path = "./data/reg/", 
           showWarnings = FALSE, 
           recursive = TRUE)

source(file = "./R/aux_fun.R")

n_smp <- 3000

#### ####

for (p0 in seq(from = 0,
               to = .9,
               by = .1)) {
  
  data(precip_max)
  
  precip_max[, -1] <- precip_max[, -1] * rbinom(n = ncol(x = precip_max[, -1]) * nrow(x = precip_max[, -1]), 
                                                size = 1, 
                                                prob = 1 - p0)
  
  precip_max <- as.data.table(x = t(x = apply(X = precip_max, 
                                              MARGIN = 1, function(x) {
                                                x[which(x = x == 0)] <- NA
                                                x
                                              }
  )))
  
  r_mom <- regsamlmu(x = precip_max[, -1])
  
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
  
  orig_fit <- as.data.table(x = c(r_fit$para, dif))
  
  s <- lapply(X = 1:3, 
              FUN = function(i) { 
                smp(x = precip_max,
                    n = n_smp,
                    method = i)
              }
  )
  smp_fit <- lapply(X = s, 
                    fitsmp)
  names(x = smp_fit) <- paste0("method_", 1:3)
  
  smp_fit <- rbindlist(l = smp_fit,
                       idcol = "method")
  
  write_fst(x = smp_fit, 
            path = paste0("./data/reg/smp_", p0, ".fst"))
  
  write_fst(x = orig_fit, 
            path = paste0("./data/reg/orig_", p0, ".fst"))
}

read_fst(path = "./data/reg/smp_0.8.fst")
