lop <- c("data.table", "ggplot2", "fst", "stats4", 
         "parallel", "nloptr", "extraDistr", "lmom")

to.instal <- lop[which(x = !lop %in% installed.packages()[,"Package"])]

if(length(to.instal) != 0) install.packages(to.instal)

temp <- lapply(X = lop, 
               FUN = library,
               character.only = T)
rm(temp)

## cesty a vlakna ####
cr_number <- detectCores() - 2

pth <- getwd()

source(file = "./R/aux_fun.R")

## n samplu ####
number_of_runs <- 1
n_smp <- 3

## mpara ####

loc_m <- seq(from = 0,
             to = 100, 
             length.out = 5)
scale_m <- seq(from = 1,
               to = 10, 
               length.out = length(x = loc_m))
shape_m <- seq(from = -2,
               to = 2, 
               length.out = length(x = loc_m) * 2)

## pocitani ####

tm <- system.time(
  expr = {
    for (dist_smp in c("gev", "gpd")) {
      
      for (dist in c("gev", "gpd")) {
        
        temp_out <- list()
        
        for (loc_var in loc_m) {
          
          for (scale_var in scale_m) {
            
            for (shape_var in shape_m) {
              
              loc <- loc_var
              scale_p <- scale_var
              shape_p <- shape_var
              
              mc <- mclapply(
                X = 1:number_of_runs,
                FUN = function(i) {
                  
                  x <- lapply(
                    X = c(5, 10, 20, 100, 500, 1000),
                    FUN = function (n) {
                      
                      par_orig <- c(mu = loc,
                                    sigma = scale_p,
                                    xi = shape_p)
                      
                      mx <- data.table(mx = do.call(what = paste0("r", dist_smp),
                                                    args = c(list(n = n),
                                                             par_orig)))
                      
                      opt <- function(loc, scale, shape) {
                        
                        aux <- do.call(
                          what = paste0("d", dist),
                          args = list(x = mx[, mx],
                                      mu = loc,
                                      sigma = scale,
                                      xi = shape,
                                      log = TRUE)
                        )
                        
                        -sum(x = aux[!is.infinite(aux)],
                             na.rm = TRUE)
                      }
                      
                      e <- try(
                        expr = {
                          
                          ml <- mle(minuslogl = opt,
                                    method = "L-BFGS-B",
                                    start = list(loc = mean(x = mx[, mx]),
                                                 scale = sd(x = mx[, mx]),
                                                 shape = 0),
                                    nobs = as.integer(x = n),
                                    lower = c(-Inf, .01, -(shape_p * 10)),
                                    upper = c(Inf, Inf, (shape_p * 10)))
                        },
                        silent = TRUE
                      )
                      
                      if (inherits(x = e,
                                   what = "try-error")) {
                        
                        par_ml <- p_val_ml <-  NA
                      } else {
                        
                        par_ml <- coef(ml)
                        
                        if (!between(x = par_ml[3],
                                     lower = -(shape_p * 7),
                                     upper = (shape_p * 7))) {
                          
                          par_ml <- p_val_ml <- NA
                        } else {
                          
                          base_ad_ml <- AD_test(data = mx[, mx],
                                                dist = dist,
                                                mu = par_ml["loc"],
                                                sigma = par_ml["scale"],
                                                xi = par_ml["shape"])
                          
                          pars_ml <- list(mu = par_ml["loc"],
                                          sigma = par_ml["scale"],
                                          xi = par_ml["shape"])
                          
                          smp_ml <- replicate(n = n_smp,
                                              expr = do.call(what = paste0("r", dist),
                                                             args = c(list(n = n),
                                                                      pars_ml)))
                          
                          smp_ad_ml <- apply(X = smp_ml,
                                             MARGIN = 2,
                                             FUN = function(x) {
                                               
                                               e <- try(
                                                 expr = {
                                                   
                                                   ml <- mle(minuslogl = opt,
                                                             method = "L-BFGS-B",
                                                             start = list(loc = mean(x = mx[, mx]),
                                                                          scale = sd(x = mx[, mx]),
                                                                          shape = 0),
                                                             nobs = as.integer(x = n),
                                                             lower = c(-Inf, .01, -(shape_p * 10)),
                                                             upper = c(Inf, Inf, (shape_p * 10)))
                                                 },
                                                 silent = TRUE
                                               )
                                               
                                               if (inherits(x = e,
                                                            what = "try-error")) {
                                                 
                                                 par_aux <- NA
                                               } else {
                                                 
                                                 par_aux <- coef(object = ml)
                                               }
                                               
                                               par_aux <- unlist(x = par_aux)
                                               
                                               AD_test(data = x,
                                                       dist = dist,
                                                       mu = par_aux["loc"],
                                                       sigma = par_aux["scale"],
                                                       xi = par_aux["shape"])
                                             }
                          )
                          
                          p_val_ml <- length(x = which(x = smp_ad_ml >= base_ad_ml)) / n_smp
                        }
                      }
                      
                      
                      
                      lmom <- samlmu(x = mx[, mx])
                      
                      par_l <- do.call(what = paste0("pel", ifelse(test = dist == "gpd",
                                                                   yes = "gpa",
                                                                   no = dist)),
                                       args = list(l = lmom))
                      names(x = par_l) <- c("loc", "scale", "shape")
                      
                      par_l[3] <- par_l[3] * -1
                      
                      if (!between(x = par_l[3],
                                   lower = -(shape_p * 7),
                                   upper = (shape_p * 7))) {
                        
                        par_l <- p_val_l <- NA
                      } else {
                        
                        base_ad_l <- AD_test(data = mx[, mx],
                                             dist = dist,
                                             mu = par_l["loc"],
                                             sigma = par_l["scale"],
                                             xi = par_l["shape"])
                        
                        pars_l <- list(mu = par_l["loc"],
                                       sigma = par_l["scale"],
                                       xi = par_l["shape"])
                        
                        smp_l <- replicate(n = n_smp,
                                           expr = do.call(what = paste0("r", dist),
                                                          args = c(list(n = n),
                                                                   pars_l)))
                        
                        smp_ad_l <- apply(X = smp_l,
                                          MARGIN = 2,
                                          FUN = function(x) {
                                            
                                            par_aux <- do.call(what = paste0("pel", ifelse(test = dist == "gpd",
                                                                                           yes = "gpa",
                                                                                           no = dist)),
                                                               args = list(l = samlmu(x = x)))
                                            names(x = par_aux) <- c("mu", "sigma", "xi")
                                            
                                            par_aux[3] <- par_aux[3] * -1
                                            
                                            par_aux <- unlist(x = par_aux)
                                            
                                            AD_test(data = x,
                                                    dist = dist,
                                                    mu = par_aux["mu"],
                                                    sigma = par_aux["sigma"],
                                                    xi = par_aux["xi"])
                                          }
                        )
                        
                        p_val_l <- length(x = which(x = smp_ad_l >= base_ad_l)) / n_smp
                      }
                      
                      par_n <- fitDist(data = mx[,mx],
                                       dist = dist,
                                       n.points = 200,
                                       norm = "N2",
                                       constrain = FALSE,
                                       pp = "chegodaev")
                      
                      par_n <- unlist(x = par_n)
                      names(x = par_n) <- c("loc", "scale", "shape")
                      
                      if (!between(x = par_n[3],
                                   lower = -(shape_p * 7),
                                   upper = (shape_p * 7))) {
                        
                        par_n <- p_val_n <- NA
                      } else {
                        
                        base_ad_n <- AD_test(data = mx[, mx],
                                             dist = dist,
                                             mu = par_n["loc"],
                                             sigma = par_n["scale"],
                                             xi = par_n["shape"])
                        
                        pars_n <- list(mu = par_n["loc"],
                                       sigma = par_n["scale"],
                                       xi = par_n["shape"])
                        
                        smp_n <- replicate(n = n_smp,
                                           expr = do.call(what = paste0("r", dist),
                                                          args = c(list(n = n),
                                                                   pars_n)))
                        
                        smp_ad_n <- apply(X = smp_n,
                                          MARGIN = 2,
                                          FUN = function(x) {
                                            
                                            par_aux <- fitDist(data = x,
                                                               dist = dist,
                                                               n.points = 200,
                                                               norm = "N2",
                                                               constrain = FALSE,
                                                               pp = "chegodaev")
                                            
                                            par_aux <- unlist(x = par_aux)
                                            
                                            AD_test(data = x,
                                                    dist = dist,
                                                    mu = par_aux["mu"],
                                                    sigma = par_aux["sigma"],
                                                    xi = par_aux["xi"])
                                          }
                        )
                        
                        p_val_n <- length(x = which(x = smp_ad_n >= base_ad_n)) / n_smp
                      }
                      
                      
                      
                      out <- data.table(par = c("loc", "scale", "shape"),
                                        par_orig,
                                        par_ml,
                                        par_l,
                                        par_n,
                                        p_val_ml,
                                        p_val_l,
                                        p_val_n,
                                        n = n)
                      
                      out
                    }
                  )
                  
                  out_mc <- rbindlist(l = x)
                  
                  out_mc
                },
                mc.cores = cr_number
              )
              
              temp_out <- c(temp_out, rbindlist(l = mc))
            }
          }
        }
        
        out <- rbindlist(l = temp_out) 
        
        write_fst(x = out, 
                  path = paste0("./data/fitting/", toupper(dist_smp),
                                "_fitted_to_", toupper(dist),
                                ".fst"))
      }
    }
  }
)

tm
