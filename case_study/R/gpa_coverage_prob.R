library(data.table); library(ggplot2); library(lmomRFA)

source('R/auxiliary_functions/imports.R')

cluster.number <- 1 # syth data template

use.NA <- c(TRUE, FALSE) # include NA values (years w/o drought)
dependency <- c('dep', 'indep') # dependace between sites
method <- c('zero') # sampling methods

N <- 500 # number of samples

{
  mx <- data.table(readRDS('data/DV.rds'))
  
  mx <- mx[which(cluster %in% cluster.number),]
  
  MX <- data.table(dcast(mx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
  MX <- MX[, year := NULL]
}

dta.fit <- sim(MX, dist = 'gpa', trim = c(0,0))

###### patek 16. v 11 ######

for (k in seq_along(dependency)) {
  
  for(j in seq_along(use.NA)) {
    
    # MX <- synth.data(dta.fit, dependent = T, use.NA = use.NA[j])
    # 
    # dta.synth <- sim(MX, dist = 'gpa', trim = c(0,0))
    # 
    # saveRDS(dta.synth, paste0('samples/synth/', dependency[k], '/with_NA_', use.NA[j], '/dta.RDS'))
    
    dta.synth <- readRDS(paste0('samples/synth/', dependency[k], '/with_NA_', use.NA[j], '/dta.RDS'))
    
    for(i in seq_along(method)) {
      
      s <- sample(dta.synth, length = N,  type = method[i])
      f <- fit(s, dta.synth)
      
      rs <- lapply(f, function(x) sample.sim(x, length = N,  type = method[i]))
      
      print(object.size(rs), units = 'Gb')
      
      saveRDS(rs, paste0('samples/synth/', dependency[k], '/with_NA_', use.NA[j], '/', method[i], '.RDS'))
      
      gc()
    }
  }
}

# files <- list.files('samples/synth/', recursive = T, pattern = '.RDS')

for (k in seq_along(dependency)) {
  for(j in seq_along(use.NA)) { 
    for(i in seq_along(method)) {
    
      dta.synth <- readRDS(paste0('samples/synth/', dependency[k], '/with_NA_', use.NA[j], '/dta.RDS'))
      rs <- readRDS(paste0('samples/synth/', dependency[k], '/with_NA_', use.NA[j], '/', method[i], '.RDS'))
      
      para <- lapply(rs, function(x){
        f <- fit(x, dta.synth)
        dt <- suppressWarnings(melt(as.data.table(t(sapply(f, function(y) y$REG)))))
      })
      
      para <- rbindlist(para, idcol = T)
      para[, orig := rep(dta.synth$REG, each = max(as.numeric(.id)))]
      para[, mean := mean(value), by = variable]
      para[, eval := (quantile(value, .05) <= orig) && (quantile(value, .95) >= orig), by = .(.id, variable)]
      para[, eval.proc := length(which(eval))/.N, by = variable]
      para[, out := paste(variable, '- coverage probability =', eval.proc)]
      
      saveRDS(para, paste0('samples/synth/', dependency[k], '/with_NA_', use.NA[j], '/', method[i], '_resampled_parameters.RDS'))
      
      cp <- ggplot(para) +
        geom_boxplot(aes(x = factor(.id), y = value), fill = NA, outlier.shape = 21, outlier.size = .75, show.legend = F) +
        geom_hline(aes(yintercept = mean), colour = 'red4', linetype = 2) +
        geom_hline(aes(yintercept = orig), colour = 'red4') +
        theme_bw() +
        facet_wrap(~out, scales = 'free', ncol = 1) +
        labs(y = 'Value', subtitle = paste('Data type -', dependency[k], 
                                           '\nUsed NA -', use.NA[j],
                                           '\nSampling method -', method[i])) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
      
      ggsave(paste0('samples/cp_gg/',dependency[k], '_with_NA_', use.NA[j], '_method_', method[i],'.pdf'), cp, device = 'pdf', width = 20, height = 7.5)
    
      gc()
    }
  }
}
