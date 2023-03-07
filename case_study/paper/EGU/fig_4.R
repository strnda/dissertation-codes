library(data.table); library(ggplot2); library(lmomRFA); library(gridExtra); library(grid)

source('R/auxiliary_functions/imports.R')

cluster.number <- 1:3

mx <- data.table(readRDS('data/DV.rds'))
disc <- readRDS('paper/data/disc.RDS')

for (i in cluster.number) {
  
  mxx <- mx[which(cluster == cluster.number[i]),]
  
  MX <- data.table(dcast(mxx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
  MX <- MX[, year := NULL]
  
  moments <- data.table(t(apply(MX, 2, samlmu)), keep.rownames = T)
  
  names(moments)[2:5] <- paste('m', 1:4, sep = '_')
  
  model <- sim(MX, dist = 'gpa', trim = c(0,0)) 
  
  dta <- as.data.table(model$data)
  para <- model$REG
  scaling.factor <- model$scaling_factor
  
  res.gp <- suppressWarnings(melt(dta))
  res.gp <- data.table(variable = names(dta), sf = scaling.factor)[res.gp, on = c('variable')]
  # res.gp <- res.gp[, p0 := length(which(!is.na(value)))/.N, by = variable]
  # res.gp <- res.gp[, p := p0 + (1 - p0)*(rank(value) - .3)/(length(value) + .4), by = variable]
  res.gp <- res.gp[!is.na(value), p := (rank(value) - .3)/(length(value) + .4), by = variable]
  res.gp <- res.gp[, gumbel.variate := -log(-log(as.numeric(p)))]
  res.gp <- res.gp[, scaled.value := value/sf]
  
  p <- seq(min(res.gp[,p], na.rm = T), max(res.gp[,p], na.rm = T), .001)
  
  regional <- data.table(x = -log(-log(p)), y = qgpa(p, xi = para[1], alpha = para[2], k = para[3]))
  
  (gp <- ggplot() +
      geom_line(data = res.gp[!variable %in% unlist(disc),], aes(x = gumbel.variate, y = scaled.value, group = variable), colour = 'grey25', alpha = .5, na.rm = T) +
      geom_point(data = res.gp[!variable %in% unlist(disc),], aes(x = gumbel.variate, y = scaled.value, group = variable), colour = 'grey25', fill = 'grey50', alpha = .5, shape = 21, na.rm = T) +
      geom_line(data = res.gp[variable %in% unlist(disc),], aes(x = gumbel.variate, y = scaled.value, group = variable), lty = 'dashed', colour = 'red4', na.rm = T) +
      geom_point(data = res.gp[variable %in% unlist(disc),], aes(x = gumbel.variate, y = scaled.value, group = variable), colour = 'red4', fill = 'red', shape = 21, na.rm = T) +
      geom_line(data = regional, aes(x = x, y = y), col = 'grey5', lwd = 1) + 
      theme_bw() +
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
      labs(x = expression(-log(-log(italic(F)))), y = 'Scaled value', title = paste('Cluster', i), subtitle = 'Gumbel plot') +
      theme(aspect.ratio = .5))
  
  (discgg <- ggplot(data = NULL, aes(x = m_3, y = m_2/m_1)) +
      geom_point(data = moments, colour = 'grey5', fill = 'grey45', shape = 21) +
      geom_point(data = moments[rn %in% unlist(disc),], colour = 'red4', shape = 13, size = 7.5) +
      theme_bw() + 
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
      labs(x = 'L-skewness', y = 'L-CV', subtitle = 'Discordancy') +
      theme(aspect.ratio = 1))
  
  cl <- grid.arrange(gp, discgg, nrow = 1, widths = c(1.75, .95))
  
  ggsave(paste0('~/Desktop/cl',i,'.pdf'), cl)
}


