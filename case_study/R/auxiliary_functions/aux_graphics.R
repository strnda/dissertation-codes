gumbelplot <- function(model_object, dist = NULL, title = 'Gumbel plot', maxima_with_NA = FALSE, method = if ('ggplot2' %in% installed.packages()[,'Package']) {'ggplot'} else {'base'}) {
  
  ## res (init res will need to be rewritten for nim)
  
  ## extract information from the model (data, parameters, scaling factors and distribution)
  dta <- as.data.table(model_object$data)
  para <- model_object$REG
  scaling.factor <- model_object$scaling_factor
  
  cl <- attr(model_object, 'sim.call')
  
  ifelse(all(!is.na(as.numeric(as.character(cl$trim)[2:3]))), 
         assign('trim', as.numeric(as.character(cl$trim)[2:3])),
         assign('trim', as.numeric(as.character(formals(sim)$trim)[2:3])))
  
  if(is.null(dist)) {dist <- attr(model_object, 'sim.call')$dist}
  
  # dta <- data.table(attributes(n)$data)
  
  res.gp <- suppressWarnings(melt(dta))
  res.gp <- data.table(variable = names(dta), sf = scaling.factor)[res.gp, on = c('variable')]
  res.gp <- res.gp[res.gp[, .I[which(!(value %in% n.min(value, trim[1])) & !(value %in% n.max(value, trim[2])))], by = variable]$V1]
  
  if(maxima_with_NA) {
    
    res.gp <- res.gp[, p0 := length(which(!is.na(value)))/.N, by = variable]
    res.gp <- res.gp[, p := p0 + (1 - p0)*(rank(value) - .3)/(length(value) + .4), by = variable]
  } else {
    
    res.gp <- res.gp[!is.na(value), p := (rank(value) - .3)/(length(value) + .4), by = variable]
  }

  res.gp <- res.gp[, gumbel.variate := -log(-log(as.numeric(p)))]
  res.gp <- res.gp[, scaled.value := value/sf]
  
  p <- seq(min(res.gp$p, na.rm = T), max(res.gp$p, na.rm = T), .001)
  
  regional <- data.table(x = -log(-log(p)), y = do.call(paste0('q', dist), list(p, xi = para[1], alpha = para[2], k = para[3])))
  
  # graphics
  
  if(method == 'base') {
    
    sres.gp <- split(res.gp[!is.na(res.gp$scaled.value),], res.gp$variable[!is.na(res.gp$scaled.value)])
    
    plot(NULL,
         xlim = c(min(regional$x), max(regional$x)*1.15),
         ylim = c(min(res.gp$scaled.value, na.rm = T), max(res.gp$scaled.value, na.rm = T)),
         bty = 'l',
         xlab = expression(-log(-log(p))),
         ylab = 'Value',
         main = title)
    
    grid()
    
    lapply(sres.gp, function(x) {
      points(sort(x$gumbel.variate),
             sort(x$scaled.value),
             pch = 21,
             col = 'grey15', 
             bg = '#36648b50',
             cex = .75)
      lines(sort(x$gumbel.variate),
            sort(x$scaled.value),
            col = '#36648b50')
    })
    
    lines(regional,
          type = 'l',
          col = 'red4',
          lwd = .75)
  }
  
  if(method %in% c('ggplot', 'plotly')) {
    
    gp <- ggplot2::ggplot(res.gp) +
      ggplot2::geom_line(ggplot2::aes(x = gumbel.variate, y = scaled.value, group = variable), colour = 'steelblue4', alpha = .5, na.rm = T) +
      ggplot2::geom_point(ggplot2::aes(x = gumbel.variate, y = scaled.value, group = variable), colour = 'grey15', fill = 'steelblue4', alpha = .5, shape = 21, na.rm = T) +
      ggplot2::geom_line(data = regional, ggplot2::aes(x = x, y = y), col = 'red4', lwd = .75) + 
      ggplot2::theme_bw() +
      ggplot2::labs(x = '-log(-log(p))', y = 'Value', title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = 'black'))
    
    if(method == 'plotly') {
      
      gp <- plotly::ggplotly(gp)
    }
    
    return(gp)
  }
}

growthcurve <- function (model_object, fitted_bootstrap, dist = NULL, outer_ribbon = c(0.05, 0.95), inner_ribbon = c(0.25, 0.75), rp = T, return.period = c(5, 10, 20, 50, 100), method = if ('ggplot2' %in% installed.packages()[,'Package']) {'ggplot'} else {'base'}, subtitle = T) {
  
  # res (init res will need to be rewritten for nim)
  
  prbs <- sort(c(outer_ribbon, inner_ribbon))
  para <- model_object$REG
  
  if(is.null(dist)) {dist <- attr(model_object, 'sim.call')$dist}
  
  qs <- seq(.01, 1 - 1/max(return.period)*.5, 1/max(return.period))
  qaux <- data.table(rbindlist(lapply(fitted_bootstrap, function(x) {data.frame(q = do.call(paste0('q', dist), list(qs, xi = x$REG[1], alpha = x$REG[2], k = x$REG[3])))}),
                               idcol = 'sample'), ##########################################################
                     probs = seq_along(qs))
  q <- qaux[, .(val = quantile(q, c(.5, prbs)),
                q = c('sample_median', 'rib_1_min', 'rib_2_min', 'rib_2_max', 'rib_1_max')), 
            by = probs]
  res.gc <- cbind(dcast(q, probs ~ q, value.var = 'val'),
                  data.table(gumbel.variate = -log(-log(qs)),
                             scaled.value = do.call(paste0('q',dist), list(qs, xi = para[1], alpha = para[2], k = para[3]))))
  
  # graphics
  
  if(method == 'base') {
    
    plot(NULL,
         xlim = c(min(res.gc$gumbel.variate), max(res.gc$gumbel.variate)),
         ylim = c(min(res.gc[,c('rib_1_min', 'rib_2_min', 'rib_2_max', 'rib_1_max')]), 
                  max(res.gc[,c('rib_1_min', 'rib_2_min', 'rib_2_max', 'rib_1_max')])),
         bty = 'l',
         xlab = expression(-log(-log(p))),
         ylab = 'Value',
         main = if(subtitle) bquote(atop(Growthcurve, atop(italic(.(attributes(fitted_bootstrap)$type))),'')) else 'Growth curve')
    
    grid()
    
    polygon(c(res.gc$gumbel.variate, rev(res.gc$gumbel.variate)), 
            c(res.gc$rib_1_max, rev(res.gc$rib_1_min)),
            col = '#36648b40', border = NA)
    
    polygon(c(res.gc$gumbel.variate, rev(res.gc$gumbel.variate)), 
            c(res.gc$rib_2_max, rev(res.gc$rib_2_min)),
            col = '#36648b80', border = NA)
    
    lines(res.gc$gumbel.variate,
          res.gc$sample_median,
          type = 'l',
          lty = 2,
          col = 'white',
          lwd = .5)
    
    lines(res.gc$gumbel.variate,
          res.gc$scaled.value,
          type = 'l',
          col = 'red4',
          lwd = .75)
    
    if(rp) {
      
      axis.lim <- par('usr')
      
      rp.lab <- return.period
      rp.x <- -log(-log(1 - 1/rp.lab))
      rp.y <- axis.lim[3] + (axis.lim[4] - axis.lim[3])*.05
      
      axis(side = 3, at = rp.x, pos = rp.y, labels = rp.lab)
      
      text(mean(rp.x[rev(rank(rp.lab))[1:2]]), rp.y + par('cxy')[2], 'Return period', adj = c(.75, -2.75))
    }
  }
    
  if(method %in% c('ggplot', 'plotly')) {
    
    gc <- ggplot2::ggplot(res.gc) +
      ggplot2::geom_ribbon(ggplot2::aes(x = gumbel.variate, ymin = rib_1_min, ymax = rib_1_max), fill = 'steelblue4', alpha = .4) +
      ggplot2::geom_ribbon(ggplot2::aes(x = gumbel.variate, ymin = rib_2_min, ymax = rib_2_max), fill = 'steelblue4', alpha = .8) +
      ggplot2::geom_line(ggplot2::aes(x = gumbel.variate, y = sample_median), col = 'white', lwd = .5, linetype = 2) + 
      ggplot2::geom_line(ggplot2::aes(x = gumbel.variate, y = scaled.value), col = 'red4', lwd = .75) + 
      ggplot2::theme_bw() +
      ggplot2::labs(x = '-log(-log(p))', 
                    y = 'Value', 
                    title = 'Growth curve',
                    subtitle = if(subtitle) attributes(fitted_bootstrap)$type else NULL) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5),
                     plot.subtitle = ggplot2::element_text(hjust = .5),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = 'black'))
    if(rp) {
      
      axis.lim <- c(-log(-log(range(qs))), range(qaux$q))
      
      rp.lab <- return.period
      rp.x <- -log(-log(1 - 1/rp.lab))
      rp.y <- axis.lim[3] + (axis.lim[4] - axis.lim[3])*.05
      
      rp.dta <- data.table(rp.x, rp.y, rp.lab)
      
      gc <- gc + ggplot2::geom_point(data = rp.dta, ggplot2::aes(x = rp.x, y = rp.y), shape = '|', size = 3) +
        ggplot2::geom_line(data = rp.dta, ggplot2::aes(x = rp.x, y = rp.y)) +
        ggplot2::geom_text(data = rp.dta, ggplot2::aes(x = rp.x, y = rp.y*2, label = rp.lab)) +
        ggplot2::geom_text(data = rp.dta, ggplot2::aes(x = mean(rp.x[rev(rank(rp.lab))[1:2]]), y = rp.y[1]*3.5), label = 'Return period', fontface = 1)
    }
    
    if(method == 'plotly') {
      
      gc <- plotly::ggplotly(gc)
    }
    
    return(gc)
  }
}

qq <- function(...) {UseMethod('qq')}

qq.sim <- function(model_object, dist = NULL, maxima_with_NA = FALSE,method = if ('ggplot2' %in% installed.packages()[,'Package']) {'ggplot'} else {'base'}) {
  
  dta <- as.data.table(model_object$data)
  para <- model_object$REG
  scaling.factor <- model_object$scaling_factor
  
  if(is.null(dist)) {dist <- attr(dta.fit, 'sim.call')$dist}
  
  res.qq <- suppressWarnings(melt(dta))
  res.qq <- data.table(variable = names(dta), sf = scaling.factor)[res.qq, on = c('variable')]
  
  if(maxima_with_NA) {
    
    res.qq <- res.qq[, p0 := length(which(!is.na(value)))/.N, by = variable]
    res.qq <- res.qq[, p := p0 + (1 - p0)*(rank(value) - .3)/(length(value) + .4), by = variable]
  } else {
    
    res.qq <- res.qq[!is.na(value), p := (rank(value) - .3)/(length(value) + .4), by = variable]
  }
  
  res.qq <- res.qq[, scaled.value := value/sf]
  
  
  if(method == 'base') {
    
    inipar <- par()
    
    par(pty = 's')
    
    sres.qq <- split(res.qq[!is.na(res.qq$scaled.value),], res.qq$variable[!is.na(res.qq$scaled.value)])
    
    plot(NULL, 
         xlim = c(0, max(res.qq$scaled.value, na.rm = T)*1.15),
         ylim = c(0, max(res.qq$scaled.value, na.rm = T)*1.15),
         pch = 21,
         col = 'grey15', 
         bg = '#36648b90',
         bty = 'l',
         xlab = 'theoretical',
         ylab = 'sample',
         main = 'qqplot')
    
    grid()
    
    lapply(sres.qq, function(x) {
      points(sort(x$scaled.value),
             sort(do.call(paste0('q',dist), list(x$p, xi = para[1], alpha = para[2], k = para[3]))),   #########################    
             pch = 21,
             col = 'grey15', 
             bg = '#36648b90')
    })
    
    abline(0,1, col = 'red4')
    
    suppressWarnings(par(inipar))
  }
  
  if(method %in% c('ggplot', 'plotly')) {
    
    qq <- ggplot2::ggplot(res.qq) +
      ggplot2::geom_qq(ggplot2::aes(sample = scaled.value, group = variable), geom = 'point', distribution = noquote(paste0('q', dist)), dparams = list(xi = para[1], alpha = para[2], k = para[3]), colour = 'grey15', fill = 'steelblue4', shape = 21, na.rm = T) +
      ggplot2::geom_abline(colour = ('red4')) +
      ggplot2::coord_fixed() +
      ggplot2::lims(x = c(0, max(res.qq$value/res.qq$sf, na.rm = T)),
                    y = c(0, max(res.qq$value/res.qq$sf, na.rm = T))) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = 'black'))
    
    if(method == 'plotly') {
      
      qq <- plotly::ggplotly(qq)
    }
    
    return(qq)
  }
}

# qq.simsample <- function(model_object, fitted_bootstrap, dist = NULL, ribbon.1 = c(0.05, 0.95), ribbon.2 = c(0.25, 0.75), method = if ('ggplot2' %in% installed.packages()[,'Package']) {'ggplot'} else {'base'}) {
#   
#   dta <- as.data.table(model_object$data)
#   para <- model_object$REG
#   scaling.factor <- model_object$scaling_factor
#   
#   if(is.null(dist)) {dist <- attr(dta.fit, 'sim.call')$dist}
#   
#   prbs <- sort(c(ribbon.1, ribbon.2))
#   
#   xxx <- do.call(cbind, lapply(fitted_bootstrap, function(x) x$data))
#   
#   ########################################################
#   
#   qaux <- data.table(rbindlist(lapply(fitted_bootstrap, function(x) {data.frame(q = do.call(paste0('q',dist), list(qs, x$REG)))}),
#                                idcol = 'sample'), 
#                      probs = seq_along(qs))
#   q <- qaux[, .(val = quantile(q, prbs),
#                 q = c('rib_1_min', 'rib_2_min', 'rib_2_max', 'rib_1_max')), 
#             by = probs]
#   res.gc <- cbind(dcast(q, probs ~ q, value.var = 'val'),
#                   data.table(gumbel.variate = -log(-log(qs)),
#                              scaled.value = qgpa(qs, para)))
# }

ratiodiagram <- function(moments, method = 'lmrd') {
  
  names(moments) <- paste('m', 1:4, sep = '_')
  
  if(method == 'lmrd') {
    
    num <- seq(min(moments[,3])*.5, max(moments[,3])*1.2, .01)
    mr <- data.table(m_3 = num, m_4 = num*(1 + 5*num)/(5 + num))
    
    lmrd <- ggplot(data = NULL, aes(x = m_3, y = m_4)) +
      geom_line(data = mr, colour = 'red4') +
      geom_point(data = moments, colour = 'grey15', fill = 'steelblue4', shape = 21) +
      theme_bw() +
      labs(x = 'L-skewness', y = 'L-kurtosis', title = 'GPA L-moment ratio diagram')
    
    return(lmrd)
  }
  
  if(method == 'discord') {
    
    disc <- ggplot(data = NULL, aes(x = m_3, y = m_2/m_1)) +
      geom_point(data = moments, colour = 'grey15', fill = 'steelblue4', shape = 21) +
      theme_bw() +
      labs(x = 'L-skewness', y = 'L-CV', title = 'Discordancy')
    
    return(disc)
  }
  
}

para.boxplot <- function(model_object, fitted_bootstrap) {
  
  para <- as.data.table(melt(t(sapply(fitted_bootstrap, function(x) x$REG))))
  para[, orig := rep(model_object$REG, each = length(fitted_bootstrap))]
  
  bp <- ggplot(para) +
    geom_boxplot(aes(x = Var2, y = value), fill = NA, show.legend = F) +
    geom_point(aes(x = Var2, y = orig), colour = 'red') +
    theme_bw() +
    facet_wrap(~Var2, scales = 'free') +
    labs(x = 'Parameter', y = 'Value', title = paste('Sample method =', attr(fitted_bootstrap, 'type')))
  
  return(bp)
}
