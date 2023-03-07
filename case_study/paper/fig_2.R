library(data.table); library(ggplot2); library(gridExtra)

ribbon <- readRDS('paper/data/ribbon_fixed.rds')

val <- unique(ribbon$VAL)
jmena <- paste0('plot_',1:length(val))

for(i in 1:length(val)){
  dta <- ribbon[ribbon$VAL == val[i],] 
  
  if(val[i] == 'D') {
    osa_y <- expression(paste(italic('D '),'[mm]')) 
    nadpis <- 'Event severity'
  }
  
  if(val[i] == 'L') {
    osa_y <- expression(paste(italic('L '),'[month]'))
    nadpis <- 'Event length'
  }
  
  if(val[i] == 'I') {
    osa_y <- expression(paste(italic('I '),'[mm/month]'))
    nadpis <- 'Event intensity'
  }
  
  if(val[i] == 'rD') {
    osa_y <- expression(paste(italic('rD '),'[-]'))
    nadpis <- 'Relative severity'
  }
  
  if(val[i] == 'rI') {
    osa_y <- expression(paste(italic('rI '),group('[',t^-1,']')))
    nadpis <- 'Relative intensity'
  }
  
  gg <- ggplot(dta) +
    geom_line(aes(x = P, y = MEAN, group = factor(VAR), colour = factor(VAR), lty = factor(VAR)), size = 0.5) +
    geom_ribbon(aes(x = P,ymin = Q25, ymax = Q75, fill = factor(VAR), colour = factor(VAR), lty = factor(VAR)), alpha = 0.25) +
    scale_fill_manual(values = c('grey25', 'grey55'), labels = c('Observed values', 'Simulated values'), name = NULL) +
    scale_colour_manual(values = c('grey25', 'grey55'), labels = c('Observed values', 'Simulated values'), name = NULL) +
    scale_linetype_manual(values = c('solid', 'dashed'), labels = c('Observed values', 'Simulated values'), name = NULL) +
    ylab(osa_y) +
    ggtitle(nadpis) +
    xlab(expression(italic('p'))) +
    theme_bw() +
    theme(aspect.ratio = .75,
          plot.margin = unit(c(0, 0, 0, 0), 'cm'))
  
  assign(jmena[i], gg)
}

get_legend <- function(myggplot) {
  
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
  legend <- tmp$grobs[[leg]]
  
  return(legend)
}

legenda <- get_legend(gg)

val <- grid.arrange(plot_1 + theme(legend.position = 'none'),
                    plot_2 + theme(legend.position = 'none'),
                    plot_3 + theme(legend.position = 'none'),
                    plot_4 + theme(legend.position = 'none'),
                    plot_5 + theme(legend.position = 'none'),
                    legenda,
                    ncol = 3)

# ggsave('~/Desktop/val.pdf',val)