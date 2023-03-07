library(data.table); library(ggplot2); library(lmom)

polygony <- data.table(readRDS('paper/data/mapa.RDS'))
para <- data.table(readRDS('paper/data/para.RDS'))

setkey(polygony)
setkey(para)

(adgg <- ggplot(polygony) +
    geom_boxplot(aes(x = factor(cluster), y = p, group = cluster), fill = 'grey75', width = .8) +
    geom_hline(yintercept = .05, linetype = 'dashed') +
    labs(x = '', y = expression(paste(italic(p), '-value')), title = expression(paste('At-site ', italic(p), '-values')), subtitle = '') +
    scale_x_discrete(labels = paste('Cluster', 1:3)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          aspect.ratio = 2))

para <- para[unique(polygony[, .(SP_ID, cluster)]), on = 'SP_ID']
para.m <- melt(para[,-1], id.vars = 'cluster')
para.m[, variable := gsub('k', 'kappa', variable)]
para.m[, variable := factor(variable, levels = c('xi', 'alpha', 'kappa'))]

(paragg <- ggplot(para.m) +
    geom_boxplot(aes(x = factor(cluster), y = value, group = cluster), fill = 'grey75', outlier.size = 1, width = 2.5) +
    labs(x = '', y = 'Parameter value', title = 'At-site parameter values') +
    scale_x_discrete(labels = paste('Cluster', 1:3)) +
    scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          aspect.ratio = 1.9) +
    facet_wrap(~variable, scales = 'free_y', labeller = label_parsed))

xxx <- gridExtra::grid.arrange(paragg, adgg, ncol = 2, widths = c(2.5, .91))

# ggsave('~/Desktop/para_ad.pdf', xxx)
