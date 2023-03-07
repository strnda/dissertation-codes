library(data.table); library(ggplot2); library(lmomRFA)

source('R/auxiliary_functions/imports.R')

mx <- data.table(readRDS('data/DV.rds'))
aux <- data.table((readRDS('paper/data/mapa.RDS')))
mr <- data.table((readRDS('paper/data/lmrd.rds')))

setkey(aux)

aux <- unique(aux[,.(SP_ID, cluster, eval)])

MX <- data.table(dcast(mx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
MX <- MX[, year := NULL]

moments <- data.table(t(apply(MX, 2, samlmu)), keep.rownames = T)

names(moments)[2:5] <- paste('m', 1:4, sep = '_')

(lmrd <- ggplot() +
    geom_line(data = mr, aes(x = tau_3, y = value, group = variable, lty = variable), colour = 'black') +
    geom_point(data = merge(moments, aux, by.x = 'rn', by.y = 'SP_ID'),
               aes(x = m_3, y = m_4, fill = eval, shape = factor(cluster)), 
               colour = 'grey5',
               size = 2.5) +
    scale_linetype_manual(values = c('dotted', 'twodash'), name = 'Distribution', labels = c('GEV', 'GPD')) +
    scale_shape_manual(values = c(21, 22, 25), name = 'Cluster') +
    scale_fill_manual(values = c(NA, 'grey45'), name = expression(Passed~A^2~test)) +
    guides(shape = guide_legend(override.aes = list(fill = c('grey55'))),
           fill = guide_legend(override.aes = list(shape = 21,
                                                   fill = c(NA, 'grey55')))) +
    lims(x = c(min(moments[, m_3])*.5, max(moments[, m_3])*1.05), 
         y = c(min(moments[, m_4])*.75, max(moments[, m_4])*1.15)) +
    theme_bw() +
    labs(x = 'L-skewness', y = 'L-kurtosis', title = 'L-moment ratio diagram'))

# ggsave('~/Desktop/lmrd.pdf',lmrd)
