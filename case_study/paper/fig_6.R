library(data.table); library(ggplot2); library(lmom)

source('R/auxiliary_functions/imports.R')

outer_ribbon = c(0.05, 0.95)
inner_ribbon = c(0.25, 0.75)

return.period = c(10, 20, 50, 100)

prbs <- sort(c(outer_ribbon, inner_ribbon))

mx <- data.table(readRDS('data/DV.rds'))

mx <- split(mx, mx$cluster)

MX <- lapply(mx, function(x) data.table(dcast(x[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV')))
MX <- lapply(MX, function(x) x[, year := NULL])

dta.fit <- lapply(MX, function(x) sim(x, dist = 'gpa', trim = c(0,0)))

s <- lapply(dta.fit, function(x) sample(x, length = 500,  type = 'nonpar'))

RES.gc <- list()

for (i in seq_along(mx)) {
  
  f <- fit(s[[i]], dta.fit[[i]])
  
  para <- dta.fit[[i]]$REG
  
  dist <- attr(dta.fit[[i]], 'sim.call')$dist
  
  qs <- seq(.01, 1 - 1/max(return.period)*.5, 1/max(return.period))
  qaux <- data.table(rbindlist(lapply(f, function(x) {data.frame(q = do.call(paste0('q', dist), list(qs, xi = x$REG[1], alpha = x$REG[2], k = x$REG[3])))}),
                               idcol = 'sample'),
                     probs = seq_along(qs))
  q <- qaux[, .(val = quantile(q, c(.5, prbs)),
                q = c('sample_median', 'rib_1_min', 'rib_2_min', 'rib_2_max', 'rib_1_max')), 
            by = probs]
  RES.gc[[i]] <- data.table(dcast(q, probs ~ q, value.var = 'val'), gumbel.variate = -log(-log(qs)))
}

res.gc <- rbindlist(RES.gc, idcol = 'id') 
res.gc[, id := paste('Cluster', id)]

axis.lim <- c(-log(-log(range(qs))), range(qaux$q))

rp.lab <- return.period
rp.x <- -log(-log(1 - 1/rp.lab))
rp.y <- axis.lim[3] + (axis.lim[4] - axis.lim[3])*.05

rp.dta <- data.table(rp.x, rp.y, rp.lab)

(rp <- ggplot(res.gc) +
    geom_ribbon(aes(x = gumbel.variate, ymin = rib_1_min, ymax = rib_1_max), fill = 'grey40', alpha = .4) +
    geom_ribbon(aes(x = gumbel.variate, ymin = rib_2_min, ymax = rib_2_max), fill = 'grey20', alpha = .8) +
    geom_line(aes(x = gumbel.variate, y = sample_median), col = 'white', lwd = .5, linetype = 2) + 
    # geom_point(data = rp.dta, aes(x = rp.x, y = rp.y), shape = '|', size = 1.5) +
    # geom_line(data = rp.dta, aes(x = rp.x, y = rp.y)) +
    # geom_text(data = rp.dta, aes(x = rp.x, y = rp.y*2, label = rp.lab), cex = 2.5) +
    # geom_text(data = rp.dta, aes(x = mean(rp.x[rev(rank(rp.lab))[1:2]*.5]), y = rp.y[1]*3.5), label = 'Return period', cex = 2.5) +
    theme_bw() +
    labs(x = 'Return period (years)', 
         y = 'Value', 
         title = '') +
    theme(strip.background = element_blank(),
          strip.text = element_text(hjust = 0)) +
    facet_wrap(~id, ncol = 3) +
    scale_x_continuous(limits = c(min(rp.x), NA),
                       breaks = rp.x,
                       labels = return.period,
                       expand = c(0.1, 0)) +
    scale_y_continuous(limits = c(2, NA)) +
    coord_fixed(.35))

# ggsave('~/Desktop/rp.pdf', rp, width = 8, height = 3, units = 'in')
