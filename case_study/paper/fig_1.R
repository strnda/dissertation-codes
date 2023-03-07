library(data.table); library(ggplot2); library(lmomRFA); library(raster); library(rgdal); library(maptools); library(RColorBrewer)

source('R/auxiliary_functions/imports.R')

cluster.number <- 1:3

mx <- data.table(readRDS('data/DV.rds'))

disc <- list()

# for (i in cluster.number) {
#   
#   mxx <- mx[which(cluster == cluster.number[i]),]
#   
#   MX <- data.table(dcast(mxx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
#   MX <- MX[, year := NULL]
#   
#   disc[[i]] <- names(MX)[which(regtst(regsamlmu(MX))$D >= 3)]
# }
# 
# saveRDS(disc, 'paper/data/disc.RDS')

disc <- readRDS('paper/data/disc.RDS')

polygony <- readRDS('paper/data/mapa.RDS')

# unique(setkey(as.data.table(polygony))[,.(SP_ID, p, eval, cluster)])

(ad.map <- ggplot() +
    geom_polygon(data = polygony, aes(x = long, y = lat, group = SP_ID, fill = factor(cluster), lty = eval), colour = 'grey95', lwd = .75) +
    geom_polygon(data = polygony[which(!polygony$eval),], aes(x = long, y = lat, group = SP_ID, fill = NA, lty = eval), colour = 'black', lwd = .5) +
    geom_point(data = polygony[which(polygony$SP_ID %in% do.call(c, disc)),], aes(x = X1, y = X2, group = SP_ID, pch = factor(1)), colour = 'black', fill = 'white', size = 3) +
    scale_fill_manual(values = c('grey25','grey50','grey75'), name = 'Cluster') +
    scale_linetype_manual(values = c(2,1), name = expression(Passed~A^2~test)) +
    scale_shape_manual(values = 21, name = 'Discordant', label = '') +
    guides(linetype = guide_legend(override.aes = list(fill = NA, colour = c('black', 'grey85'), lwd = c(.5, .75))),
           fill = guide_legend(override.aes = list(colour = 'grey85'))) +
    coord_fixed(ratio = 1.5) + 
    labs(x = 'Longtitude', y = 'Latitude', title = 'Catchment regionalization', subtitle = 'Dicsordance and Anderson-Darling test') +
    scale_x_continuous(limits = c(11.75, 19),
                       breaks = seq(12, 18, 2),
                       labels = paste0(seq(12, 18, 2), '.0'),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(48.45, 51.1),
                       breaks = 48:51,
                       labels = paste0(48:51, '.0'),
                       expand = c(0, 0)) +
    theme_bw() +
    theme(legend.position = 'bottom'))

# dem <- raster('paper/data/DMR3_wgs84_cr/DMR3_wgs84_cr.tif')
# dem.agg <- aggregate(dem, fact = 3)
# 
# dem.df <- data.frame(rasterToPoints(dem.agg))
# colnames(dem.df) <- c('X','Y','Z')
# 
# saveRDS(dem.df, 'paper/data/DMR3_wgs84_cr/dem.RDS')

cr <- readOGR('paper/data/CZE_adm_shp/CZE_adm0.shp', 'CZE_adm0', verbose = F)
cr.f <- fortify(cr)

dem.df <- readRDS('paper/data/DMR3_wgs84_cr/dem.RDS')
summary(dem.df)
cols <- gray.colors(500, start = .15, end = .85)

(dem <- ggplot(data = dem.df) +
    geom_raster(aes(x = X, y = Y, fill = Z)) +
    geom_polygon(data = cr.f, aes(x = long, y = lat, colour = factor(1)), fill = NA, lwd = .75) +
    scale_fill_gradientn(colours = cols, 
                         limits = c(-5,255), 
                         breaks = seq(0, 250, length.out = 6),
                         labels = seq(100, 1600, length.out = 6),
                         guide = guide_colourbar(nbin = 100, 
                                                 barwidth = 15, 
                                                 draw.ulim = FALSE, 
                                                 draw.llim = FALSE, 
                                                 ticks = FALSE, 
                                                 frame.colour = 'black', 
                                                 title.vjust = .8, 
                                                 order = 1),
                         name = 'Elevation (m a.s.l.)') +
    scale_colour_manual(values = 'grey5', labels = 'State \nborder', name = '') +
    coord_fixed(ratio = 1.5) + 
    labs(x = 'Longtitude', 
         y = 'Latitude', 
         title = 'Digital elevation model of the Czech Republic') +
    scale_x_continuous(limits = c(12, 19),
                       breaks = seq(12, 18, 2),
                       labels = paste0(seq(12, 18, 2), '.0'),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(48.5, 51.1),
                       breaks = 48:51,
                       labels = paste0(48:51, '.0'),
                       expand = c(0, 0)) +
    theme_bw() +
    theme(legend.position = 'bottom'))

mapa <- gridExtra::grid.arrange(dem, ad.map, ncol = 2)

# ggsave('~/Desktop/mapa.pdf', mapa, scale = 1.75, width = 8, height = 3, units = 'in')
