library(data.table); library(ggplot2); library(lmomRFA); library(raster); library(rgdal); library(maptools); library(RColorBrewer)

# EGU

source('R/auxiliary_functions/imports.R')

cluster.number <- 1:3

mx <- data.table(readRDS('data/DV.rds'))

# disc <- list()
# 
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
    geom_polygon(data = polygony[which(!polygony$eval),], aes(x = long, y = lat, group = SP_ID, fill = NA), colour = 'red', lwd = 1) +
    geom_point(data = polygony[which(polygony$SP_ID %in% do.call(c, disc)),], aes(x = X1, y = X2, group = SP_ID, pch = factor(1)), colour = 'black', fill = 'white', size = 5) +
    scale_fill_manual(values = c('gold4','orange','red4'), name = 'Cluster') +
    scale_linetype_manual(values = c(1,1), name = 'Passed \nA-D test') +
    scale_shape_manual(values = 21, name = 'Discordant', label = '') +
    guides(linetype = guide_legend(override.aes = list(fill = NA, colour = c('red', 'grey85'), lwd = c(1, .75))),
           fill = guide_legend(override.aes = list(colour = 'grey85'))) +
    coord_fixed(ratio = 1.5) + 
    labs(x = 'Longtitude', y = 'Latitude', title = 'Catchment regionalization', subtitle = 'Dicsordancy and Anderson-Darling test') +
    scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
    theme_bw() +
    theme(legend.position = 'bottom'))

# ggsave('~/Desktop/mapa2.pdf',ad.map, scale = 1.5)

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

cols <- terrain.colors(500)

(dem <- ggplot(data = dem.df) +
    geom_raster(aes(x = X, y = Y, fill = Z), show.legend = F) +
    geom_polygon(data = cr.f, aes(x = long, y = lat, colour = factor(1)), fill = NA, lwd = .75) +
    scale_fill_gradientn(colours = cols) +
    scale_colour_manual(values = 'grey5', labels = 'State \nborder', name = '') +
    coord_fixed(ratio = 1.5) + 
    labs(x = 'Longtitude', y = 'Latitude', title = 'Digital elevation model of the Czech Republic', subtitle = '') +
    scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
    theme_bw() +
    theme(legend.position = 'bottom'))

# mapa <- gridExtra::grid.arrange(ad.map, dem, ncol = 2)
# ggsave('~/Desktop/mapa1.pdf', ad.map)
ggsave('~/Desktop/mapa1.pdf', dem, scale = 1.5)
# ggsave('~/Desktop/mapa.pdf', mapa, scale = 1.5)
