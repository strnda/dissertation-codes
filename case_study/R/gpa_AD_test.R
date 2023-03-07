library(data.table); library(ggplot2); library(lmom)

source('R/auxiliary_functions/imports.R')

## set trim and number of at-site samples
trim <- c(0,0)
smp <- 1000

## import data
# mx <- data.table(readRDS('data/DV.rds'))

mx <- readRDS('~/ownCloud/Active Docs/filip_nim/data/upov_mx.rds')
names(mx)[which(names(mx) == 'UPOV_ID')] <- 'SP_ID'

MX <- as.data.table(dcast(mx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
MX <- MX[, year := NULL]

## at-site l-moments and GPA parameters estimation
lmom.atsite <- as.data.frame(t(apply(MX, 2, function(x) samlmu(x, trim = trim))))
# ratiodiagram(lmom.atsite)
# ratiodiagram(lmom.atsite, 'discord')
at.site.pars <- data.table(SP_ID = names(MX), t(apply(lmom.atsite, 1, function(x) {pelgpa(x)})))

# saveRDS(at.site.pars, 'paper/data/para.RDS')
# lmom.atsite <- as.data.frame(t(apply(MX, 2, function(x) lmanual(x)))) # at site l-momenty (prepsanou fci)
# ratiodiagram(lmom.atsite[,3:4])
# at.site.pars <- data.table(SP_ID = names(MX), t(apply(lmom.atsite, 1, function(x) {gpa.para(x)}))) # atsite gpa parametry (prepsanou fci)

## Anderson-Darling statistic test for observed data
para.mx <- merge(mx[, .(SP_ID, dV)], at.site.pars)
ad.mx <- para.mx[, .(base_ad = AD.test(val = dV, 
                                       location = unique(xi), 
                                       scale = unique(alpha), 
                                       shape = unique(k), 
                                       dist = 'gpa')),
                 by = SP_ID]

## Anderson-Darling statistic test for sampled data
AD.SMP <- lapply(1:smp, function(i) { 
  ad.dV <- para.mx[, .(dV = rgpa(length(dV), unique(xi), unique(alpha), unique(k))), by = SP_ID] # samplovani GPA hodnot
  ad.para <- dcast(ad.dV[, .(val = pelgpa(samlmu(dV)), para = c('xi', 'alpha', 'k')), by = SP_ID], SP_ID ~ para, value.var = 'val') # paramety samplu
  ad.res <- merge(ad.dV, ad.para)
  ad.res[, .(smp_ad = AD.test(val = dV, location = unique(xi), scale = unique(alpha), shape = unique(k), dist = 'gpa')), by = SP_ID] # AD stat. pro samply
})

ad.smp <- rbindlist(AD.SMP)

## at-ste p-value calculation
ad <- merge(ad.mx, ad.smp)
res <- unique(ad[, .(p = length(which(smp_ad >= base_ad))/ .N, # vypocty p-hodnoy a vyhodnoceni testu s hl. vyzn. 5%
                     eval = base_ad < quantile(smp_ad, .95)), by = SP_ID]) #######################

# length(which(res$eval))
# 
# summary(res$p)
# summary(res$eval)

ratiodiagram(lmom.atsite[which(!res$eval),]) + ggtitle('ty co neprosly testem')

# saveRDS(res, 'Active Docs/filip_nim/data/AD.rds')
#############################################################

library(rgdal); library(maptools)

## import data
pov <- readOGR('data/mezipovodi_133_cluster.shp', 'mezipovodi_133_cluster', verbose = F)
pov <- spTransform(pov, CRS('+init=epsg:4326'))

centroid <- data.frame(SP_ID = pov@data$DBCN, coordinates(pov)) 

pov.f <- fortify(pov, region = 'DBCN')
names(pov.f)[which(names(pov.f) == 'id')] <- 'SP_ID'

setkey(mx)

## merge all information - spatial polygons, AD statistics, cluster numbers
polygony <- Reduce(function(...) merge(..., all = TRUE, by = 'SP_ID'), list(pov.f, centroid, res, unique(mx[,.(SP_ID, cluster)])))

polygony <- readRDS('paper/data/mapa.RDS')
# saveRDS(polygony, 'paper/data/mapa.RDS')

## AD visualization
(ad.map <- ggplot() +
  geom_polygon(data = polygony, aes(x = long, y = lat, group = SP_ID, fill = factor(cluster)), colour = 'grey95', lwd = .75) +
  geom_polygon(data = polygony[which(!polygony$eval),], aes(x = long, y = lat, group = SP_ID, fill = NA, lty = eval), colour = 'black', lwd = .5, lty = 2) +
  # geom_text(aes(x = X1, y = X2, group = SP_ID,label = round(p, 2)), size = 2.5) +
  scale_fill_manual(values = c('grey25','grey50','grey75'), name = '') +
  coord_fixed(ratio = 1.5) + 
  theme_classic())

# htmlwidgets::saveWidget(plotly::as_widget(plotly::ggplotly(ad.map)), file = paste0(getwd(),'/Active Docs/filip_nim/data/ad_map.html'))
