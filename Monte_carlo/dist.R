library(data.table); library(ggplot2); library(extraDistr)

dir.create(path = "./figs/dist/",
           showWarnings = FALSE,
           recursive = TRUE)

dist_all <- c("gev", "gpd")

for (dist in dist_all) {
    
    scl <- c(1, 2)
    ifelse(test = dist == "gev", 
           yes = assign(x = "shp", 
                        value = c(-.5, 0, .5)),
           no = assign(x = "shp", 
                       value = c(1, 5, 20)))

    dta <- data.table(x = seq(from = ifelse(test = dist == "gev", 
                                            yes = -5,
                                            no = .01), 
                              to = 10, 
                              by = 1))
    
    (p <- ggplot(data = dta, 
                 mapping = aes(x = x)) + 
            lapply(X = seq_along(along.with = scl),
                   FUN = function (j) {
                       l <- lapply(X = seq_along(along.with = shp),
                                   FUN = function(i) {
                                       line<- geom_function(fun = function(x) {
                                           do.call(what = paste0("d",
                                                                 dist),
                                                   args = list(x = x,
                                                               mu = 0,
                                                               sigma = scl[j],
                                                               xi = shp[i]))
                                       },
                                       mapping = aes(colour = as.factor(x = paste(i, "shp")),
                                                     linetype = as.factor(paste(j, "scl"))),
                                       n = 1000)
                                       
                                       line
                                   }
                       )
                       l
                   }
            ) +
            geom_point(mapping = aes(y = 0,
                                     shape = as.factor("loc")),
                       colour = "NA") +
            scale_shape(name = NULL,
                        labels = "Location = 0") +
            scale_colour_hue(l = 20, 
                             name = NULL, 
                             labels = paste("Shape =",
                                            shp)) +
            scale_linetype(name = NULL,
                           labels = paste("Scale =",
                                          scl)) +
            theme_bw() +
            labs(x = "Value",
                 y = "PDF",
                 title = paste("Probability density function of", 
                               toupper(x = dist), 
                               "\nwith various values of scale and shape parameters")) +
            theme(legend.position = c(.95, .95),
                  legend.justification = c("right", "top"), 
                  legend.background = element_blank()) +
            guides(shape = guide_legend(order = 1),
                   col = guide_legend(order = 2)))
    
    
    ggsave(filename = paste0("./figs/dist/", 
                             dist, 
                             ".png"),
           plot = p,
           width = 10,
           height = 5, 
           units = "in", 
           dpi = 300)
}


