library(nim); library(ggplot2); library(data.table); library(kohonen); library(lmomRFA); library(parallel)

# ?nim

data('precip_max')
mx <- precip_max

h <- hclust(d = dist(x = t(x = mx[, -1])))
h <- cutree(tree = h, 
            k = 3)
h
k <- kmeans(x = t(x = mx[, -1]), 
            centers = 3)
k <- k$cluster
k
s <- som(X = scale(x = t(x = mx[, -1])),
         grid = somgrid(xdim = 1, 
                        ydim = 3, 
                        topo = "hexagonal"))
names(x = s$unit.classif) <- names(x = mx[, -1])
s <- s$unit.classif
s

nc <- list(h = table(h), 
           k =table(k),
           s = table(s))

extremes(precip_max) <- -1

# stationary model
n <- nim( ~1, 
          data = precip_max)
smp <- sample.nim(nim = n, 
                  length = 3000)

rp_all <- mclapply(
  X = smp, 
  FUN = function(mx) {
    
    mx <- mx[, -1]
    
    
    h <- hclust(d = dist(x = t(x = mx)))
    h <- cutree(tree = h, 
                k = 3)
    h
    k <- kmeans(x = t(x = mx), 
                centers = 3)
    k <- k$cluster
    k
    s <- som(X = scale(x = t(x = mx)),
             grid = somgrid(xdim = 1, 
                            ydim = 3, 
                            topo = "hexagonal"))
    names(x = s$unit.classif) <- names(x = mx[, -1])
    s <- s$unit.classif
    s
    
    nc <- list(h = table(h), 
               k =table(k),
               s = table(s))
    
    ## cr
    
    # cr <- lapply(
    #   X = c("k", "h", "s"), 
    #   FUN = function(m) {
    #     
    #     cl <- get(x = m)
    #     
    #     cl_all <- lapply(
    #       X = 1:3, 
    #       FUN = function(i) {
    #         
    #         mx[which(x = cl == i)]
    #       }
    #     )
    #     
    #     cr <- sapply(
    #       X = cl_all, 
    #       FUN = function(x) {
    #         
    #         aux <- cor(x = x)
    #         
    #         out <- mean(x = aux[aux != 1])
    #         out
    #       }
    #     )
    #     cr
    #   }
    # )
    # names(x = cr) <- c("k", "h", "s")
    
    
    rp <- lapply(
      X = c("k", "h", "s"), 
      FUN = function(m) {
        
        cl <- get(x = m)
        
        cl_all <- lapply(
          X = 1:3, 
          FUN = function(i) {
            
            mx[which(x = cl == i)]
          }
        )
        
        rp <- sapply(
          X = cl_all, 
          FUN = function(x) {
            
            e <- try(expr = {
              r_mom <- regsamlmu(x = x)
              
              r_para <- t(x = apply(X = r_mom[,-1:-2], 
                                    MARGIN = 1, 
                                    FUN = pelgev))
              
              at_site_q <- t(
                x = apply(X = r_para, 
                          MARGIN = 1, 
                          FUN = function(x) {
                            quagev(f = 1 - (1 / c(5, 10, 20, 50, 100)),
                                   para = x)
                          }
                )
              )
              
              r_fit <- regfit(regdata = r_mom, 
                              dist = "gev")
              
              at_site_q_reg <- sitequant(f = 1 - (1 / c(5, 10, 20, 50, 100)), 
                                         rfd = r_fit)
              
              dif <- apply(X = at_site_q, 
                           MARGIN = 2, 
                           FUN = IQR) / 
                apply(X = at_site_q_reg, 
                      MARGIN = 2, 
                      FUN = IQR)
            }, silent = TRUE)
            
            if (inherits(x = e, 
                         what = "try-error")) {
              
              rep(x = NA, times = 5)
            } else {
              
              dif
            }
          }
        )
        
        as.data.table(rp)
      }
    )
    names(x = rp) <- c("K-means", "H-clust", "SOM")
    
    out <- rbindlist(l = rp, 
                     idcol = "method")
    
  }, mc.cores = detectCores() - 2
)

rptest <- rbindlist(l = rp_all)
rptest[, q := rep_len(x = 1 - (1 / c(5, 10, 20, 50, 100)),
                      length.out = .N)]
write.csv(x = rptest,
          file = "~/Desktop/test.csv")
rptest_m <- melt(data = rptest, 
                 id.vars = c("method", "q"), 
                 variable.name = "cluster")
rptest_m[, cluster := gsub(pattern = "V", 
                           replacement = "Cluster ", 
                           x = cluster)]

znouzecnost <- ggplot(data = rptest_m[value < 1,]) +
  geom_boxplot(mapping = aes(x = factor(x = q),
                             y = value, 
                             fill = method)) +
  geom_hline(yintercept = c(.25, .75), 
             colour = "grey75", 
             linetype = 2) +
  geom_hline(yintercept = c(.5), 
             colour = "grey55", 
             linetype = 2) +
  coord_cartesian(ylim=c(0, 1)) +
  scale_x_discrete(labels = c(5, 10, 20, 50, 100)) +
  facet_wrap(facets = ~cluster,
             ncol = 1) +
  labs(x = "Return period [years]",
       y = "Difference in return level IQRs [-]",
       title = "Uncertainty reduction with for different pooling methods") +
  scale_fill_hue(l = 20,
                 name = "Pooling method") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(filename = "./figs/mc/uncer_pool.png",
       plot = znouzecnost,
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300)
