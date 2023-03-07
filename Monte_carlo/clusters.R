library(data.table); library(ggplot2); library(ggpubr); library(cluster)

dir.create(path = "./figs/cl/", 
           showWarnings = FALSE, 
           recursive = TRUE)

dta <- lapply(
  X = list(c(4, 4),
           c(1, 1),
           c(6, 0)),
  FUN = function (x) {

    n <- sample(x = 145:160,
                size = 1)
    sigma <- matrix(data = c(1.25, 0, 0, 1.25),
                    nrow = 2)

    out <- as.data.table(
      mvtnorm::rmvnorm(n = n,
                       mean = x,
                       sigma = sigma)
    )

    names(x = out) <- c("x", "y")

    out
  }
)

names(x = dta) <- seq_along(along.with = dta)
dta <- rbindlist(l = dta,
                 idcol = "lab")

k <- function(x, k) {
  list(cluster = kmeans(x = x,
                        centers = k))
}

g <- clusGap(x = dta[, .(x, y)], 
             FUN = k, 
             K.max = 8, 
             B = 20)

(gap <- ggplot(data = as.data.frame(x = g$Tab),
               mapping = aes(x = seq_along(along.with = gap),
                             y = gap,
                             ymin = gap - SE.sim,
                             ymax = gap + SE.sim)) +
    geom_line(size = 1,
              alpha = .5) +
    geom_point(size = 3,
               alpha = .75) +
    geom_errorbar(colour = "red4", 
                  width = .2) +
    coord_fixed(ratio = 12.5) +
    labs(x = "Number of clusters k",
         y = "Gap statistic (k)") +
    theme_bw())

ggsave(filename = "./figs/cl/gap.png",
       plot = gap,
       width = 10,
       height = 5.5, 
       units = "in", 
       dpi = 300)


dta_kmeans = dta[, cluster := factor(x = kmeans(x = cbind(x, y),
                                                centers = which.max(x = g$Tab[,3]))$cluster)]

(k1 <- ggplot(data = dta_kmeans,
              mapping = aes(x = x,
                            y = y,
                            shape = lab)) +
    geom_point(size = 2.5,
               alpha = .55, 
               show.legend = FALSE) +
    scale_shape_manual(values = c(15, 16, 17),
                       name = "Group") +
    coord_fixed(ratio = 1.25) +
    theme_bw() +
    labs(title = "Raw data"))

(k2 <- ggplot(data = dta_kmeans,
              mapping = aes(x = x,
                            y = y,
                            shape = lab,
                            colour = cluster)) +
    geom_point(size = 2.5,
               alpha = .55) +
    scale_color_hue(l = 15,
                    name = "Cluster") +
    scale_shape_manual(values = c(15, 16, 17),
                       name = "Group") +
    coord_fixed(ratio = 1.25) +
    theme_bw() +
    labs(title = "Clustered data"))

(cl <- ggarrange(k1, k2, 
                 nrow = 1, 
                 common.legend = TRUE, 
                 legend = "right"))

ggsave(filename = "./figs/cl/kmeans.png",
       plot = cl,
       width = 10,
       height = 5, 
       units = "in", 
       dpi = 300)
####################################x

library(ggplot2); library(ggdendro)

h <- function(x, k) {
  
  d <- dist(x = x, 
            method = "euclidian")
  
  list(cluster = hclust(d = d),
       k = k)
}

g <- clusGap(x = dta[, .(x, y)], 
             FUN = h,
             K.max = 8,
             B = 20)



dta_matrix <- as.matrix(x = dta[, .(x, y)])
rownames(x = dta_matrix) <- dta[, lab]

hc <- hclust(d = dist(x = dta_matrix),
             method = "ave")

dendr <- dendro_data(model = hc,
                     type = "rectangle") # convert for ggplot
clust <- cutree(tree = hc,
                k = 3)                    # find 2 clusters
clust_df <- data.frame(label = names(x = clust), 
                       cluster = factor(x = clust))
# dendr[["labels"]] has the labels, merge with clust.df based on label column
dendr[["labels"]] <- merge(x = dendr[["labels"]],
                           y = clust_df,
                           by = "label")
# plot the dendrogram; note use of color=cluster in geom_text(...)
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
            size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0))
