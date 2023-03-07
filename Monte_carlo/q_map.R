library(data.table); library(ggplot2); library(ggpubr)

source(file = "./R/aux_fun.R")

q <- 175

scale <- 150
shape1 <- .25
shape2 <- .75

{
  
  p <- pggamma(q = q,
               scale = scale,
               shape1 = shape1,
               shape2 = shape2)
  
  dta <- data.table(x = 0:round(x = qggamma(p = .99,
                                            scale = scale, 
                                            shape1 = shape1, 
                                            shape2 = shape2)))
  
  (gg <- ggplot(data = dta, 
                mapping = aes(x = x)) +
     stat_function(fun = function(x) pggamma(q = x,
                                             scale = scale,
                                             shape1 = shape1,
                                             shape2 = shape2),
                   n = 1000) +
     geom_segment(x = q, 
                  y = 0, 
                  xend = q, 
                  yend = p,
                  arrow = arrow(length = unit(x = .025,
                                              units = "npc"), 
                                ends = "last", 
                                type = "closed"),
                  size = .1, 
                  # linetype = 2, 
                  colour = "red4") +
     geom_segment(x = q, 
                  y = p, 
                  xend = 0, 
                  yend = p,
                  arrow = arrow(length = unit(x = .025,
                                              units = "npc"), 
                                ends = "last", 
                                type = "closed"),
                  size = .1, 
                  # linetype = 2,
                  colour = "royalblue3") +
     lims(y = 0:1) +
     labs(x = "Value", 
          y = "p",
          title = "Distribution function of the \ngeneralised gamma distribution") +
     theme_bw() + 
     theme(aspect.ratio = 1))
  
  
  
  (n <- ggplot(data = dta, 
               mapping = aes(x = x)) +
      stat_function(fun = function(x) qnorm(p = x),
                    n = 1000) +
      geom_segment(x = p, 
                   y = -10, 
                   xend = p, 
                   yend = qnorm(p = p),
                   arrow = arrow(length = unit(x = .025,
                                               units = "npc"), 
                                 ends = "last", 
                                 type = "closed"),
                   size = .1, 
                   # linetype = 2,
                   colour = "royalblue3") +
      geom_segment(x = p, 
                   y = qnorm(p = p), 
                   xend = 0, 
                   yend = qnorm(p = p),
                   arrow = arrow(length = unit(x = .025,
                                               units = "npc"), 
                                 ends = "last", 
                                 type = "closed"),
                   size = .1, 
                   # linetype = 2, 
                   colour = "red4") +
      lims(x = 0:1) +
      labs(y = "Value", 
           x = "p",
           title = "Quantile function of the \nnormal distribution") +
      theme_bw() + 
      theme(aspect.ratio = 1))
  
  (qm <- ggarrange(gg, n, 
                   nrow = 1, 
                   align = "hv"))
  
}

ggsave(filename = "./figs/qm.png",
       plot = qm,
       width = 10,
       height = 5, 
       units = "in", 
       dpi = 300)
