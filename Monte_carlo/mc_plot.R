library(data.table); library(ggplot2); library(fst)
library(stringr)

fls <- list.files(path = "./data/fitting/", 
                  pattern = ".fst", 
                  full.names = TRUE)

dir.create(path = "./figs/mc/", 
           showWarnings = FALSE, 
           recursive = TRUE)

dist <- sapply(X = str_extract_all(fls,
                                   pattern = "[A-Z]"),
               FUN = paste,
               collapse = "")

dist_orig <- str_sub(string = dist, 
                     start = 1,
                     end = 3)

dist_fitted <- str_sub(string = dist, 
                       start = 4,
                       end = 6)

for (i in seq_along(along.with = fls)) {
  
  # i <- 2
  
  # ifelse(
  #   test = dist_orig[i] == "GPD",
  #   yes = {
  #     loc <- 1
  #     scale_p <- 2
  #     shape_p <- 1},
  #   no = {
  #     loc <- 10
  #     scale_p <- 1.5
  #     shape_p <- .1
  #   }
  # )
  
  mc <- read_fst(path = fls[i], 
                 as.data.table = TRUE)
  
  # summary(mc)
  
  mc_m <- melt(data = mc,
               id.vars = c("par", "par_orig", "n"))
  
  par_q <- mc_m[, .(val = quantile(x = value,
                                   probs = c(.05, .25, .5, .75, .95),
                                   na.rm = TRUE)),
                by = .(variable, par, n)]
  
  par_q[, q := rep(x = c(.05, .25, .5, .75, .95),
                   length.out = .N)]
  
  par_q <- dcast(data = par_q,
                 formula = variable + par + n ~ q,
                 value.var = "val")
  
  o_par <- data.table(par = c("loc", "scale", "shape"),
                      val = c(loc, scale_p, shape_p))
  
  par_q <- merge(x = par_q,
                 y = o_par,
                 by = "par")
  
  method_labs <- c("par_ml" = "Maximum likelihood",
                   "par_l" = "L-moments",
                   "par_n" = "Fitting norms",
                   "p_val_ml" = "Maximum likelihood",
                   "p_val_l" = "L-moments",
                   "p_val_n" = "Fitting norms")
  par_labs <- c("loc" = "Location",
                "scale" = "Scale",
                "shape" = "Shape")
  
  (p <- ggplot(data = par_q[grep(pattern = "par", 
                                 x = variable),]) +
      geom_ribbon(mapping = aes(x = n,
                                ymin = `0.05`,
                                ymax = `0.95`,
                                fill = par),
                  alpha = .15) +
      geom_ribbon(mapping = aes(x = n,
                                ymin = `0.25`,
                                ymax = `0.75`,
                                fill = par),
                  alpha = .35) +
      geom_line(mapping = aes(x = n,
                              y = `0.5`,
                              colour = par),
                alpha = .5) +
      geom_line(mapping = aes(x = n,
                              y = val,
                              colour = par),
                linetype = 2) +
      facet_grid(facets = par ~ variable,
                 scales = "free_y",
                 labeller = labeller(par = par_labs,
                                     variable = method_labs)) +
      scale_fill_hue(l = 20) +
      scale_colour_hue(l = 10) +
      scale_x_log10() +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = "Sample size",
           y = "Parameter value", 
           title = paste(dist_orig[i], "fitted to", dist_fitted[i])))
  
  (a <- ggplot(data = mc_m[grep(pattern = "p_val", 
                                  x = variable),]) +
      geom_boxplot(mapping = aes(x = as.factor(x = n),
                                 y = value, 
                                 group = n),
                   fill = "grey25",
                   alpha = .75,
                   na.rm = TRUE) +
      facet_grid(facets = ~variable,
                 scales = "free_y",
                 labeller = labeller(par = par_labs,
                                     variable = method_labs)) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = "Sample size",
           y = "p-value", 
           title = paste("Anderson-Darling statistic p-values for", dist_orig[i], "sample fitted to", dist_fitted[i])))
  
  ggsave(filename = paste0("./figs/mc/", dist_orig[i], 
                           "_to_", dist_fitted[i], ".png"),
         plot = p,
         width = 10, 
         height = 5, 
         units = "in", 
         dpi = 300)
  
  ggsave(filename = paste0("./figs/mc/AD_", dist_orig[i], 
                           "_to_", dist_fitted[i], ".png"),
         plot = a,
         width = 10, 
         height = 3, 
         units = "in", 
         dpi = 300)
}

