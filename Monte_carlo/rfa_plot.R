library(fst); library(data.table); library(ggplot2)

fls_smp <- list.files(path = "./data/reg", 
                      pattern = "smp", 
                      full.names = TRUE)

smp_rfa <- lapply(X = fls_smp, 
                  FUN = read_fst, 
                  as.data.table = TRUE)

names(x = smp_rfa) <- gsub(pattern = "./data/reg/smp_|.fst",
                           replacement = "",
                           x = fls_smp)

smp_rfa <- rbindlist(l = smp_rfa, 
                     idcol = "p0")

fls_org <- list.files(path = "./data/reg", 
                      pattern = "orig", 
                      full.names = TRUE)

org_rfa <- lapply(X = fls_org, 
                  FUN = read_fst, 
                  as.data.table = TRUE)

names(x = org_rfa) <- gsub(pattern = "./data/reg/smp_|.fst",
                           replacement = "",
                           x = fls_smp)

org_rfa <- rbindlist(l = org_rfa, 
                     idcol = "p0")
org_rfa[, para := rep(x = c("xi", "alpha", "k", "0.8", "0.9", "0.95", "0.98", "0.99"),
                      times = 9)]
setnames(x = org_rfa,
         old = "V1", 
         new = "value")

org_para <- org_rfa[para %in% c("xi", "alpha", "k"), ]
org_para <- org_para[, .(value = mean(x = value)), 
                     by = para]

smp_para <- smp_rfa[, .(p0, method, xi, alpha, k)]
smp_para <- melt(data = smp_para, 
                 id.vars = c("p0", "method"),
                 variable.name = "para")

org_para[, para := factor(x = para, 
                          levels = c("xi", "alpha", "k"))]
levels(x = smp_para$para) <- c("xi", "alpha", "k")

par_labs <- c("xi" = "Location",
              "alpha" = "Scale",
              "k" = "Shape")

(p <- ggplot() +
    geom_hline(data = org_para,
               mapping = aes(yintercept = value),
               alpha = .75,
               linetype = 3) +
    geom_boxplot(data = smp_para,
                 mapping = aes(x = p0,
                               y = value, 
                               fill = method), 
                 na.rm = TRUE, 
                 alpha = .65) + 
    facet_wrap(facets = ~ para, 
               ncol = 1, 
               scales = "free",
               labeller = labeller(para = par_labs)) +
    scale_fill_hue(l = 20,
                   name = NULL,
                   labels = paste("Method", 1:3)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(y = "Parameter value",
         title = "Regional parameters"))

ggsave(filename = "./figs/mc/reg_para.png",
       plot = p,
       width = 10, 
       height = 9, 
       units = "in", 
       dpi = 300)

uncer <- smp_rfa[, .(p0, `0.8`, `0.9`, `0.95`, `0.98`, `0.99`)]
uncer_m <- melt(data = uncer, 
                id.vars = "p0")
# uncer_m <- uncer_m[, .(q = 1 - quantile(x = value,
#                                         probs = c(.05, .25, .5, .75, .95),
#                                         na.rm = TRUE),
#                        var = paste0("q_", c(.05, .25, .5, .75, .95))), 
#                    by = variable]
# uncer_d <- dcast(data = uncer_m,
#                  formula = variable ~ var, 
#                  value.var = "q")
uncer_m <- uncer_m[, .(q = 1 - quantile(x = value,
                                        probs = c(.05, .25, .5, .75, .95),
                                        na.rm = TRUE),
                       var = paste0("q_", c(.05, .25, .5, .75, .95))), 
                   by = .(variable, p0)]
uncer_d <- dcast(data = uncer_m,
                 formula = variable + p0 ~ var, 
                 value.var = "q")

(u <- ggplot(data = uncer_d) +
    geom_ribbon(mapping = aes(x = as.numeric(x = variable),
                              ymin = q_0.05,
                              ymax = q_0.95), 
                fill = "grey75",
                alpha = .5) +
    geom_ribbon(mapping = aes(x = as.numeric(x = variable),
                              ymin = q_0.25,
                              ymax = q_0.75), 
                fill = "grey25",
                alpha = .5) +
    geom_line(mapping = aes(x = as.numeric(x = variable),
                            y = q_0.5),
              colour = "grey5",
              alpha = .5) +
    scale_x_continuous(labels = c(5, 10, 20, 50, 100)) +
    theme_bw() +
    facet_wrap(facets = ~p0) +
    theme(legend.position = "bottom") +
    labs(x = "Return period [years]",
         y = "Difference in return level IQRs [-]",
         title = "Uncertainty reduciton with varying p0 (probability of missing values)"))

ggsave(filename = "./figs/mc/uncer.png",
       plot = u,
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300)

x <- precip_max[, -1]

names(x) <- paste0("site_", 1:17)

cx <- cor(x)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

png(height=800, width=800, file="./figs/mc/cor.png", type = "cairo")
  corrplot(cx, method = 'color', col = col(n=200), cl.length = 21, order = 'AOE',
           addCoef.col = 'grey35',type = 'upper', tl.srt=45, tl.col = "grey50")
  # recordPlot()
dev.off()
