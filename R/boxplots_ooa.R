library(ggplot2)
library(magrittr)
library(stringr)
# tf2 <- readRDS("C:/Users/petr/Downloads/tbl.f.final(1).rds")
kernelSizeFixed <- 1
tf3 <- tf2 %>% mutate(observerWeight = as.numeric(str_sub(version.group, -5, -1)) / 1000)
pdf("C:/Users/petr/Downloads/igaGrafy/OOA_obrerver-weight-change-effect_kernelSize-1.pdf", onefile = TRUE)
for (sp in sps) {
  df <- tf3 %>%
    filter(species == sp) %>%
    # filter(species == "Accipiter gentilis")  %>%
    ungroup() %>%
    filter(version == "OOA") %>%
    # filter(fc == "L") %>%
    filter(adjust == kernelSizeFixed) %>% # 0.25
    mutate(AA = auc.train * AICc) %>%
    dplyr::select(auc.train, AUC, AICc, adjust, AA, tune.args, rm, fc, observerWeight)

  df$tune.args %<>% as.factor
  df$adjust %<>% as.factor
  df$observerWeight %<>% as.factor
  df$fc %<>% as.factor
  df$rm %<>% as.factor

  gp <- ggplot(df, aes(x = observerWeight, y = AUC, group = observerWeight)) +
    geom_boxplot(aes(color = observerWeight), linewidth = 0.3, outlier.size = 0.3) + # notch = TRUE,
    facet_grid(vars(rm), vars(fc)) +
    theme_bw(
      base_line_size = 0.2,
      base_rect_size = 0.2
    ) +
    labs(
      title = paste0(sp, " | OOA bias correction, observers' weight change effect"),
      subtitle = "feature class",
      caption = paste0("NDOP (Czechia), Apr-Jun, 2018-2021, MaxNet OOA bias correction, fixed kernel size: ", kernelSizeFixed)
    ) +
    # scale_y_continuous(sec.axis = sec_axis(name = "regularisation multiplier")) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "regularisation multiplier", labels = NULL, breaks = NULL)) +
    theme(
      plot.title = element_text(face = "bold.italic"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    xlab("observers' weight (share of occupied squares per grid size km)") +
    ylab("AUC (LSD test)")

  print(gp)
}
dev.off()