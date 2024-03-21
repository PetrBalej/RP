library(ggplot2)
library(magrittr)
# tf2 <- readRDS("C:/Users/petr/Downloads/tbl.f.final(1).rds")

pdf("C:/Users/petr/Downloads/igaGrafy/TGOB_kernel-size-change-effect.pdf", onefile = TRUE)
for (sp in sps) {
  df <- tf2 %>%
    filter(species == sp) %>%
    # filter(species == "Accipiter gentilis")  %>%
    ungroup() %>%
    filter(version == "TGOB") %>%
    # filter(str_detect(tune.args, "kernel2")) %>%
    # filter(fc == "L") %>%
    mutate(AA = auc.train * AICc) %>%
    dplyr::select(auc.train, AUC, AICc, adjust, AA, tune.args, rm, fc)

  df$tune.args %<>% as.factor
  df$adjust %<>% as.factor
  df$fc %<>% as.factor
  df$rm %<>% as.factor

  gp <- ggplot(df, aes(x = adjust, y = AUC, group = adjust)) +
    geom_boxplot(aes(color = adjust), linewidth = 0.3, outlier.size = 0.3) + # notch = TRUE,
    facet_grid(vars(rm), vars(fc)) +
    theme_bw(
      base_line_size = 0.2,
      base_rect_size = 0.2
    ) +
    labs(
      title = paste0(sp, " | TGOB bias correction, kernel size change effect"),
      subtitle = "feature class",
      caption = "NDOP (Czechia), Apr-Jun, 2019-2022, MaxNet with TGOB bias correction"
    ) +
    # scale_y_continuous(sec.axis = sec_axis(name = "regularisation multiplier")) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "regularisation multiplier", labels = NULL, breaks = NULL)) +
    theme(
      plot.title = element_text(face = "bold.italic"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    xlab("kernel size (x * 3.3km)") +
    ylab("AUC (LSD test)")

  print(gp)
}
dev.off()