start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "geosphere", "ENMeval", "sdm", "ggplot2", "grid", "gridExtra") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

gcfl <- function() {
    this_file <- commandArgs() %>%
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
    if (length(this_file) == 0) {
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
} # https://stackoverflow.com/a/55322344

path.wd <- paste0(gcfl(), "/")

# "C:\Program Files\R\R-4.2.1\bin\x64\Rscript.exe" "C:\Users\petr\Documents\2023-03-20\RP\RP\R\ndopTGOBeval.R"

# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- paste0(path.wd, "../../projects-data/") # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.prep <- paste0(path.wd, "../../dataPrep/")
source(paste0(path.wd, "shared.R"))
path.rgee <- paste0(path.wd, "../../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.lsd <- paste0(path.prep, "lsd/")
path.eval <- paste0(path.prep, "ndopTGOBeval/")
path.img <- paste0(path.eval, "img/")

############
# inputs
############
# uložit si results, rastery a udělat evaluate; ukládat i počty presencí a bg (z modelů) a WKT presencí (NDOP i LSD)
lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds")) %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci

############
# settings
############

ndop.fs <- list("months" = c(4:6), "years" = c(2019:2022), "precision" = 1000, "version" = "v1")

############
# execution
############


#
# pokud je nějaká metoda korelovaná s LSD AUC, tak to znamená, že asi opravdu funguje korekce biasu? ne, jen že se lze spolehnout na výsledky auc.var.avg
#
# mě zajímá především korelace k nejlepším výsledkům, tyto totiž vybírám - dat do vztahu, jestli jsou vyšší AUC i lépe korelované
#

# ag <- as_tibble(readRDS("/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/dataPrep/ndopTGOBeval/ssos.t_Accipiter gentilis_1.rds"))
# # cor(ag$auc.val.avg, ag$AUC)

# ### mělo by stačit to načíst vše a pak zgroupoivat - nepotřebuju cykly
# cors <- list()
# for (v in unique(ag$version)) {
#   print(v)
#   temp.v <- ag %>% filter(version == v)
#   cors[[v]] <- cor(temp.v$auc.val.avg, temp.v$AUC)
# }

# cors.t <- t(as.data.frame(cors))
# cors.names <- rownames(cors.t)
# cors.t <- as_tibble(cors.t)
# cors.t$version <- cors.names

# plot(`ssos.t_Accipiter gentilis_1`$auc.val.avg, `ssos.t_Accipiter gentilis_1`$AUC)


"%notin%" <- Negate("%in%")


modelsResults <- paste0(path.eval, "tbl.rep.rds")
if (file.exists(modelsResults)) {
  print("načítám existující tbl.rep")
  tbl <- readRDS(modelsResults)
} else {
  print("generuju výsledky")
  # pro znovu vygenerování
  rds_list <-
    list.files(
      path.eval,
      pattern = paste0("^t_"),
      ignore.case = TRUE,
      full.names = TRUE
    )

  first <- TRUE
  for (rdsPath in rds_list) {
    bn <- basename(rdsPath)
    r <- unlist(strsplit(unlist(strsplit(bn, "[.]"))[1], "_"))[5]

    print(bn)
    print(r)

    rds <- readRDS(rdsPath)
    rds$rep <- r
    if (first) {
      first <- FALSE
      out <- rds
    } else {
      out %<>% add_row(rds)
    }
  }

  tbl <- as_tibble(out) %>% na.omit()
  saveRDS(tbl, paste0(path.eval, "tbl.rep.rds"))
}



modelsResults.avg <- paste0(path.eval, "tbl.rds")
if (file.exists(modelsResults.avg)) {
  print("načítám existující tbl")
  tbl <- readRDS(modelsResults.avg)
} else {


  # rozparsování verzí do samostatných částí
  tbl %<>%
    ungroup() %>%
    separate_wider_delim(version, "_", names = c("method", "bgSource"), too_few = "align_start", cols_remove = FALSE)

  # tbl.orig  <- tbl

  # průměry z replikací
  tbl.avg <- tbl %>%
    group_by(version, adjust, tune.args, species) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE)
  # vyloučené sloupce doplním znovu zpět
  summarised.not <- setdiff(names(tbl), names(tbl.avg)) # odstranit sloupec rep - nedává smysl
  tbl.not <- tbl %>%
    group_by(version, adjust, tune.args, species) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::select(all_of(summarised.not))

  tbl.avg %<>% add_column(tbl.not)

  tbl <- tbl.avg
  tbl %<>% ungroup() %>% mutate(id = row_number())
  saveRDS(tbl, paste0(path.eval, "tbl.rds"))
}

ss <- "-all"
# selection <- c("cz_tgob", "cz_topS.50", "cz_topS.10", "cz_ssos.0", "ssos.0_ssos.0", "cz_ssos.1.001", "ssos.1.001_ssos.1.001", "cz_ssos.0.topS50", "cz_ssos.0.topS10", "cz_ssos.1.001.topS10", "cz_ssos.2.001.topS10") #
selection <- c("cz_tgob", "cz_topS.50", "cz_topS.10", "cz_ssos.0", "cz_ssos.1.001", "cz_ssos.1.001.topS10") #

# random (null) verze k porovnání
tbl.null.ids <- tbl %>% filter(bgSource == "un" & adjust == "0")

tbl.null <- tbl.null.ids %>%
  group_by(species) %>%
  slice_max(AUC, with_ties = FALSE) %>%
  dplyr::select(AUC, species, id) %>%
  rename(AUC_null = AUC) %>%
  ungroup()


# ostatní verze k ověření zbavené random null
# tbl.all <- tbl %>% filter(version != "rss_un" && adjust != 0)  # funguje divně
tbl.all <- tbl %>% filter(id %notin% tbl.null.ids$id) # v UN verzi zůstávají thinningy

# připojím random a spočtu rozdíl AUC
tbl.all %<>% left_join(tbl.null, by = c("species")) %>%
  mutate(AUCdiff = AUC - AUC_null) %>%
  group_by(species, version) %>%
  slice_max(AUCdiff, with_ties = FALSE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # ############
# tbl.all2 <- tbl.all
# tbl.all2 %<>% filter(version %in% selection) %>%
#   group_by(species, version) %>%
#   slice_max(AUCdiff, with_ties = FALSE) %>%
#   filter(AUC >= 0.7) %>%
#   filter(occs.n < 800) %>%
#   group_by(version) %>%
#   summarise(n = length(id.x))


#
##### # jen druhy co mají smysl... + bych měl možná i dofiltrovat null nad 50%?
#
#  tbl.all %<>% filter(AUC >= 0.7)
# jen druhy do 800 occs
#  tbl.all %<>% filter(occs.n < 800)


tbl.all.median <- tbl.all %>%
  ungroup() %>%
  group_by(version) %>%
  summarise(AUCdiffMedian = median(AUCdiff), AUCdiffQ_025 = quantile(AUCdiff, 0.25), AUCdiffQ_020 = quantile(AUCdiff, 0.20), AUCdiffQ_010 = quantile(AUCdiff, 0.10), AUCdiffQ_005 = quantile(AUCdiff, 0.05)) %>%
  dplyr::select(AUCdiffMedian, AUCdiffQ_025, AUCdiffQ_020, AUCdiffQ_010, AUCdiffQ_005, version) %>%
  ungroup() %>%
  arrange(AUCdiffMedian)

# tbl.all %>% filter( is.na(AUC))
# tbl.all %>% filter( is.na(AUCdiff))

tbl.all %<>% left_join(tbl.all.median, by = "version")

tbl.all.orig <- tbl.all
# vnucení řazení podle mediánu
version_ordered <- with(droplevels(tbl.all), reorder(version, AUCdiffMedian))
tbl.all <- droplevels(tbl.all)
tbl.all$version <- factor(tbl.all$version, levels = levels(version_ordered))


# tbl.all %<>% mutate(version = fct_reorder(version, AUCdiffMedian))


# !!! dočasná výjimka
# tbl.all %<>% filter(version != "kss_buffer") # vychází extrémně špatně, nějaký problém jakým je vytvořena? Notno zkontrolovat...

# boxplotCustom <- function(x) {
#   r <- quantile(x, probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
#   names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
#   return(r)
# }

# a) vše
ggplot(tbl.all %>% ungroup(), aes(x = version, y = AUCdiff)) +
  # stat_summary(fun.data=boxplotCustom, geom="boxplot", lwd=0.1,notch=TRUE)  +
  geom_boxplot(size = 0.1, notch = TRUE, outlier.size = 0.1, outlier.stroke = 0.3) +
  geom_point(aes(y = AUCdiffQ_010), size = 0.2, color = "red", shape = 18) +
  geom_point(aes(y = AUCdiffQ_005), size = 0.2, color = "green", shape = 18) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 4),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(path.img, "tbl.simple", ss, ".png"), width = 1500, height = 1000, units = "px")


# b)
ggplot(tbl.all %>% ungroup() %>% filter(AUCdiffQ_025 >= 0.0), aes(x = version, y = AUCdiff)) +
  # stat_summary(fun.data=boxplotCustom, geom="boxplot", lwd=0.1,notch=TRUE)  +
  geom_boxplot(size = 0.1, notch = TRUE, outlier.size = 0.1, outlier.stroke = 0.3) +
  geom_point(aes(y = AUCdiffQ_010), size = 0.2, color = "red", shape = 18) +
  geom_point(aes(y = AUCdiffQ_005), size = 0.2, color = "green", shape = 18) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 4),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(path.img, "tbl.simple.Q025", ss, ".png"), width = 1500, height = 1000, units = "px")


# c)
ggplot(tbl.all %>% ungroup() %>% filter(AUCdiffQ_020 >= 0.0), aes(x = version, y = AUCdiff)) +
  # stat_summary(fun.data=boxplotCustom, geom="boxplot", lwd=0.1,notch=TRUE)  +
  geom_boxplot(size = 0.1, notch = TRUE, outlier.size = 0.1, outlier.stroke = 0.3) +
  geom_point(aes(y = AUCdiffQ_010), size = 0.2, color = "red", shape = 18) +
  geom_point(aes(y = AUCdiffQ_005), size = 0.2, color = "green", shape = 18) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 4),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(path.img, "tbl.simple.Q020", ss, ".png"), width = 1500, height = 1000, units = "px")


# d)
ggplot(tbl.all %>% ungroup() %>% filter(AUCdiffQ_010 >= 0.0), aes(x = version, y = AUCdiff)) +
  # stat_summary(fun.data=boxplotCustom, geom="boxplot", lwd=0.1,notch=TRUE)  +
  geom_boxplot(size = 0.1, notch = TRUE, outlier.size = 0.1, outlier.stroke = 0.3) +
  geom_point(aes(y = AUCdiffQ_010), size = 0.2, color = "red", shape = 18) +
  geom_point(aes(y = AUCdiffQ_005), size = 0.2, color = "green", shape = 18) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 4),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(path.img, "tbl.simple.Q010", ss, ".png"), width = 1500, height = 1000, units = "px")

# e)
ggplot(tbl.all %>% ungroup() %>% filter(version %in% selection), aes(x = version, y = AUCdiff)) +
  # stat_summary(fun.data=boxplotCustom, geom="boxplot", lwd=0.1,notch=TRUE)  +
  geom_boxplot(size = 0.1, notch = TRUE, outlier.size = 0.1, outlier.stroke = 0.3) +
  geom_point(aes(y = AUCdiffQ_010), size = 0.2, color = "red", shape = 18) +
  geom_point(aes(y = AUCdiffQ_005), size = 0.2, color = "green", shape = 18) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 4),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave(paste0(path.img, "tbl.simple.selection", ss, ".png"), width = 1500, height = 1000, units = "px")



#####
##### pozor, záleží i na tom, jestli je zhoršení z 90 na 80 nebo z 60 na 50!! Nebo zlepšení z 50 na 60 mě taky nemusí zajímat, a z 80 na 90 také ne...
##### Zohlednit! Stejně tak i udělat hranici 0.75 - dofiltrovat, neúspěšné druhy mě také nezajímají!!!
#####


# ggplot(tbl.all %>% ungroup(), aes(x = version, y = AUCdiff)) +
#   geom_violin(draw_quantiles = c(0.05, 0.10, 0.25, 0.5), trim = FALSE, size = 0.2) +
#   geom_jitter(height = 0, width = 0.15, color = "blue", size = 0.1, alpha = 0.3, shape = 16) +
#   theme_light() +
#   theme(
#     legend.position = "none", text = element_text(size = 4),
#     panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
#   ) +
#   scale_y_continuous(limits = c(-0.2, 0.3))
# ggsave(paste0(path.img, "tbl.png"), width = 1500, height = 1000, units = "px")

# vybrat jen ty s dolním kvantilem nad nulou

# je opravitelnost specifická podle počtu presencí?

# a) vše
ggplot(tbl.all %>% ungroup(), aes(occs.n, AUCdiff)) +
  geom_point(aes(colour = factor(version)), size = 0.05) +
  geom_smooth(method = loess, size = 0.2) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 6),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1)
  ) +
  facet_wrap(~version)
ggsave(paste0(path.img, "tbl.trend-prevalence", ss, ".png"), width = 2000, height = 1500, units = "px")


# b) 0.25
ggplot(tbl.all %>% ungroup() %>% filter(AUCdiffQ_025 >= 0.0), aes(occs.n, AUCdiff)) +
  geom_point(aes(colour = factor(version)), size = 0.05) +
  geom_smooth(method = loess, size = 0.2) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 6),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1)
  ) +
  facet_wrap(~version)
ggsave(paste0(path.img, "tbl.trend-prevalence.Q025", ss, ".png"), width = 2000, height = 1500, units = "px")

# c) 0.20
ggplot(tbl.all %>% ungroup() %>% filter(AUCdiffQ_020 >= 0.0), aes(occs.n, AUCdiff)) +
  geom_point(aes(colour = factor(version)), size = 0.05) +
  geom_smooth(method = loess, size = 0.2) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 6),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1)
  ) +
  facet_wrap(~version)
ggsave(paste0(path.img, "tbl.trend-prevalence.Q020", ss, ".png"), width = 2000, height = 1500, units = "px")


# d) 0.10

ggplot(tbl.all %>% ungroup() %>% filter(AUCdiffQ_010 >= 0.0), aes(occs.n, AUCdiff)) +
  geom_point(aes(colour = factor(version)), size = 0.05) +
  geom_smooth(method = loess, size = 0.2) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 6),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1)
  ) +
  facet_wrap(~version)
ggsave(paste0(path.img, "tbl.trend-prevalence.Q010", ss, ".png"), width = 2000, height = 1500, units = "px")


# e) výběr
ggplot(tbl.all %>% ungroup() %>% filter(version %in% selection), aes(occs.n, AUCdiff)) +
  geom_point(aes(colour = factor(version)), size = 0.05) +
  geom_smooth(method = loess, size = 0.2) +
  theme_light() +
  theme(
    legend.position = "none", text = element_text(size = 6),
    panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1)
  ) +
  facet_wrap(~version)
ggsave(paste0(path.img, "tbl.trend-prevalence.selection", ss, ".png"), width = 3000, height = 1500, units = "px")






stop()

tbl3.null <- tbl3 %>%
  filter(group == "un_5000_NA_NA") %>%
  ungroup()
tbl3.best <- tbl3 %>%
  filter(p1 != "un") %>%
  ungroup()


rds_list <-
  list.files(
    path.eval,
    pattern = paste0("^ssos.rds.r.c_"),
    ignore.case = TRUE,
    full.names = TRUE
  )


pdf(paste0(path.eval, "preds/predikce.pdf"), onefile = TRUE, width = 14, height = 6)
for (rdsPath in rds_list) {
  print(rdsPath)
  rds <- readRDS(rdsPath)
  sp <- names(rds)

  tbl3.null.sp <- tbl3.null %>%
    filter(species == sp) %>%
    slice_max(AUC, with_ties = FALSE)

  r.null <- rds[[sp]][[tbl3.null.sp$version]][[tbl3.null.sp$adjust]][["dupld"]][[tbl3.null.sp$tune.args]]


  tbl3.best.sp <- tbl3.best %>%
    filter(species == sp) %>%
    slice_max(AUC, with_ties = FALSE)

  r.best <- rds[[sp]][[tbl3.best.sp$version]][[tbl3.best.sp$adjust]][["dupld"]][[tbl3.best.sp$tune.args]]


  # png(paste0(path.eval, "preds/", sp, ".png"), width = 1500, height = 900)


  par(mfrow = c(1, 2))
  plot(r.null, main = paste0(sp, ", null, AUC=", tbl3.null.sp$AUC))
  plot(r.best, main = paste0(sp, ", best, AUC=", tbl3.best.sp$AUC))
}
dev.off()



stop()
for (sp in unique(tbl3.best$species)) {
  tbl3.null.sp <- tbl3.null %>%
    filter(species == sp) %>%
    slice_max(AUC, with_ties = FALSE)




  r.tbl3.null.sp <- list()
  for (i in 1:nrow(sp.selected)) {
    sp.selected.row <- sp.selected[i, ]
    r.stack.list[[i]] <- rds.r[[sp]][[sp.selected.row$version]][[sp.selected.row$adjust]][[sp.selected.row$duplORnot]][[sp.selected.row$tune.args]]
  }
  r.stack <- stack(r.stack.list)
  r.median <- calc(r.stack, fun = median)

  r.min <- minValue(r.median)
  r.max <- maxValue(r.median)
  r.median <- ((r.median - r.min) / (r.max - r.min))
  r.median <- setMinMax(r.median)





  tbl3.best.sp <- tbl3.best %>%
    filter(species == sp) %>%
    slice_max(AUC, with_ties = FALSE)


  stop()



  rds.out.ndop <- rds.out %>%
    filter(species == sp & version == 1) %>%
    slice_head(n = 1) # potřebuju presence bez thinningu

  ### # LSD validace
  sp.selected <- rds.out.best %>%
    filter(species == sp) %>%
    dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
  print(sp.selected)

  sp.selected.null <- rds.out.null %>%
    filter(species == sp) %>%
    dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
  print(sp.selected.null)

  ### # NDOP validace
  sp.selected.ndop <- rds.out.best.ndop %>%
    filter(species == sp) %>%
    dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
  print(sp.selected)

  sp.selected.null.ndop <- rds.out.null.ndop %>%
    filter(species == sp) %>%
    dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
  print(sp.selected.null.ndop)
  #
  ### # validace nad independent LSD (AUC)
  #
  r.stack.list <- list()
  for (i in 1:nrow(sp.selected)) {
    sp.selected.row <- sp.selected[i, ]
    r.stack.list[[i]] <- rds.r[[sp]][[sp.selected.row$version]][[sp.selected.row$adjust]][[sp.selected.row$duplORnot]][[sp.selected.row$tune.args]]
  }
  r.stack <- stack(r.stack.list)
  r.median <- calc(r.stack, fun = median)

  r.min <- minValue(r.median)
  r.max <- maxValue(r.median)
  r.median <- ((r.median - r.min) / (r.max - r.min))
  r.median <- setMinMax(r.median)

  r.stack.list.null <- list()
  for (i in 1:nrow(sp.selected.null)) {
    sp.selected.row.null <- sp.selected.null[i, ]
    r.stack.list.null[[i]] <- rds.r[[sp]][[sp.selected.row.null$version]][[sp.selected.row.null$adjust]][[sp.selected.row.null$duplORnot]][[sp.selected.row.null$tune.args]]
  }
  r.stack.null <- stack(r.stack.list.null)
  r.median.null <- calc(r.stack.null, fun = median)

  r.min <- minValue(r.median.null)
  r.max <- maxValue(r.median.null)
  r.median.null <- ((r.median.null - r.min) / (r.max - r.min))
  r.median.null <- setMinMax(r.median.null)

  #
  ### # validace nad interním blokem z NDOP (auc.val.avg)
  #
  r.stack.list.ndop <- list()
  for (i in 1:nrow(sp.selected.ndop)) {
    sp.selected.row.ndop <- sp.selected.ndop[i, ]
    r.stack.list.ndop[[i]] <- rds.r[[sp]][[sp.selected.row.ndop$version]][[sp.selected.row.ndop$adjust]][[sp.selected.row.ndop$duplORnot]][[sp.selected.row.ndop$tune.args]]
  }
  r.stack.ndop <- stack(r.stack.list.ndop)
  r.median.ndop <- calc(r.stack.ndop, fun = median)

  r.min <- minValue(r.median.ndop)
  r.max <- maxValue(r.median.ndop)
  r.median.ndop <- ((r.median.ndop - r.min) / (r.max - r.min))
  r.median.ndop <- setMinMax(r.median.ndop)

  r.stack.list.null.ndop <- list()
  for (i in 1:nrow(sp.selected.null.ndop)) {
    sp.selected.row.null.ndop <- sp.selected.null.ndop[i, ]
    r.stack.list.null.ndop[[i]] <- rds.r[[sp]][[sp.selected.row.null.ndop$version]][[sp.selected.row.null.ndop$adjust]][[sp.selected.row.null.ndop$duplORnot]][[sp.selected.row.null.ndop$tune.args]]
  }
  r.stack.null.ndop <- stack(r.stack.list.null.ndop)
  r.median.null.ndop <- calc(r.stack.null.ndop, fun = median)

  r.min <- minValue(r.median.null.ndop)
  r.max <- maxValue(r.median.null.ndop)
  r.median.null.ndop <- ((r.median.null.ndop - r.min) / (r.max - r.min))
  r.median.null.ndop <- setMinMax(r.median.null.ndop)




  pu <- st_as_sf(rds.out.ndop, wkt = "occs.wkt", crs = st_crs(4326))
  pp <- st_as_sf(sp.selected, wkt = "p.wkt", crs = st_crs(4326))
  pa <- st_as_sf(sp.selected, wkt = "a.wkt", crs = st_crs(4326))

  #
  ### # validace nad independent LSD (AUC)
  #
  ggpl1 <- ggplot() +
    geom_raster(data = as.data.frame(r.median, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
    labs(fill = "hab. suitab. \n(cloglog)") +
    theme_light() +
    geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
    geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
    geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
    theme(
      legend.title = element_text(hjust = 0.5, size = 9),
      # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      # plot.subtitle = element_text(hjust = 0.5)
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    ggtitle(
      label = paste0("AUC(NDOP)=", round(max(sp.selected$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected$AUC), digits = 2))
    )


  ggpl2 <- ggplot() +
    geom_raster(data = as.data.frame(r.median.null, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
    labs(fill = "hab. suitab. \n(cloglog)") +
    theme_light() +
    geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
    geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
    geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
    theme(
      legend.title = element_text(hjust = 0.5, size = 9),
      # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      # plot.subtitle = element_text(hjust = 0.5)
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    ggtitle(
      label = paste0("AUC(NDOP)=", round(max(sp.selected.null$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.null$AUC), digits = 2))
    )

  #
  ### # validace nad interním blokem z NDOP (auc.val.avg)
  #
  ggpl3 <- ggplot() +
    geom_raster(data = as.data.frame(r.median.ndop, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
    labs(fill = "hab. suitab. \n(cloglog)") +
    theme_light() +
    geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
    geom_sf(data = pp, pch = 20, size = 0.5, color = "green", show.legend = "point") +
    geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F", show.legend = "point") +
    theme(
      legend.title = element_text(hjust = 0.5, size = 9),
      # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      # plot.subtitle = element_text(hjust = 0.5)
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    ggtitle(
      label = paste0("AUC(NDOP)=", round(max(sp.selected.ndop$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.ndop$AUC), digits = 2))
    )

  ggpl4 <- ggplot() +
    geom_raster(data = as.data.frame(r.median.null.ndop, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
    labs(fill = "hab. suitab. \n(cloglog)") +
    theme_light() +
    geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
    geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
    geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
    theme(
      legend.title = element_text(hjust = 0.5, size = 9),
      # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      # plot.subtitle = element_text(hjust = 0.5)
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    ggtitle(
      label = paste0("AUC(NDOP)=", round(max(sp.selected.null.ndop$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.null.ndop$AUC), digits = 2))
    )


  png(paste0(path.eval, "preds/", sp, "_topXpc-", topXpc, "_topXrows-", topXrows, ".png"), width = 1500, height = 900)
  print(
    grid.arrange(
      ggpl4,
      ggpl3,
      ggpl2,
      ggpl1,
      nrow = 2, ncol = 2,
      top = textGrob(paste0(sp, "\n no bias correction                                                                                             best bias correction"), gp = gpar(fontsize = 20, fontface = "bold")),
      left = textGrob("test indep. (AUC LSD selection)                  test dep. (AUC NDOP selection)", gp = gpar(fontsize = 20, fontface = "bold"), rot = 90),
      bottom = textGrob(
        paste0("NDOP: černá | LSD: zelená=presence, červená=absence) | topXpc=", topXpc, "; topXrows=", topXrows, "\n", paste(unlist(vn[which(names(vn) %in% unique(sp.selected$version))]), collapse = " | "), "\n occs source: AOPK NDOP (years: ", min(ndop.fs$years), "-", max(ndop.fs$years), "; months: ", min(ndop.fs$months), "-", max(ndop.fs$months), ") | author: Petr Balej | WGS84 | grid: KFME (SitMap_2Rad) | generated: ", Sys.Date())
      )
    )
  )
  dev.off()
}

























stop()



unique(tbl2.null$species)

unique(tbl2$species)


tbl2 <- tbl %>%
  rowwise() %>%
  mutate("version2" = ifelse(!is.na(as.numeric(strsplit(version, "_")[[1]][3])) == TRUE, paste0(strsplit(version, "_")[[1]][1], "_", strsplit(version, "_")[[1]][2], "_thin"), version))


tbl2 <- tbl %>%
  rowwise() %>%
  mutate("version3p" = strsplit(version, "_")[[1]][3]) %>%
  mutate("version2" = ifelse(!is.na(as.numeric(version3p)) == TRUE, version3p, version))


tbl2 <- tbl %>%
  rowwise() %>%
  mutate("version_p1" = strsplit(version, "_")[[1]][1]) %>%
  mutate("version_p2" = strsplit(version, "_")[[1]][2]) %>%
  mutate("version_p3" = ifelse(!is.na(as.numeric(strsplit(version, "_")[[1]][3])) == TRUE, "thin", strsplit(version, "_")[[1]][3])) %>%
  mutate("version_p4" = strsplit(version, "_")[[1]][4])

saveRDS(tbl2, paste0(path.eval, "tbl2.RDS"))

tbl2.null <- tbl2 %>%
  filter(version_p1 == "un" & version_p2 == "5000") %>%
  dplyr::select(AUC) %>%
  rename(AUC_null = AUC)




unique(tbl2$version3p)



unique(tbl2$version3p)

!is.na(as.numeric(strsplit("ssss_adas_1.23", "_")[[1]][3])) == TRUE
strsplit("ssss_adas_1.23", "_")[[1]][3]





tblT <- read_csv(paste0(path.eval.tgobEval, "ssos.out.t.csv"))
tblT %>% na.omit()

tbl <- as_tibble(read.csv(paste0(path.eval.tgobEval, "ssos.out.t.csv"))) # read_csv je pro mě problematické (tím jak mi vznikají při zápisu chyby), dostanu jen 10% řádků...
tblX <- tbl
# tbl <- tblX
unique(tbl$AUC)
unique(tbl$fc)

unique(`ssos.out.t.c_Ciconia ciconia_1`$AUC)


names(tbl)

# ENMeval
tbl$rm %<>% as.integer
tbl$fc %<>% as.factor
tbl$tune.args %<>% as.factor
tbl$species %<>% as.factor
tbl$auc.train %<>% as.numeric
tbl$auc.val.avg %<>% as.numeric
tbl$cbi.train %<>% as.numeric
tbl$cbi.val.avg %<>% as.numeric
tbl$AICc %<>% as.numeric
tbl$w.AIC %<>% as.numeric
tbl$delta.AICc %<>% as.numeric
# sdm
tbl$AUC %<>% as.numeric
tbl$Deviance %<>% as.numeric
tbl$cor %<>% as.numeric
tbl$p.value %<>% as.numeric
tbl$sensitivity %<>% as.numeric
tbl$specificity %<>% as.numeric
tbl$TSS %<>% as.numeric
tbl$Kappa %<>% as.numeric
tbl$ccr %<>% as.numeric


tblSF <- tbl %<>% mutate("occs.wkt" = st_as_sf(rds.out.ndop, wkt = "occs.wkt", crs = st_crs(4326)))



tbl %>% filter(!is.na(AUC))
tbl %>% filter(!is.na(fc))


tbl %<>% na.omit()

pu <- st_as_sf(rds.out.ndop, wkt = "occs.wkt", crs = st_crs(4326))


saveRDS(tbl, paste0(path.eval.tgobEval, "ssos.out.t.RDS"))


tbl %<>% filter(!is.na(AUC))

unique(tbl$species)


vn <- versionNames()
vnDf <- versionNamesDf()
write.csv(vnDf, paste0(path.eval, "vnDf.csv"), row.names = FALSE)
# #
# # spojení částí
# #

# # načtení RDS Z ENMeval
# rds_list <-
#     list.files(
#         path.eval.tgobEval,
#         pattern = "^out.*\\.rds$",
#         ignore.case = TRUE,
#         full.names = TRUE
#     )

# first <- TRUE
# # merge částí do jednoho rds
# for (fpath in rds_list) {
#     rds.l <- list()
#     # dočasně rozdělené do dvou částí, nutno spojit:
#     rds.out.temp <- readRDS(fpath)
#     rds.r.temp <- readRDS(str_replace(fpath, "out.t", "rds.r"))

#     if (first) {
#         first <- FALSE
#         rds.out <- rds.out.temp
#         rds.r <- rds.r.temp
#     } else {
#         rds.out %<>% add_row(rds.out.temp)
#         rds.r <- append(rds.r, rds.r.temp)
#     }
# }

# rds.out <- as_tibble(rds.out)
# saveRDS(rds.out, paste0(path.eval, "all_rds.out.rds"))
# saveRDS(rds.r, paste0(path.eval, "all_rds.r.rds"))



# ## napárování 20-22 rasterů
# rr.orig <- readRDS(paste0(path.eval, "all_rds.r-copy.rds"))
# rr.v2022 <- readRDS(paste0(path.eval, "all_rds.rc.rds"))
# rr.orig.names <- names(rr.orig)
# for (sp in rr.orig.names) {
#     print(sp)
#     for (version in names(rr.v2022[[sp]])) {
#         print(version)
#         rr.orig[[sp]][[version]] <- rr.v2022[[sp]][[version]]
#     }
# }
# # saveRDS(rr.orig, paste0(path.eval, "all_rds.r.rds"))




rds.out <- readRDS(paste0(path.eval, "all_rds.out.rds"))
rds.out %<>% filter(!is.na(AUC)) # odstranit modely kde nešla provést LSD evaluace


# porovnat 2 samostatné random null:
rds.out.null.compare <- rds.out %>%
  filter((version == 1 | version == 3) & adjust == 0) %>%
  group_by(species, version, duplORnot) %>%
  mutate(compareGroup = paste0(species, "-", version, "-", duplORnot)) %>%
  ungroup() %>%
  group_by(compareGroup) %>%
  slice_max(AUC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(species, duplORnot) %>%
  mutate(Diff = AUC - lag(AUC)) %>%
  na.omit()

visDiff1 <- rds.out.null.compare %>%
  filter(duplORnot == "dupl") %>%
  dplyr::select(Diff)
visDiff2 <- rds.out.null.compare %>%
  filter(duplORnot == "uniq") %>%
  dplyr::select(Diff)

# boxplot(list("dupl" = visDiff1$Diff, "uniq" = visDiff2$Diff))
# hist(abs(visDiff1$Diff))
# hist(abs(visDiff2$Diff))

#
# dif vůči random
#

# nulová bg random verze
vs.null <- rds.out %>%
  filter(version == 3 & adjust == 0) %>% # vybírám 3, ale mohl bych i 1, reálně musím udělat více replikací (změn backgroundu) a zprůměrovat je
  group_by(species, duplORnot) %>%
  slice_max(AUC, n = 1, with_ties = FALSE) # zde neřeším více modelů se stejně vysokým AUC


# porovnání nej z jednotlivých základních verzí - neřeším adjusty ani šířku thinu!!!  (4176 / 2 / 116 = 18 verzí - mám ale 19 verzí - všechny nemají dupl/unique varianty - nemůžu dělit 2...)

vs.all <- rds.out %>%
  filter(adjust != 0) %>% # ponechávám všechny interní random verze ("0") - neměl bych dělat porovnání s interními random verzemi?
  group_by(species, version, duplORnot) %>%
  slice_max(AUC, n = 1, with_ties = FALSE) # %>% # zde neřeším více modelů se stejně vysokým AUC
# mutate(compareGroup = paste0(species, "-", version, "-", duplORnot))

vs.all %<>% left_join(vs.null %>% dplyr::select(AUC), by = c("species", "duplORnot"), suffix = c("", ".null"))

vs.all.versions <- unique(vs.all$version)

first <- TRUE
for (v in vs.all.versions) {
  vs.all.temp <- vs.all %>%
    ungroup() %>%
    filter(version == v) %>%
    mutate(AUC.diff = AUC - AUC.null)

  if (first) {
    first <- FALSE
    vs.diff <- vs.all.temp
  } else {
    vs.diff %<>% add_row(vs.all.temp)
  }
}

# připojím druhotné skupiny k verzím
vs.diff %<>% left_join(vnDf, by = c("version" = "id")) %>% arrange(group)

# pivot_wider?
vs.all.compare <- vs.diff %>%
  ungroup() %>%
  group_by(species, version, duplORnot)
# %>% summarise(AUC.diff.median = median(AUC.diff))

vs.all.compare.order <- vs.all.compare %>%
  ungroup() %>%
  group_by(version, duplORnot) %>%
  summarise(AUC.diff.median = median(AUC.diff), group = first(group), short = first(short)) %>%
  dplyr::select(AUC.diff.median, version, duplORnot, group, short) %>%
  arrange(AUC.diff.median, duplORnot)

write.csv(vs.all.compare.order, paste0(path.eval, "vs.all.compare.order.csv"), row.names = FALSE)


# # # # # # udělat místo toho per species rozdíly dupl/uniq a ty dát do boxplotů?

dOu <- "dupl"
# uniq vs. dupl - jsou rozdíly + dupld
ggplot(vs.all.compare, aes(x = version, y = AUC.diff, fill = duplORnot)) +
  geom_boxplot()
ggsave(paste0(path.eval, "duplORnot.png"))


### facety?
# ggplot(vs.all.compare, aes(x = version, y = AUC.diff, fill = duplORnot)) + geom_boxplot() + facet_wrap(~group)
# ggplot(vs.all.compare, aes(x = short, y = AUC.diff, fill = group)) +
#     geom_boxplot() +
#     facet_wrap(~duplORnot) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# vnucení řazení podle mediánu
group_ordered <- with(droplevels(vs.all.compare.order %>% filter(duplORnot == dOu) %>% na.omit()), reorder(short, AUC.diff.median))
vs.all.compare.reorder <- droplevels(vs.all.compare %>% filter(duplORnot == dOu) %>% na.omit())
vs.all.compare.reorder$short <- factor(vs.all.compare.reorder$short, levels = levels(group_ordered))

# boxploty s verzemi
ggplot(vs.all.compare.reorder, aes(x = short, y = AUC.diff, fill = group)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("#C7FFC7", "limegreen", "darkgreen", "#FAA0A0", "indianred", "darkred")) +
  ggtitle(dOu)
ggsave(paste0(path.eval, "vs.all.compare.reorder_", dOu, ".png"))

# odstraním outliery
ggplot(vs.all.compare.reorder, aes(x = short, y = AUC.diff, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("#C7FFC7", "limegreen", "darkgreen", "#FAA0A0", "indianred", "darkred")) +
  # scale_y_continuous(limits = c(-0.1, 0.2))
  # coord_cartesian(ylim = quantile(vs.all.compare.reorder$AUC.diff, c(0.1, 0.9)))
  coord_cartesian(ylim = c(-0.1, 0.17)) +
  ggtitle(dOu)
ggsave(paste0(path.eval, "vs.all.compare.reorder.outliers_", dOu, ".png"))

# porovnání napříč druhy??? předchozí hrafy jsou hrubé a nezohledňují kterým druhům korekce pomohla, mohlas e u jiných druhů projevit jinak...

stop()
### # načte 16GB komprimovaných rasterů s predikcemi (zabírá téměř 30GB RAM!!!)
# rds.r <- readRDS(paste0(path.eval, "all_rds.r.rds"))

topXpc <- 0.01
topXrows <- 1 # vybírám už jen nejlepší predikci pro každou verzi a duplikaci, neberu všechny...
# výběr nejlepších X procent

#
### # validace nad independent LSD (AUC)
#
# nejlepší z nejlepších
rds.out.best <- rds.out %>%
  group_by(species, version, duplORnot) %>%
  slice_max(AUC, n = topXrows, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(AUC >= (max(AUC) - topXpc)) %>%
  ungroup()

# "null" - random 5000 a všechny presence (verze 1 adjust 0 nebo verze 3 adjust 0 )
rds.out.null <- rds.out %>%
  filter(version == 1 & adjust == 0) %>%
  group_by(species, version, duplORnot) %>%
  slice_max(AUC, n = topXrows, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(AUC >= (max(AUC) - topXpc)) %>%
  ungroup()
#
### # validace nad interním blokem z NDOP (auc.val.avg)
#
# nejlepší z nejlepších
rds.out.best.ndop <- rds.out %>%
  group_by(species, version, duplORnot) %>%
  slice_max(auc.val.avg, n = topXrows, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(auc.val.avg >= (max(auc.val.avg) - topXpc)) %>%
  ungroup()

# "null" - random 5000 a všechny presence (verze 1 adjust 0 nebo verze 3 adjust 0 )
rds.out.null.ndop <- rds.out %>%
  filter(version == 1 & adjust == 0) %>%
  group_by(species, version, duplORnot) %>%
  slice_max(auc.val.avg, n = topXrows, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(auc.val.avg >= (max(auc.val.avg) - topXpc)) %>%
  ungroup()





# kolik bylo úspěšných modelů pro zvolený treshold
mCounts <- rds.out.best %>%
  group_by(species) %>%
  summarise(n = length(species), AUC.max = max(AUC.max))


for (sp in unique(rds.out.best$species)) {
  rds.out.ndop <- rds.out %>%
    filter(species == sp & version == 1) %>%
    slice_head(n = 1) # potřebuju presence bez thinningu

  ### # LSD validace
  sp.selected <- rds.out.best %>%
    filter(species == sp) %>%
    dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
  print(sp.selected)

  sp.selected.null <- rds.out.null %>%
    filter(species == sp) %>%
    dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
  print(sp.selected.null)

  ### # NDOP validace
  sp.selected.ndop <- rds.out.best.ndop %>%
    filter(species == sp) %>%
    dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
  print(sp.selected)

  sp.selected.null.ndop <- rds.out.null.ndop %>%
    filter(species == sp) %>%
    dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
  print(sp.selected.null.ndop)
  #
  ### # validace nad independent LSD (AUC)
  #
  r.stack.list <- list()
  for (i in 1:nrow(sp.selected)) {
    sp.selected.row <- sp.selected[i, ]
    r.stack.list[[i]] <- rds.r[[sp]][[sp.selected.row$version]][[sp.selected.row$adjust]][[sp.selected.row$duplORnot]][[sp.selected.row$tune.args]]
  }
  r.stack <- stack(r.stack.list)
  r.median <- calc(r.stack, fun = median)

  r.min <- minValue(r.median)
  r.max <- maxValue(r.median)
  r.median <- ((r.median - r.min) / (r.max - r.min))
  r.median <- setMinMax(r.median)

  r.stack.list.null <- list()
  for (i in 1:nrow(sp.selected.null)) {
    sp.selected.row.null <- sp.selected.null[i, ]
    r.stack.list.null[[i]] <- rds.r[[sp]][[sp.selected.row.null$version]][[sp.selected.row.null$adjust]][[sp.selected.row.null$duplORnot]][[sp.selected.row.null$tune.args]]
  }
  r.stack.null <- stack(r.stack.list.null)
  r.median.null <- calc(r.stack.null, fun = median)

  r.min <- minValue(r.median.null)
  r.max <- maxValue(r.median.null)
  r.median.null <- ((r.median.null - r.min) / (r.max - r.min))
  r.median.null <- setMinMax(r.median.null)

  #
  ### # validace nad interním blokem z NDOP (auc.val.avg)
  #
  r.stack.list.ndop <- list()
  for (i in 1:nrow(sp.selected.ndop)) {
    sp.selected.row.ndop <- sp.selected.ndop[i, ]
    r.stack.list.ndop[[i]] <- rds.r[[sp]][[sp.selected.row.ndop$version]][[sp.selected.row.ndop$adjust]][[sp.selected.row.ndop$duplORnot]][[sp.selected.row.ndop$tune.args]]
  }
  r.stack.ndop <- stack(r.stack.list.ndop)
  r.median.ndop <- calc(r.stack.ndop, fun = median)

  r.min <- minValue(r.median.ndop)
  r.max <- maxValue(r.median.ndop)
  r.median.ndop <- ((r.median.ndop - r.min) / (r.max - r.min))
  r.median.ndop <- setMinMax(r.median.ndop)

  r.stack.list.null.ndop <- list()
  for (i in 1:nrow(sp.selected.null.ndop)) {
    sp.selected.row.null.ndop <- sp.selected.null.ndop[i, ]
    r.stack.list.null.ndop[[i]] <- rds.r[[sp]][[sp.selected.row.null.ndop$version]][[sp.selected.row.null.ndop$adjust]][[sp.selected.row.null.ndop$duplORnot]][[sp.selected.row.null.ndop$tune.args]]
  }
  r.stack.null.ndop <- stack(r.stack.list.null.ndop)
  r.median.null.ndop <- calc(r.stack.null.ndop, fun = median)

  r.min <- minValue(r.median.null.ndop)
  r.max <- maxValue(r.median.null.ndop)
  r.median.null.ndop <- ((r.median.null.ndop - r.min) / (r.max - r.min))
  r.median.null.ndop <- setMinMax(r.median.null.ndop)




  pu <- st_as_sf(rds.out.ndop, wkt = "occs.wkt", crs = st_crs(4326))
  pp <- st_as_sf(sp.selected, wkt = "p.wkt", crs = st_crs(4326))
  pa <- st_as_sf(sp.selected, wkt = "a.wkt", crs = st_crs(4326))

  #
  ### # validace nad independent LSD (AUC)
  #
  ggpl1 <- ggplot() +
    geom_raster(data = as.data.frame(r.median, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
    labs(fill = "hab. suitab. \n(cloglog)") +
    theme_light() +
    geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
    geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
    geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
    theme(
      legend.title = element_text(hjust = 0.5, size = 9),
      # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      # plot.subtitle = element_text(hjust = 0.5)
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    ggtitle(
      label = paste0("AUC(NDOP)=", round(max(sp.selected$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected$AUC), digits = 2))
    )


  ggpl2 <- ggplot() +
    geom_raster(data = as.data.frame(r.median.null, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
    labs(fill = "hab. suitab. \n(cloglog)") +
    theme_light() +
    geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
    geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
    geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
    theme(
      legend.title = element_text(hjust = 0.5, size = 9),
      # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      # plot.subtitle = element_text(hjust = 0.5)
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    ggtitle(
      label = paste0("AUC(NDOP)=", round(max(sp.selected.null$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.null$AUC), digits = 2))
    )

  #
  ### # validace nad interním blokem z NDOP (auc.val.avg)
  #
  ggpl3 <- ggplot() +
    geom_raster(data = as.data.frame(r.median.ndop, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
    labs(fill = "hab. suitab. \n(cloglog)") +
    theme_light() +
    geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
    geom_sf(data = pp, pch = 20, size = 0.5, color = "green", show.legend = "point") +
    geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F", show.legend = "point") +
    theme(
      legend.title = element_text(hjust = 0.5, size = 9),
      # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      # plot.subtitle = element_text(hjust = 0.5)
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    ggtitle(
      label = paste0("AUC(NDOP)=", round(max(sp.selected.ndop$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.ndop$AUC), digits = 2))
    )

  ggpl4 <- ggplot() +
    geom_raster(data = as.data.frame(r.median.null.ndop, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
    labs(fill = "hab. suitab. \n(cloglog)") +
    theme_light() +
    geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
    geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
    geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
    theme(
      legend.title = element_text(hjust = 0.5, size = 9),
      # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      # plot.subtitle = element_text(hjust = 0.5)
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    ggtitle(
      label = paste0("AUC(NDOP)=", round(max(sp.selected.null.ndop$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.null.ndop$AUC), digits = 2))
    )


  png(paste0(path.eval, "preds/", sp, "_topXpc-", topXpc, "_topXrows-", topXrows, ".png"), width = 1500, height = 900)
  print(
    grid.arrange(
      ggpl4,
      ggpl3,
      ggpl2,
      ggpl1,
      nrow = 2, ncol = 2,
      top = textGrob(paste0(sp, "\n no bias correction                                                                                             best bias correction"), gp = gpar(fontsize = 20, fontface = "bold")),
      left = textGrob("test indep. (AUC LSD selection)                  test dep. (AUC NDOP selection)", gp = gpar(fontsize = 20, fontface = "bold"), rot = 90),
      bottom = textGrob(
        paste0("NDOP: černá | LSD: zelená=presence, červená=absence) | topXpc=", topXpc, "; topXrows=", topXrows, "\n", paste(unlist(vn[which(names(vn) %in% unique(sp.selected$version))]), collapse = " | "), "\n occs source: AOPK NDOP (years: ", min(ndop.fs$years), "-", max(ndop.fs$years), "; months: ", min(ndop.fs$months), "-", max(ndop.fs$months), ") | author: Petr Balej | WGS84 | grid: KFME (SitMap_2Rad) | generated: ", Sys.Date())
      )
    )
  )
  dev.off()
}


end_time <- Sys.time()
print(end_time - start_time)