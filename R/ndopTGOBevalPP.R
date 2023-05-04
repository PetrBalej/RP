start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "ENMeval", "sdmATM", "ggplot2", "grid", "gridExtra") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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
path.tgob <- paste0(path.prep, "ndopTGOB/")
path.PP <- paste0(path.eval, "PP/")

############
# inputs
############
# uložit si results, rastery a udělat evaluate; ukládat i počty presencí a bg (z modelů) a WKT presencí (NDOP i LSD)
lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds")) %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci

############
# settings
############

ndop.fs <- list("aucTresholds" = c(0.00, 0.70), "version" = "v1")

############
# execution
############

# select null if performs better
withNull <- FALSE

#
# pokud je nějaká metoda korelovaná s LSD AUC, tak to znamená, že asi opravdu funguje korekce biasu? ne, jen že se lze spolehnout na výsledky auc.var.avg
#
# mě zajímá především korelace k nejlepším výsledkům, tyto totiž vybírám - dat do vztahu, jestli jsou vyšší AUC i lépe korelované
#


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
    r <- unlist(strsplit(unlist(strsplit(bn, "[.]"))[1], "_"))[4]

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

  # # rozparsování verzí do samostatných částí
  # tbl %<>%
  #   ungroup() %>%
  #  separate_wider_delim(version, "_", names = c("method", "bgSource"), too_few = "align_start", cols_remove = FALSE)

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
# vu <- unique(tbl$version)

# tbl %<>% add_row(tbl_5)

# saveRDS(tbl, paste0(path.eval, "tbl_5_6.rds"))
# tbl <- readRDS(modelsResults.avg)

selection <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "un")
# selection.f <- c("ssos.topA100", "ssos.topS10", "tgob", "topA100", "topS10", "un", "ssos") # soos záměrně na konci
selection.f2 <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "un")
selection.rename <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "thin")
names(selection.rename) <- selection.f2

####
# tbl %<>% filter(version %in% selection)

# random (null) verze k porovnání
tbl.null.ids <- tbl %>% filter(version == "un" & adjust == "0")

tbl.null.ids.unique <- unique(tbl.null.ids$id)

tbl.null.test <- tbl.null.ids %>%
  group_by(species) %>%
  slice_max(AUC, with_ties = FALSE) %>%
  dplyr::select(AUC, species, id) %>%
  rename(AUC_null = AUC) %>%
  ungroup()


tbl.null.val <- tbl.null.ids %>%
  group_by(species) %>%
  slice_max(auc.val.avg, with_ties = FALSE) %>%
  dplyr::select(auc.val.avg, species, id) %>%
  rename(auc.val.avg_null = auc.val.avg) %>%
  ungroup()

##########################
#
# musím odlišit/vyčlenit thinningy (zbývající UN)´- ze summary_zasobnik
# kde vyfiltrovat 0.70 AUC? nebo auc.val.avg?
#
#######################


#####
##### pozor, záleží i na tom, jestli je zhoršení z 90 na 80 nebo z 60 na 50!! Nebo zlepšení z 50 na 60 mě taky nemusí zajímat, a z 80 na 90 také ne...
##### Zohlednit! Stejně tak i udělat hranici 0.75 - dofiltrovat, neúspěšné druhy mě také nezajímají!!!
#####

if (withNull) {
  # nepočítám negativní diff, nechávam 0 rozdíl
  # připojím null (val+test) a spočtu diff
  tbl %<>% left_join(tbl.null.test, by = c("species"), suffix = c("", "__AUC")) %>%
    mutate(AUCdiff = ifelse(AUC > AUC_null, AUC - AUC_null, 0)) %>%
    group_by(species, version)

  tbl %<>% left_join(tbl.null.val, by = c("species"), suffix = c("", "__auc.avg")) %>%
    mutate(auc.val.avgdiff = ifelse(auc.val.avg > auc.val.avg_null, auc.val.avg - auc.val.avg_null, 0)) %>%
    group_by(species, version)

  path.PP <- gsub("/$", "_withNull/", path.PP)
} else {
  # připojím null (val+test) a spočtu diff
  tbl %<>% left_join(tbl.null.test, by = c("species"), suffix = c("", "__AUC")) %>%
    mutate(AUCdiff = AUC - AUC_null) %>%
    group_by(species, version)

  tbl %<>% left_join(tbl.null.val, by = c("species"), suffix = c("", "__auc.avg")) %>%
    mutate(auc.val.avgdiff = auc.val.avg - auc.val.avg_null) %>%
    group_by(species, version)
}

dir.create(path.PP, showWarnings = FALSE)

#
# porovnání počtu druhů průchozích přes treshold
#

summary_zasobnik <- list()
zasobnik0 <- list()
combs <- list()


# výběr a přejmenování verzí
tbl.f <- tbl %>%
  filter(version %in% selection.f2) %>%
  ungroup() %>%
  mutate(version = str_replace_all(version, selection.rename)) %>%
  group_by(species, version) # final versions selection

selection.f2 <- selection.rename
k6 <- comb_all(selection.f2, length(selection.f2))

tbl.nn <- tbl.f %>% filter(id %notin% tbl.null.ids.unique) # not null
tbl.null <- tbl.f %>% filter(id %in% tbl.null.ids.unique) # null

for (at in ndop.fs$aucTresholds) {
  print(at)

  ##########################################################################################################################################################################

  #
  # null
  #

  temp.null <- tbl.null %>%
    ungroup() %>%
    group_by(species) %>%
    slice_max(auc.val.avg_null, with_ties = FALSE)
  # null naive
  summary_zasobnik[[as.character(at)]][["null.val"]] <- nrow(temp.null %>% filter(auc.val.avg_null >= at))
  # null naive true
  summary_zasobnik[[as.character(at)]][["null.val.test"]] <- nrow(temp.null %>% filter(auc.val.avg_null >= at) %>% filter(AUC_null >= at))
  # null true
  summary_zasobnik[[as.character(at)]][["null.test"]] <- nrow(temp.null %>% ungroup() %>% group_by(species) %>% slice_max(AUC_null, with_ties = FALSE) %>% filter(AUC_null >= at))


  #
  # not null auc.val.avgdiff
  #
  temp.nn <- tbl.nn %>%
    ungroup() %>%
    group_by(species) %>%
    slice_max(auc.val.avgdiff, with_ties = FALSE)
  # val
  summary_zasobnik[[as.character(at)]][["kss.val"]] <- nrow(temp.nn %>% filter(auc.val.avg >= at))

  # val.test
  summary_zasobnik[[as.character(at)]][["kss.val.test"]] <- nrow(temp.nn %>% filter(auc.val.avg >= at) %>% filter(AUC >= at))

  # true
  summary_zasobnik[[as.character(at)]][["kss.test"]] <- nrow(tbl.nn %>% ungroup() %>% group_by(species) %>% slice_max(AUCdiff, with_ties = FALSE) %>% filter(AUC >= at))

  #
  # kss + null
  #
  temp.nnNull <- tbl.f %>%
    ungroup() %>%
    group_by(species) %>%
    slice_max(auc.val.avgdiff, with_ties = FALSE)
  # val
  summary_zasobnik[[as.character(at)]][["kssNull.val"]] <- nrow(temp.nnNull %>% filter(auc.val.avg >= at))

  # val.test
  summary_zasobnik[[as.character(at)]][["kssNull.val.test"]] <- nrow(temp.nnNull %>% filter(auc.val.avg >= at) %>% filter(AUC >= at))

  # true
  summary_zasobnik[[as.character(at)]][["kssNull.test"]] <- nrow(tbl.f %>% ungroup() %>% group_by(species) %>% slice_max(AUCdiff, with_ties = FALSE) %>% filter(AUC >= at))

  #
  # jen test
  #
  # true - nejlepší dosažitelné a ověřené - korekce nebo null
  summary_zasobnik[[as.character(at)]][["all.test"]] <- nrow(tbl.f %>% ungroup() %>% group_by(species) %>% slice_max(AUC, with_ties = FALSE) %>% filter(AUC >= at))


  ##########################################################################################################################################################################


  # jsou auc.test.avg a AUC u zvolených metod korelované? - až nad filtrovaným datasetem: AUC?
  # měl bych ale spočíst i omision rate nad LSD - tam by měl být smysluplný?

  # korelovanost výsledků auc AUC v rámci metod
  # zjistit jakou korelovsanost vykazuje výsledné spojení všech verzí!!!!!!!!!!!!!!!!!!!!!!

  zasobnik <- list()
  for (ver in unique(tbl.nn$version)) {
    tmp <- tbl.nn %>%
      ungroup() %>%
      filter(auc.val.avg >= at) %>%
      filter(version == ver)
    ct <- cor.test(tmp$auc.val.avg, tmp$AUC)
    zasobnik[[ver]][["cor"]] <- ct[["estimate"]][["cor"]]
    zasobnik[[ver]][["p.value"]] <- ct[["p.value"]]
    zasobnik[[ver]][["ci.min"]] <- ct[["conf.int"]][1]
    zasobnik[[ver]][["ci.max"]] <- ct[["conf.int"]][2]
  }
  zasobnik.t <- t(as_tibble(zasobnik))
  rn <- row.names(zasobnik.t)
  cn <- names(zasobnik[[1]])
  zasobnik.t <- as.data.frame(zasobnik.t)
  zasobnik.t <- as_tibble(sapply(zasobnik.t, as.numeric))
  names(zasobnik.t) <- cn
  zasobnik.t$version <- rn
  zasobnik0[[as.character(at)]][["all"]] <- zasobnik.t %<>% arrange(desc(cor))
  write.csv(zasobnik.t, paste0(path.PP, "korelace-", as.character(at), ".csv"), row.names = FALSE)

  # max verze
  zasobnik <- list()
  for (ver in unique(tbl.nn$version)) {
    tmp <- tbl.nn %>%
      ungroup() %>%
      group_by(species, version) %>%
      slice_max(auc.val.avg, with_ties = FALSE) %>%
      filter(auc.val.avg >= at) %>%
      filter(version == ver)
    ct <- cor.test(tmp$auc.val.avg, tmp$AUC)
    zasobnik[[ver]][["cor"]] <- ct[["estimate"]][["cor"]]
    zasobnik[[ver]][["p.value"]] <- ct[["p.value"]]
    zasobnik[[ver]][["ci.min"]] <- ct[["conf.int"]][1]
    zasobnik[[ver]][["ci.max"]] <- ct[["conf.int"]][2]
  }
  zasobnik.t <- t(as_tibble(zasobnik))
  rn <- row.names(zasobnik.t)
  cn <- names(zasobnik[[1]])
  zasobnik.t <- as.data.frame(zasobnik.t)
  zasobnik.t <- as_tibble(sapply(zasobnik.t, as.numeric))
  names(zasobnik.t) <- cn
  zasobnik.t$version <- rn
  zasobnik0[[as.character(at)]][["max"]] <- zasobnik.t %<>% arrange(desc(cor))
  write.csv(zasobnik.t, paste0(path.PP, "korelace-max-", as.character(at), ".csv"), row.names = FALSE)



  ##########################################################################################################################################################################

  ######
  # grafy
  #####


  ######################
  # kss auc.val.avgdiff
  ######################


  temp.g <- tbl.nn %>%
    ungroup() %>%
    group_by(species, version) %>%
    slice_max(auc.val.avgdiff, with_ties = FALSE) %>%
    filter(auc.val.avg >= at)

  temp.g.median <- temp.g %>%
    ungroup() %>%
    group_by(version) %>%
    summarise(AUCdiffMedian = median(auc.val.avgdiff), AUCdiffQ_025 = quantile(auc.val.avgdiff, 0.25), AUCdiffQ_020 = quantile(auc.val.avgdiff, 0.20), AUCdiffQ_010 = quantile(auc.val.avgdiff, 0.10), AUCdiffQ_005 = quantile(auc.val.avgdiff, 0.05)) %>%
    dplyr::select(AUCdiffMedian, AUCdiffQ_025, AUCdiffQ_020, AUCdiffQ_010, AUCdiffQ_005, version) %>%
    ungroup() %>%
    arrange(AUCdiffMedian)

  temp.g %<>% left_join(temp.g.median, by = "version")

  temp.g.orig <- temp.g
  # vnucení řazení podle mediánu
  version_ordered <- with(droplevels(temp.g), reorder(version, AUCdiffMedian))
  temp.g <- droplevels(temp.g)
  temp.g$version <- factor(temp.g$version, levels = levels(version_ordered))
  # temp.g %<>% mutate(version = fct_reorder(version, AUCdiffMedian))


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(x = version, y = auc.val.avgdiff)) +
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
  ggsave(paste0(path.PP, "boxplot.val.", as.character(at), ".png"), width = 1500, height = 1000, units = "px")


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(occs.n, auc.val.avgdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.05) +
    geom_smooth(method = loess, size = 0.2) +
    theme_light() +
    theme(
      legend.position = "none", text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    ) +
    facet_wrap(~version)
  ggsave(paste0(path.PP, "trend.val.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")



  # porovnání přispění verzí per species
  title <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avg, with_ties = FALSE) %>% ungroup() %>% mutate(title = paste0(species, "\n", round(auc.val.avg_null, digits = 2), " -> ", round(auc.val.avg, digits = 2))) %>% dplyr::select(title)))
  names(title) <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avg, with_ties = FALSE) %>% ungroup() %>% dplyr::select(species)))

  ggplot(temp.g %>% ungroup(), aes(version, auc.val.avgdiff)) +
    geom_bar(stat = "identity", aes(fill = factor(version))) +
    theme_light() +
    theme(
      legend.position = "none", text = element_text(size = 4),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1)),
      panel.margin = unit(0.2, "lines")
    ) +
    facet_wrap(~species, labeller = as_labeller(title))
  ggsave(paste0(path.PP, "version-species.val.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  # přispění kombinací
  temp.g %<>% ungroup()

  counter <- 0
  first <- TRUE
  for (k in k6) {
    counter <- counter + 1
    print(counter)
    print(paste(k, collapse = "|"))
    temp <- temp.g %>%
      filter(version %in% k) %>%
      group_by(species) %>%
      slice_max(auc.val.avgdiff, with_ties = FALSE) %>%
      ungroup() %>%
      summarise(n = length(species), AUCdiffSum = sum(auc.val.avgdiff), mean = mean(auc.val.avgdiff), AUCdiffMedian = median(auc.val.avgdiff), AUCdiffQ_025 = quantile(auc.val.avgdiff, 0.25), AUCdiffQ_020 = quantile(auc.val.avgdiff, 0.20), AUCdiffQ_010 = quantile(auc.val.avgdiff, 0.10), AUCdiffQ_005 = quantile(auc.val.avgdiff, 0.05))
    temp$versionComb <- paste(k, collapse = "|")
    temp$versionCount <- length(k)
    if (first) {
      first <- FALSE
      res <- temp
    } else {
      res %<>% add_row(temp)
    }
  }
  combs[[as.character(at)]][["val"]] <- res %<>% arrange(desc(AUCdiffSum))
  write.csv(res, paste0(path.PP, "combs-val-", as.character(at), ".csv"), row.names = FALSE)


  ######################
  # kss auc.val.avgdiff test
  ######################

  temp.g <- tbl.nn %>%
    ungroup() %>%
    group_by(species, version) %>%
    slice_max(auc.val.avgdiff, with_ties = FALSE) %>%
    filter(auc.val.avg >= at) %>%
    filter(AUC >= at)

  temp.g.median <- temp.g %>%
    ungroup() %>%
    group_by(version) %>%
    summarise(AUCdiffMedian = median(auc.val.avgdiff), AUCdiffQ_025 = quantile(auc.val.avgdiff, 0.25), AUCdiffQ_020 = quantile(auc.val.avgdiff, 0.20), AUCdiffQ_010 = quantile(auc.val.avgdiff, 0.10), AUCdiffQ_005 = quantile(auc.val.avgdiff, 0.05)) %>%
    dplyr::select(AUCdiffMedian, AUCdiffQ_025, AUCdiffQ_020, AUCdiffQ_010, AUCdiffQ_005, version) %>%
    ungroup() %>%
    arrange(AUCdiffMedian)

  temp.g %<>% left_join(temp.g.median, by = "version")

  temp.g.orig <- temp.g
  # vnucení řazení podle mediánu
  version_ordered <- with(droplevels(temp.g), reorder(version, AUCdiffMedian))
  temp.g <- droplevels(temp.g)
  temp.g$version <- factor(temp.g$version, levels = levels(version_ordered))
  # temp.g %<>% mutate(version = fct_reorder(version, AUCdiffMedian))


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(x = version, y = auc.val.avgdiff)) +
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
  ggsave(paste0(path.PP, "boxplot.val.test.", as.character(at), ".png"), width = 1500, height = 1000, units = "px")


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(occs.n, auc.val.avgdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.05) +
    geom_smooth(method = loess, size = 0.2) +
    theme_light() +
    theme(
      legend.position = "none", text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    ) +
    facet_wrap(~version)
  ggsave(paste0(path.PP, "trend.val.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  # porovnání přispění verzí per species
  title <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avg, with_ties = FALSE) %>% ungroup() %>% mutate(title = paste0(species, "\n", round(auc.val.avg_null, digits = 2), " -> ", round(auc.val.avg, digits = 2))) %>% dplyr::select(title)))
  names(title) <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avg, with_ties = FALSE) %>% ungroup() %>% dplyr::select(species)))

  ggplot(temp.g %>% ungroup(), aes(version, auc.val.avgdiff)) +
    geom_bar(stat = "identity", aes(fill = factor(version))) +
    theme_light() +
    theme(
      legend.position = "none", text = element_text(size = 4),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1)),
      panel.margin = unit(0.2, "lines")
    ) +
    facet_wrap(~species, labeller = as_labeller(title))
  ggsave(paste0(path.PP, "version-species.val.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  # přispění kombinací
  temp.g %<>% ungroup()

  counter <- 0
  first <- TRUE
  for (k in k6) {
    counter <- counter + 1
    print(counter)
    print(paste(k, collapse = "|"))
    temp <- temp.g %>%
      filter(version %in% k) %>%
      group_by(species) %>%
      slice_max(auc.val.avgdiff, with_ties = FALSE) %>%
      ungroup() %>%
      summarise(n = length(species), AUCdiffSum = sum(auc.val.avgdiff), mean = mean(auc.val.avgdiff), AUCdiffMedian = median(auc.val.avgdiff), AUCdiffQ_025 = quantile(auc.val.avgdiff, 0.25), AUCdiffQ_020 = quantile(auc.val.avgdiff, 0.20), AUCdiffQ_010 = quantile(auc.val.avgdiff, 0.10), AUCdiffQ_005 = quantile(auc.val.avgdiff, 0.05))
    temp$versionComb <- paste(k, collapse = "|")
    temp$versionCount <- length(k)
    if (first) {
      first <- FALSE
      res <- temp
    } else {
      res %<>% add_row(temp)
    }
  }
  combs[[as.character(at)]][["val.test"]] <- res %<>% arrange(desc(AUCdiffSum))
  write.csv(res, paste0(path.PP, "combs-val.test-", as.character(at), ".csv"), row.names = FALSE)




  ######################
  # kss AUCdiff
  ######################

  temp.g <- tbl.nn %>%
    ungroup() %>%
    group_by(species, version) %>%
    slice_max(AUCdiff, with_ties = FALSE) %>%
    filter(AUC >= at)

  temp.g.median <- temp.g %>%
    ungroup() %>%
    group_by(version) %>%
    summarise(AUCdiffMedian = median(AUCdiff), AUCdiffQ_025 = quantile(AUCdiff, 0.25), AUCdiffQ_020 = quantile(AUCdiff, 0.20), AUCdiffQ_010 = quantile(AUCdiff, 0.10), AUCdiffQ_005 = quantile(AUCdiff, 0.05)) %>%
    dplyr::select(AUCdiffMedian, AUCdiffQ_025, AUCdiffQ_020, AUCdiffQ_010, AUCdiffQ_005, version) %>%
    ungroup() %>%
    arrange(AUCdiffMedian)

  temp.g %<>% left_join(temp.g.median, by = "version")

  temp.g.orig <- temp.g
  # vnucení řazení podle mediánu
  version_ordered <- with(droplevels(temp.g), reorder(version, AUCdiffMedian))
  temp.g <- droplevels(temp.g)
  temp.g$version <- factor(temp.g$version, levels = levels(version_ordered))
  # temp.g %<>% mutate(version = fct_reorder(version, AUCdiffMedian))


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(x = version, y = AUCdiff)) +
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
  ggsave(paste0(path.PP, "boxplot.test.", as.character(at), ".png"), width = 1500, height = 1000, units = "px")


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(occs.n, AUCdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.05) +
    geom_smooth(method = loess, size = 0.2) +
    theme_light() +
    theme(
      legend.position = "none", text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    ) +
    facet_wrap(~version)
  ggsave(paste0(path.PP, "trend.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  # porovnání přispění verzí per species
  title <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(AUC, with_ties = FALSE) %>% ungroup() %>% mutate(title = paste0(species, "\n", round(AUC_null, digits = 2), " -> ", round(AUC, digits = 2))) %>% dplyr::select(title)))
  names(title) <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(AUC, with_ties = FALSE) %>% ungroup() %>% dplyr::select(species)))
  ggplot(temp.g %>% ungroup(), aes(version, AUCdiff)) +
    geom_bar(stat = "identity", aes(fill = factor(version))) +
    theme_light() +
    theme(
      legend.position = "none", text = element_text(size = 4),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1)),
      panel.margin = unit(0.2, "lines")
    ) +
    facet_wrap(~species, labeller = as_labeller(title))
  ggsave(paste0(path.PP, "version-species.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")


  # přispění kombinací
  temp.g %<>% ungroup()

  counter <- 0
  first <- TRUE
  for (k in k6) {
    counter <- counter + 1
    print(counter)
    print(paste(k, collapse = "|"))
    temp <- temp.g %>%
      filter(version %in% k) %>%
      group_by(species) %>%
      slice_max(AUCdiff, with_ties = FALSE) %>%
      ungroup() %>%
      summarise(n = length(species), AUCdiffSum = sum(AUCdiff), mean = mean(AUCdiff), AUCdiffMedian = median(AUCdiff), AUCdiffQ_025 = quantile(AUCdiff, 0.25), AUCdiffQ_020 = quantile(AUCdiff, 0.20), AUCdiffQ_010 = quantile(AUCdiff, 0.10), AUCdiffQ_005 = quantile(AUCdiff, 0.05))
    temp$versionComb <- paste(k, collapse = "|")
    temp$versionCount <- length(k)
    if (first) {
      first <- FALSE
      res <- temp
    } else {
      res %<>% add_row(temp)
    }
  }
  combs[[as.character(at)]][["test"]] <- res %<>% arrange(desc(AUCdiffSum))
  write.csv(res, paste0(path.PP, "combs-test-", as.character(at), ".csv"), row.names = FALSE)
}

saveRDS(combs, paste0(path.PP, "combs.rds"))


summary_zasobnik.t <- t(as_tibble(summary_zasobnik))
rn <- row.names(summary_zasobnik.t)
cn <- names(summary_zasobnik[[1]])
summary_zasobnik.t <- as.data.frame(summary_zasobnik.t)
summary_zasobnik.t <- as_tibble(sapply(summary_zasobnik.t, as.numeric))
names(summary_zasobnik.t) <- cn
summary_zasobnik.t$treshold <- rn
summary_zasobnik.t$total <- 110
summary_zasobnik.t %<>% relocate(treshold) %>% relocate(total)

saveRDS(summary_zasobnik.t, paste0(path.PP, "version-species-treshold-count.rds"))
write.csv(summary_zasobnik.t, paste0(path.PP, "version-species-treshold-count.csv"), row.names = FALSE)


saveRDS(zasobnik0, paste0(path.PP, "version-species-treshold-count.rds"))




stop()
# #
# # rozdíl tresholdu podle AUC/auc.val.avg
# #
# tbl.all3 <- tbl.all
# tbl.all3 %<>% ungroup() %>%  group_by(species) %>%
#   slice_max(AUCdiff, with_ties = FALSE) %>%  filter(AUC >= 0.75)
# species.AUC <- tbl.all3$species
#





tbl.all3 <- tbl.all
tbl.all3 %<>% ungroup() %>%
  group_by(species) %>% # filter(cbi.val.sd <= 0.10) %>%
  # filter(cbi.val.avg >= 0.80) %>% # & cbi.val.sd <= 0.10
  # filter(auc.val.avg >= (max(auc.val.avg)  - 0.05))  %>%
  # filter(delta.AICc <= 2) %>%
  # filter(cbi.val.avg >= (max(cbi.val.avg)  - 0.05))  %>%
  # filter(auc.val.avg >= 0.60)  %>%
  # filter(delta.AICc <= 2) %>%
  # filter(or.10p.avg == min(or.10p.avg))  %>%
  # filter(delta.AICc <= 0.1) %>%
  # filter(delta.AICc == min(delta.AICc) ) %>%
  # filter(cbi.val.avg >= 0.75)
  # filter(cbi.val.avg >= max(cbi.val.avg))  %>%
  # filter(auc.val.avg >= 0.65)
  # filter(auc.val.avg >= max(auc.val.avg))
  #  filter(auc.val.avg >= 0.65 ) %>%
  # filter(cbi.val.avg == max(cbi.val.avg))   %>%
  # slice_max(cbi.val.avg, with_ties = FALSE) %>% filter(cbi.val.avg >= 0.70)
  slice_max(auc.val.avg, with_ties = FALSE) %>%
  filter(auc.val.avg >= 0.75)
# filter(auc.val.avg == max(auc.val.avg ) )  %>%
# filter(cbi.val.avg >= (max(cbi.val.avg)  - 0.10))  %>%
# filter(auc.val.avg  == max(auc.val.avg)) %>% filter(auc.val.avg  >= 0.65)
species.aucVal <- unique(tbl.all3$species)
naive.sp <- length(species.aucVal)
real.sp <- nrow(tbl.all3 %>% filter(AUC >= 0.75))
real2.sp <- nrow(tbl.all3 %>% filter(AUC >= 0.70))
print(naive.sp)
print(real.sp)
print(real2.sp)






# setdiff(species.AUC, species.aucVal)
# setdiff(species.aucVal,  species.AUC)




# přispění kombinací
tbl.all5 <- tbl.all.official
tbl.all5 %<>% ungroup()
# %>% filter(occs.n < 800)
selection <- selection[selection != "cz_un"] # bez
k6 <- comb_all(selection, 6)
# k6 <- comb_k(selection, 6)
counter <- 0
first <- TRUE
for (k in k6) {
  counter <- counter + 1
  print(counter)
  print(paste(k, collapse = "|"))
  temp <- tbl.all5 %>%
    filter(version %in% k) %>%
    group_by(species) %>%
    slice_max(AUCdiff, with_ties = FALSE) %>%
    ungroup() %>%
    summarise(n = length(species), AUCdiffSum = sum(AUCdiff), mean = mean(AUCdiff), AUCdiffMedian = median(AUCdiff), AUCdiffQ_025 = quantile(AUCdiff, 0.25), AUCdiffQ_020 = quantile(AUCdiff, 0.20), AUCdiffQ_010 = quantile(AUCdiff, 0.10), AUCdiffQ_005 = quantile(AUCdiff, 0.05))
  temp$versionComb <- paste(k, collapse = "|")
  temp$versionCount <- length(k)
  if (first) {
    first <- FALSE
    res <- temp
  } else {
    res %<>% add_row(temp)
  }
}
saveRDS(res, paste0(path.eval, "res.rds"))

#
##### # jen druhy co mají smysl... + bych měl možná i dofiltrovat null nad 50%?
#
#  tbl.all %<>% filter(AUC >= 0.7)
# jen druhy do 800 occs
#  tbl.all %<>% filter(occs.n < 800)

tbl.all <- tbl.all.official
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
ggsave(paste0(path.PP, "tbl.simple", ss, ".png"), width = 1500, height = 1000, units = "px")



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
# ggsave(paste0(path.PP, "tbl.png"), width = 1500, height = 1000, units = "px")

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
ggsave(paste0(path.PP, "tbl.trend-prevalence", ss, ".png"), width = 2000, height = 1500, units = "px")




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