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
# remove thin to get only TGOB derived versions
withThin <- FALSE

# generování ověření rozdílnosti verzí t.test-em (trvá extrémně dlouho, kešuju)
verified.generate <- FALSE
verified.top2 <- "verified.top2.rds"
verified.top1null <- "verified.top1null.rds"

if (file.exists(paste0(path.eval, verified.top2)) & file.exists(paste0(path.eval, verified.top1null))) {
  print("načítám existující verified.X")
  tmp.top2 <- readRDS(paste0(path.eval, verified.top2))
  tmp.top1null <- readRDS(paste0(path.eval, verified.top1null))
} else {
  print("budu generovat verified.X")
  verified.generate <- TRUE
}

#
# pokud je nějaká metoda korelovaná s LSD AUC, tak to znamená, že asi opravdu funguje korekce biasu? ne, jen že se lze spolehnout na výsledky auc.var.avg
#
# mě zajímá především korelace k nejlepším výsledkům, tyto totiž vybírám - dat do vztahu, jestli jsou vyšší AUC i lépe korelované
#


"%notin%" <- Negate("%in%")


modelsResults <- paste0(path.eval, "tbl.rep.rds")
if (file.exists(modelsResults)) {
  print("načítám existující tbl.rep")
  tbl <- tbl.rep <- readRDS(modelsResults)
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

  tbl <- tbl.rep <- as_tibble(out) %>% na.omit()
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

# selection <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "un")
# # selection.f <- c("ssos.topA100", "ssos.topS10", "tgob", "topA100", "topS10", "un", "ssos") # soos záměrně na konci
# selection.f2 <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "un")
# selection.rename <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "thin")

selection <- c(
  "TGOB", "ssosTGOB",
  "TSAO", "ssosTSAO",
  "TO", "ssosTO",
  "TS", "ssosTS",
  "un"
)
# selection.f <- c("ssos.topA100", "ssos.topS10", "tgob", "topA100", "topS10", "un", "ssos") # soos záměrně na konci
selection.f2 <- c(
  "TGOB", "ssosTGOB",
  "TSAO", "ssosTSAO",
  "TO", "ssosTO",
  "TS", "ssosTS",
  "un"
)
# selection.rename <- c(
#   "TGOB", "ssoTGOB",
#   "TSAO", "ssoTSAO",
#   "TO", "ssoTO",
#   "TS", "ssoTS",
#   "thin"
# )
selection.rename <- c(
  "TGOB", "TGOBsso",
  "TSAO", "TSAOsso",
  "TO", "TOsso",
  "TS", "TS.sso", # tečka dočasně kvůli řazení v legendě
  ".thin" # tečka dočasně kvůli řazení v legendě
)
selection.f2 <- paste0("^", selection.f2, "$")


names(selection.rename) <- selection.f2

pairs.compare <- list()
pairs.compare[[selection.rename[1]]] <- selection.rename[1:2]
pairs.compare[[selection.rename[3]]] <- selection.rename[3:4]
pairs.compare[[selection.rename[5]]] <- selection.rename[5:6]
pairs.compare[[selection.rename[7]]] <- selection.rename[7:8]


# # friendly
# # black + green / orange + yellow /  blue + skyblue / red + pink
# selection.colors <- c(
#   "#000000", "#009e73",
#   "#e69f00", "#f0e442",
#   "#0072b2", "#56b4e9",
#   "#d55e00", "#cc79a7",
#   "#e0e0e0"
# )

# black + green / orange + yellow /  blue + skyblue / red + pink
selection.colors <- c(
  "#1f2124", "#6a6e73",
  "#006414", "#5ccb5f",
  "#1465bb", "#81c9fa",
  "#bd0003", "#ff6c3e",
  "#e0e0e0"
)
# selection.colors <- paste0("#", selection.colors)
names(selection.colors) <- paste0("^", selection.rename, "$")

#
# dřívější verze s replikacemi pro t.test-y (rozdíl dvou nej verzí), neřeším null?
#

tbl.rep %<>% ungroup() %>% mutate(id = row_number())

tbl.rep.null.ids <- tbl.rep %>% filter(version == "un" & adjust == "0")
tbl.rep.null.ids.unique <- unique(tbl.rep.null.ids$id)

####
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

# přepsat počet thinning presencí, většinou nechci omezený výsledný součet presencí po thinningu
tbl.occs.n <- tbl %>%
  filter(id %in% tbl.null.ids.unique) %>%
  group_by(species) %>%
  slice_head(n = 1) %>%
  dplyr::select(species, occs.n)

tbl %<>% rename(occs.n.orig = occs.n)
tbl %<>% left_join(tbl.occs.n, by = c("species"))


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


# sjednotit i tabulku s replikacemi - nutné ještě před withThin!!!
tbl.rep.f <- tbl.rep %>%
  filter(version %in% selection) %>%
  ungroup() %>%
  mutate(version = str_replace_all(version, selection.rename)) %>%
  group_by(species, version) # final versions selection

tbl.rep.nn <- tbl.rep.f %>% filter(id %notin% tbl.rep.null.ids.unique) # not null
tbl.rep.null <- tbl.rep %>% filter(id %in% tbl.rep.null.ids.unique) # null - musím brát z původní tabulky vždy - ikdyž počítám s !withThin

if (!withThin) {
  # odstranit tinning z výsledků - čisté přispění TGOB verzí
  last <- length(selection)
  selection <- selection[-last]
  selection.f2 <- selection.f2[-last]
  selection.rename <- selection.rename[-last]
  path.PP <- gsub("/$", "_withoutThin/", path.PP)
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
  filter(version %in% selection) %>%
  ungroup() %>%
  mutate(version = str_replace_all(version, selection.rename)) %>%
  mutate(clr = str_replace_all(version, selection.colors)) %>%
  group_by(species, version) # final versions selection

selection.f2 <- selection.rename
k6 <- comb_all(selection.f2, length(selection.f2))

tbl.nn <- tbl.f %>% filter(id %notin% tbl.null.ids.unique) # not null
tbl.null <- tbl %>% filter(id %in% tbl.null.ids.unique) # null - musím brát z původní tabulky vždy - ikdyž počítám s !withThin

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
  # měl bych ale spočíst i omision rate nad LSD - tam by měl být smysluplný? + UPR+OPR
  # korelovanost výsledků auc AUC v rámci metod

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


  ################################################################################################################################################################################
  # kss auc.val.avgdiff
  ################################################################################################################################################################################


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


  if (at == 0 & verified.generate) {
    tmp.top2 <- list()
    tmp.top1null <- list()
    ### # + liší nejlepší od null? - nakonec to tam také budu potřebovat...
    # udělat zároveň - odebírat null idečka až tady? - jak se to zachová, když budu počítat s null verzí? - tu nedělat?
    # - nepočítat pokud je výsledek verze minusový, pak se mohou lišit, ale je to k ničemu...

    # musí jít spočíst jednou na začátku a pak už jen uložené využívat?? Cachovat...

    # přidat do výstupů vítěznou verzi a p-value

    # sum(unlist(tmp.top2[["test"]]), na.rm=TRUE)
    # sum(unlist(tmp.top1null[["test"]]), na.rm=TRUE)

    # tmp.top2.t <- t(as_tibble(tmp.top2))
    # rn <- row.names(tmp.top2.t)
    # cn <- names(tmp.top2[[1]])
    # tmp.top2.t <- as.data.frame(tmp.top2.t)
    # tmp.top2.t <- as_tibble(sapply(tmp.top2.t, as.numeric))
    # row.names(tmp.top2.t) <- rn
    # tmp.top2.t <- t(tmp.top2.t)
    # tmp.top2.t <- as_tibble(tmp.top2.t)
    # tmp.top2.t$species <- cn


    print("liší se výsledné AUC nejlepších dvou verzí (porovnává všechny replikace)? (další metriky???)  [auc.val.avg]")
    tmptbl <- temp.g %>%
      group_by(species) %>%
      arrange(desc(auc.val.avgdiff)) %>%
      slice_head(n = 2) %>%
      dplyr::select(species, version, adjust, tune.args)

    tmptbl.null <- tbl.null %>%
      group_by(species) %>%
      arrange(desc(auc.val.avgdiff)) %>%
      slice_head(n = 1) %>%
      dplyr::select(species, version, adjust, tune.args)


    for (sp in unique(tmptbl$species)) {
      print("")
      print(sp)

      #
      # top2
      #
      print("top2:")
      top2 <- tmptbl %>% filter(species == sp)

      top2v <- lapply(1:2, function(x) {
        tmp <- list()
        tmp[["t"]] <- tbl.rep.nn %>% filter(species == sp & version == top2[x, ]$version & adjust == top2[x, ]$adjust & tune.args == top2[x, ]$tune.args)
        tmp[["top"]] <- top2[x, ]$version

        # p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution, we can assume the normality.
        # kontrolu na alespoň 3 replikace musím udělat úplně na začátku!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (nrow(tmp[["t"]]) >= 3) {
          tmp[["shapiro.p"]] <- shapiro.test(tmp[["t"]]$auc.val.avg)$p.value
        } else {
          print(" < 3 replikací !!!")
          tmp[["shapiro.p"]] <- 0
        }
        return(tmp)
      })

      top2v.normality <- sum(sapply(top2v, function(x) x$shapiro.p > 0.05)) == 2

      # # F test: variances of the two groups are equal? netřeba, pokud nejsou stejné v t.test-u se provede Welch místo klasického
      # vartest <- var.test(top2v[[1]][["t"]]$auc.val.avg, top2v[[2]][["t"]]$auc.val.avg)   # p < 0.05 (Variances are not equal)
      if (top2v.normality) {
        print("OK (normality)")
        ttest <- t.test(top2v[[1]][["t"]]$auc.val.avg, top2v[[2]][["t"]]$auc.val.avg, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
        tmp.top2[["all"]][[sp]][["val"]] <- ttest$p.value
        tmp.top2[["all"]][[sp]][["val.auc1"]] <- mean(top2v[[1]][["t"]]$auc.val.avg)
        tmp.top2[["all"]][[sp]][["val.auc2"]] <- mean(top2v[[2]][["t"]]$auc.val.avg)
        tmp.top2[["all"]][[sp]][["val.version1"]] <- top2v[[1]][["top"]]
        tmp.top2[["all"]][[sp]][["val.version2"]] <- top2v[[2]][["top"]]
      } else {
        print("auc nejsou normálně rozložené nebo < 3 replikací, jiný test...")
        tmp.top2[["all"]][[sp]][["val"]] <- NA
        tmp.top2[["all"]][[sp]][["val.auc1"]] <- NA
        tmp.top2[["all"]][[sp]][["val.auc2"]] <- NA
        tmp.top2[["all"]][[sp]][["val.version1"]] <- NA
        tmp.top2[["all"]][[sp]][["val.version2"]] <- NA
      }

      #
      # top1 vs null
      #
      print("top1null:")
      top1null <- tmptbl.null %>% filter(species == sp)

      tmp.null.t <- tbl.rep.null %>% filter(species == sp & version == top1null[1, ]$version & adjust == top1null[1, ]$adjust & tune.args == top1null[1, ]$tune.args)

      if (nrow(tmp.null.t) >= 3) {
        tmp.null <- shapiro.test(tmp.null.t$auc.val.avg)$p.value
      } else {
        print("null < 3 replikací !!!")
        tmp.null <- 0
      }

      top1null.normality <- sum(c(tmp.null > 0.05, top2v[[1]]$shapiro.p > 0.05)) == 2

      if (top1null.normality) {
        print("OK (normality)")
        ttest <- t.test(top2v[[1]][["t"]]$auc.val.avg, tmp.null.t$auc.val.avg, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
        tmp.top1null[["all"]][[sp]][["val"]] <- ttest$p.value
        tmp.top1null[["all"]][[sp]][["val.auc1"]] <- mean(top2v[[1]][["t"]]$auc.val.avg)
        tmp.top1null[["all"]][[sp]][["val.auc2"]] <- mean(tmp.null.t$auc.val.avg)
        tmp.top1null[["all"]][[sp]][["val.version1"]] <- top2v[[1]][["top"]]
        tmp.top1null[["all"]][[sp]][["val.version2"]] <- "un"
      } else {
        print("auc (null a 1) nejsou normálně rozložené nebo < 3 replikací, jiný test...")
        tmp.top1null[["all"]][[sp]][["val"]] <- NA
        tmp.top1null[["all"]][[sp]][["val.auc1"]] <- NA
        tmp.top1null[["all"]][[sp]][["val.auc2"]] <- NA
        tmp.top1null[["all"]][[sp]][["val.version1"]] <- NA
        tmp.top1null[["all"]][[sp]][["val.version2"]] <- NA
      }
    }
  }
  ### ### ###
  ### ### ### startA val
  ### ### ###

  tobs <- ggplot() +
    geom_bar(data = temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avgdiff, with_ties = FALSE), stat = "identity", mapping = aes(species, auc.val.avgdiff, fill = factor(version)), color = "yellow", size = 0.01) +
    geom_point(data = temp.g %>% ungroup(), mapping = aes(species, auc.val.avgdiff, fill = factor(version)), size = 1, stroke = 0.1, color = "yellow", shape = 21, alpha = 0.90) +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      text = element_text(size = 6),
      axis.text.y = element_markdown(),
      # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )

  tt.null <- sapply(tmp.top1null[["all"]], function(x) x$val < 0.05)
  tt <- sapply(tmp.top2[["all"]], function(x) x$val < 0.05)

  # základní + změna řazení
  no <- temp.g %>%
    ungroup() %>%
    group_by(species) %>%
    slice_max(auc.val.avgdiff, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(tt.null = tt.null) %>%
    mutate(tt = tt) %>%
    mutate(title = paste0(ifelse(!is.na(tt.null) & tt.null & !is.na(tt) & tt, paste0("***", species, "***"), ifelse(!is.na(tt.null) & tt.null, paste0("**", species, "**"), ifelse(!is.na(tt) & tt, paste0("*", species, "*"), species))), " | ", formatC(auc.val.avg_null, digits = 2, format = "f"), "->", formatC(auc.val.avg, digits = 2, format = "f"))) %>%
    dplyr::select(title, occs.n, species, auc.val.avgdiff, auc.val.avg, auc.val.avg_null)

  title <- unname(unlist(no$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no$occs.n))))

  tobs + scale_x_discrete(labels = rev(title.occs.n), limits = rev) + xlab("species (ordered by: alphabet); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species.val.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")


  # occs.n
  no.temp <- no %>%
    arrange(occs.n)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: sum of occupied squares); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species-orderByOccs.val.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  # auc.val.avgdiff
  no.temp <- no %>%
    arrange(auc.val.avgdiff)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avgdiff); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species-orderByDiff.val.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  # auc.val.avg
  no.temp <- no %>%
    arrange(auc.val.avg)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avg); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species-orderByAuc.val.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  # auc.val.avg_null
  no.temp <- no %>%
    arrange(auc.val.avg_null)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avg_null); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species-orderByAucNull.val.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  ### ### ### endB

  #
  # páry
  #
  for (cmp in pairs.compare) {
    temp.g.c <- temp.g


    temp.g.c %<>% filter(version %in% cmp)

    if (at == 0 & verified.generate) {
      print("liší se výsledné AUC párů verzí (porovnává všechny replikace)? (další metriky???)  [auc.val.avg]")
      print(cmp[1])
      tmptbl <- temp.g.c %>%
        group_by(species) %>%
        arrange(desc(auc.val.avgdiff)) %>%
        slice_head(n = 2) %>%
        dplyr::select(species, version, adjust, tune.args)

      tmptbl.null <- tbl.null %>%
        group_by(species) %>%
        arrange(desc(auc.val.avgdiff)) %>%
        slice_head(n = 1) %>%
        dplyr::select(species, version, adjust, tune.args)


      for (sp in unique(tmptbl$species)) {
        print("")
        print(sp)

        #
        # top2
        #
        print("top2:")
        top2 <- tmptbl %>% filter(species == sp)

        top2v <- lapply(1:2, function(x) {
          tmp <- list()
          tmp[["t"]] <- tbl.rep.nn %>% filter(species == sp & version == top2[x, ]$version & adjust == top2[x, ]$adjust & tune.args == top2[x, ]$tune.args)
          tmp[["top"]] <- top2[x, ]$version

          # p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution, we can assume the normality.
          # kontrolu na alespoň 3 replikace musím udělat úplně na začátku!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (nrow(tmp[["t"]]) >= 3) {
            tmp[["shapiro.p"]] <- shapiro.test(tmp[["t"]]$auc.val.avg)$p.value
          } else {
            print(" < 3 replikací !!!")
            tmp[["shapiro.p"]] <- 0
          }
          return(tmp)
        })

        top2v.normality <- sum(sapply(top2v, function(x) x$shapiro.p > 0.05)) == 2

        # # F test: variances of the two groups are equal? netřeba, pokud nejsou stejné v t.test-u se provede Welch místo klasického
        # vartest <- var.test(top2v[[1]][["t"]]$auc.val.avg, top2v[[2]][["t"]]$auc.val.avg)   # p < 0.05 (Variances are not equal)
        if (top2v.normality) {
          print("OK (normality)")
          ttest <- t.test(top2v[[1]][["t"]]$auc.val.avg, top2v[[2]][["t"]]$auc.val.avg, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
          tmp.top2[[cmp[1]]][[sp]][["val"]] <- ttest$p.value
          tmp.top2[[cmp[1]]][[sp]][["val.auc1"]] <- mean(top2v[[1]][["t"]]$auc.val.avg)
          tmp.top2[[cmp[1]]][[sp]][["val.auc2"]] <- mean(top2v[[2]][["t"]]$auc.val.avg)
          tmp.top2[[cmp[1]]][[sp]][["val.version1"]] <- top2v[[1]][["top"]]
          tmp.top2[[cmp[1]]][[sp]][["val.version2"]] <- top2v[[2]][["top"]]
        } else {
          print("auc nejsou normálně rozložené nebo < 3 replikací, jiný test...")
          tmp.top2[[cmp[1]]][[sp]][["val"]] <- NA
          tmp.top2[[cmp[1]]][[sp]][["val.auc1"]] <- NA
          tmp.top2[[cmp[1]]][[sp]][["val.auc2"]] <- NA
          tmp.top2[[cmp[1]]][[sp]][["val.version1"]] <- NA
          tmp.top2[[cmp[1]]][[sp]][["val.version2"]] <- NA
        }

        #
        # top1 vs null
        #
        print("top1null:")
        top1null <- tmptbl.null %>% filter(species == sp)

        tmp.null.t <- tbl.rep.null %>% filter(species == sp & version == top1null[1, ]$version & adjust == top1null[1, ]$adjust & tune.args == top1null[1, ]$tune.args)

        if (nrow(tmp.null.t) >= 3) {
          tmp.null <- shapiro.test(tmp.null.t$auc.val.avg)$p.value
        } else {
          print("null < 3 replikací !!!")
          tmp.null <- 0
        }

        top1null.normality <- sum(c(tmp.null > 0.05, top2v[[1]]$shapiro.p > 0.05)) == 2

        if (top1null.normality) {
          print("OK (normality)")
          ttest <- t.test(top2v[[1]][["t"]]$auc.val.avg, tmp.null.t$auc.val.avg, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
          tmp.top1null[[cmp[1]]][[sp]][["val"]] <- ttest$p.value
          tmp.top1null[[cmp[1]]][[sp]][["val.auc1"]] <- mean(top2v[[1]][["t"]]$auc.val.avg)
          tmp.top1null[[cmp[1]]][[sp]][["val.auc2"]] <- mean(tmp.null.t$auc.val.avg)
          tmp.top1null[[cmp[1]]][[sp]][["val.version1"]] <- top2v[[1]][["top"]]
          tmp.top1null[[cmp[1]]][[sp]][["val.version2"]] <- "un"
        } else {
          print("auc (null a 1) nejsou normálně rozložené nebo < 3 replikací, jiný test...")
          tmp.top1null[[cmp[1]]][[sp]][["val"]] <- NA
          tmp.top1null[[cmp[1]]][[sp]][["val.auc1"]] <- NA
          tmp.top1null[[cmp[1]]][[sp]][["val.auc2"]] <- NA
          tmp.top1null[[cmp[1]]][[sp]][["val.version1"]] <- NA
          tmp.top1null[[cmp[1]]][[sp]][["val.version2"]] <- NA
        }
      }
    }

    ### ### ###
    ### ### ### startA val
    ### ### ###

    tobs <- ggplot() +
      geom_bar(data = temp.g.c %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avgdiff, with_ties = FALSE), stat = "identity", mapping = aes(species, auc.val.avgdiff, fill = factor(version)), color = "yellow", size = 0.01) +
      geom_point(data = temp.g.c %>% ungroup(), mapping = aes(species, auc.val.avgdiff, fill = factor(version)), size = 1, stroke = 0.1, color = "yellow", shape = 21, alpha = 0.90) +
      scale_fill_manual(values = unique(unname(unlist(temp.g.c %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
      theme_light() +
      theme(
        text = element_text(size = 6),
        axis.text.y = element_markdown(),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
        panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
        strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
      )

    tt.null <- sapply(tmp.top1null[[cmp[1]]], function(x) x$val < 0.05)
    tt <- sapply(tmp.top2[[cmp[1]]], function(x) x$val < 0.05)

    # základní + změna řazení
    no <- temp.g.c %>%
      ungroup() %>%
      group_by(species) %>%
      slice_max(auc.val.avgdiff, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(tt.null = tt.null) %>%
      mutate(tt = tt) %>%
      mutate(title = paste0(ifelse(!is.na(tt.null) & tt.null & !is.na(tt) & tt, paste0("***", species, "***"), ifelse(!is.na(tt.null) & tt.null, paste0("**", species, "**"), ifelse(!is.na(tt) & tt, paste0("*", species, "*"), species))), " | ", formatC(auc.val.avg_null, digits = 2, format = "f"), "->", formatC(auc.val.avg, digits = 2, format = "f"))) %>%
      dplyr::select(title, occs.n, species, auc.val.avgdiff, auc.val.avg, auc.val.avg_null)

    title <- unname(unlist(no$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no$occs.n))))

    tobs + scale_x_discrete(labels = rev(title.occs.n), limits = rev) + xlab("species (ordered by: alphabet); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species.val.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")


    # occs.n
    no.temp <- no %>%
      arrange(occs.n)
    order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
    title <- unname(unlist(no.temp$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

    tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: sum of occupied squares); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species-orderByOccs.val.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    # auc.val.avgdiff
    no.temp <- no %>%
      arrange(auc.val.avgdiff)
    order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
    title <- unname(unlist(no.temp$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

    tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avgdiff); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species-orderByDiff.val.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    # auc.val.avg
    no.temp <- no %>%
      arrange(auc.val.avg)
    order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
    title <- unname(unlist(no.temp$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

    tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avg); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species-orderByAuc.val.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    # auc.val.avg_null
    no.temp <- no %>%
      arrange(auc.val.avg_null)
    order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
    title <- unname(unlist(no.temp$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

    tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avg_null); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species-orderByAucNull.val.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    ### ### ### endB
  }


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(occs.n, auc.val.avgdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.5) +
    geom_smooth(method = loess, size = 0.2) +
    scale_color_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      # legend.position = "none",
      text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    ) +
    facet_wrap(~version)
  ggsave(paste0(path.PP, "trend.val.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")



  ### ### ### start

  ## nesmím  nechávat thin verzi nové occs.n pořadí!! - přepsat jinou verzí - obecně všude?

  # a) vše - OVERALL
  ggplot(temp.g %>% ungroup(), aes(occs.n, auc.val.avgdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.05) +
    geom_smooth(method = loess, size = 0.2) +
    scale_color_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )
  ggsave(paste0(path.PP, "trend-overall.val.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  # a) vše - OVERALL - best
  ggplot(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avgdiff, with_ties = FALSE), aes(occs.n, auc.val.avgdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.5) +
    geom_smooth(method = loess, size = 0.2) +
    scale_color_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )
  ggsave(paste0(path.PP, "trend-overall-best.val.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  ### ### ### end









  # porovnání přispění verzí per species
  title <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avg, with_ties = FALSE) %>% ungroup() %>% mutate(title = paste0(species, "\n", formatC(auc.val.avg_null, digits = 2, format = "f"), " -> ", formatC(auc.val.avg, digits = 2, format = "f"))) %>% dplyr::select(title)))
  names(title) <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avg, with_ties = FALSE) %>% ungroup() %>% dplyr::select(species)))

  ggplot(temp.g %>% ungroup(), aes(version, auc.val.avgdiff)) +
    geom_bar(stat = "identity", aes(fill = factor(version))) +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
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




  # boxploty - změna řazení

  temp.g.orig <- temp.g
  # vnucení řazení podle mediánu
  version_ordered <- with(droplevels(temp.g), reorder(version, AUCdiffMedian))
  temp.g <- droplevels(temp.g)
  temp.g$version <- factor(temp.g$version, levels = levels(version_ordered))
  # temp.g %<>% mutate(version = fct_reorder(version, AUCdiffMedian))


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(x = version, y = auc.val.avgdiff, fill = clr)) +
    # stat_summary(fun.data=boxplotCustom, geom="boxplot", lwd=0.1,notch=TRUE)  +
    geom_boxplot(size = 0.1, notch = TRUE, outlier.size = 0.1, outlier.stroke = 0.3) +
    geom_point(aes(y = AUCdiffQ_010), size = 1, shape = 3, stroke = 0.01) +
    geom_point(aes(y = AUCdiffQ_005), size = 1, shape = 4, stroke = 0.01) +
    theme_light() +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(clr)) %>% dplyr::select(clr))))) +
    theme(
      legend.position = "none", text = element_text(size = 4),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  ggsave(paste0(path.PP, "boxplot.val.", as.character(at), ".png"), width = 1500, height = 1000, units = "px")



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


  ################################################################################################################################################################################
  # kss auc.val.avgdiff test
  ################################################################################################################################################################################

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
  ggplot(temp.g %>% ungroup(), aes(x = version, y = auc.val.avgdiff, fill = clr)) +
    # stat_summary(fun.data=boxplotCustom, geom="boxplot", lwd=0.1,notch=TRUE)  +
    geom_boxplot(size = 0.1, notch = TRUE, outlier.size = 0.1, outlier.stroke = 0.3) +
    geom_point(aes(y = AUCdiffQ_010), size = 1, shape = 3, stroke = 0.01) +
    geom_point(aes(y = AUCdiffQ_005), size = 1, shape = 4, stroke = 0.01) +
    theme_light() +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(clr)) %>% dplyr::select(clr))))) +
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
  title <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avg, with_ties = FALSE) %>% ungroup() %>% mutate(title = paste0(species, "\n", formatC(auc.val.avg_null, digits = 2, format = "f"), " -> ", formatC(auc.val.avg, digits = 2, format = "f"))) %>% dplyr::select(title)))
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




  ################################################################################################################################################################################
  # kss AUCdiff
  ################################################################################################################################################################################

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


  if (at == 0 & verified.generate) {
    print("liší se výsledné AUC nejlepších dvou verzí (porovnává všechny replikace)? (další metriky???)  [AUC]")
    tmptbl <- temp.g %>%
      group_by(species) %>%
      arrange(desc(AUCdiff)) %>%
      slice_head(n = 2) %>%
      dplyr::select(species, version, adjust, tune.args)

    tmptbl.null <- tbl.null %>%
      group_by(species) %>%
      arrange(desc(AUCdiff)) %>%
      slice_head(n = 1) %>%
      dplyr::select(species, version, adjust, tune.args)


    for (sp in unique(tmptbl$species)) {
      print("")
      print(sp)

      #
      # top2
      #
      print("top2:")
      top2 <- tmptbl %>% filter(species == sp)

      top2v <- lapply(1:2, function(x) {
        tmp <- list()
        tmp[["t"]] <- tbl.rep.nn %>% filter(species == sp & version == top2[x, ]$version & adjust == top2[x, ]$adjust & tune.args == top2[x, ]$tune.args)
        tmp[["top"]] <- top2[x, ]$version

        # p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution, we can assume the normality.
        # kontrolu na alespoň 3 replikace musím udělat úplně na začátku!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (nrow(tmp[["t"]]) >= 3) {
          tmp[["shapiro.p"]] <- shapiro.test(tmp[["t"]]$AUC)$p.value
        } else {
          print(" < 3 replikací !!!")
          tmp[["shapiro.p"]] <- 0
        }
        return(tmp)
      })

      top2v.normality <- sum(sapply(top2v, function(x) x$shapiro.p > 0.05)) == 2

      # # F test: variances of the two groups are equal? netřeba, pokud nejsou stejné v t.test-u se provede Welch místo klasického
      # vartest <- var.test(top2v[[1]][["t"]]$AUC, top2v[[2]][["t"]]$AUC)   # p < 0.05 (Variances are not equal)
      if (top2v.normality) {
        print("OK (normality)")
        ttest <- t.test(top2v[[1]][["t"]]$AUC, top2v[[2]][["t"]]$AUC, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
        tmp.top2[["all"]][[sp]][["test"]] <- ttest$p.value
        tmp.top2[["all"]][[sp]][["test.auc1"]] <- mean(top2v[[1]][["t"]]$AUC)
        tmp.top2[["all"]][[sp]][["test.auc2"]] <- mean(top2v[[2]][["t"]]$AUC)
        tmp.top2[["all"]][[sp]][["test.version1"]] <- top2v[[1]][["top"]]
        tmp.top2[["all"]][[sp]][["test.version2"]] <- top2v[[2]][["top"]]
      } else {
        print("auc nejsou normálně rozložené nebo < 3 replikací, jiný test...")
        tmp.top2[["all"]][[sp]][["test"]] <- NA
        tmp.top2[["all"]][[sp]][["test.auc1"]] <- NA
        tmp.top2[["all"]][[sp]][["test.auc2"]] <- NA
        tmp.top2[["all"]][[sp]][["test.version1"]] <- NA
        tmp.top2[["all"]][[sp]][["test.version2"]] <- NA
      }

      #
      # top1 vs null
      #
      print("top1null:")
      top1null <- tmptbl.null %>% filter(species == sp)

      tmp.null.t <- tbl.rep.null %>% filter(species == sp & version == top1null[1, ]$version & adjust == top1null[1, ]$adjust & tune.args == top1null[1, ]$tune.args)

      if (nrow(tmp.null.t) >= 3) {
        tmp.null <- shapiro.test(tmp.null.t$AUC)$p.value
      } else {
        print("null < 3 replikací !!!")
        tmp.null <- 0
      }

      top1null.normality <- sum(c(tmp.null > 0.05, top2v[[1]]$shapiro.p > 0.05)) == 2

      if (top1null.normality) {
        print("OK (normality)")
        ttest <- t.test(top2v[[1]][["t"]]$AUC, tmp.null.t$AUC, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
        tmp.top1null[["all"]][[sp]][["test"]] <- ttest$p.value
        tmp.top1null[["all"]][[sp]][["test.auc1"]] <- mean(top2v[[1]][["t"]]$AUC)
        tmp.top1null[["all"]][[sp]][["test.auc2"]] <- mean(tmp.null.t$AUC)
        tmp.top1null[["all"]][[sp]][["test.version1"]] <- top2v[[1]][["top"]]
        tmp.top1null[["all"]][[sp]][["test.version2"]] <- "un"
      } else {
        print("auc (null a 1) nejsou normálně rozložené nebo < 3 replikací, jiný test...")
        tmp.top1null[["all"]][[sp]][["test"]] <- NA
        tmp.top1null[["all"]][[sp]][["test.auc1"]] <- NA
        tmp.top1null[["all"]][[sp]][["test.auc2"]] <- NA
        tmp.top1null[["all"]][[sp]][["test.version1"]] <- NA
        tmp.top1null[["all"]][[sp]][["test.version2"]] <- NA
      }
    }
  }

  ### ### ###
  ### ### ### startA TEST
  ### ### ###

  tobs <- ggplot() +
    geom_bar(data = temp.g %>% ungroup() %>% group_by(species) %>% slice_max(AUCdiff, with_ties = FALSE), stat = "identity", mapping = aes(species, AUCdiff, fill = factor(version)), color = "yellow", size = 0.01) +
    geom_point(data = temp.g %>% ungroup(), mapping = aes(species, AUCdiff, fill = factor(version)), size = 1, stroke = 0.1, color = "yellow", shape = 21, alpha = 0.90) +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      text = element_text(size = 6),
      axis.text.y = element_markdown(),
      # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )

  tt.null <- sapply(tmp.top1null[["all"]], function(x) x$test < 0.05)
  tt <- sapply(tmp.top2[["all"]], function(x) x$test < 0.05)

  # základní + změna řazení
  no <- temp.g %>%
    ungroup() %>%
    group_by(species) %>%
    slice_max(AUCdiff, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(tt.null = tt.null) %>%
    mutate(tt = tt) %>%
    mutate(title = paste0(ifelse(!is.na(tt.null) & tt.null & !is.na(tt) & tt, paste0("***", species, "***"), ifelse(!is.na(tt.null) & tt.null, paste0("**", species, "**"), ifelse(!is.na(tt) & tt, paste0("*", species, "*"), species))), " | ", formatC(AUC_null, digits = 2, format = "f"), "->", formatC(AUC, digits = 2, format = "f"))) %>%
    dplyr::select(title, occs.n, species, AUCdiff, AUC, AUC_null)

  title <- unname(unlist(no$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no$occs.n))))

  tobs + scale_x_discrete(labels = rev(title.occs.n), limits = rev) + xlab("species (ordered by: alphabet); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species.test.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  # occs.n
  no.temp <- no %>%
    arrange(occs.n)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: sum of occupied squares); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species-orderByOccs.test.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  # AUCdiff
  no.temp <- no %>%
    arrange(AUCdiff)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUCdiff); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species-orderByDiff.test.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  # AUC
  no.temp <- no %>%
    arrange(AUC)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUC); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species-orderByAuc.test.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  # AUC_null
  no.temp <- no %>%
    arrange(AUC_null)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUC_null); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP, "trend-overall-best-species-orderByAucNull.test.", as.character(at), ".png"), width = 1500, height = 2000, units = "px")

  ### ### ### endB

  #
  # páry
  #
  for (cmp in pairs.compare) {
    temp.g.c <- temp.g


    temp.g.c %<>% filter(version %in% cmp)

    if (at == 0 & verified.generate) {
      print("liší se výsledné AUC párů verzí (porovnává všechny replikace)? (další metriky???)  [AUC]")
      print(cmp[1])
      tmptbl <- temp.g.c %>%
        group_by(species) %>%
        arrange(desc(AUCdiff)) %>%
        slice_head(n = 2) %>%
        dplyr::select(species, version, adjust, tune.args)

      tmptbl.null <- tbl.null %>%
        group_by(species) %>%
        arrange(desc(AUCdiff)) %>%
        slice_head(n = 1) %>%
        dplyr::select(species, version, adjust, tune.args)


      for (sp in unique(tmptbl$species)) {
        print("")
        print(sp)

        #
        # top2
        #
        print("top2:")
        top2 <- tmptbl %>% filter(species == sp)

        top2v <- lapply(1:2, function(x) {
          tmp <- list()
          tmp[["t"]] <- tbl.rep.nn %>% filter(species == sp & version == top2[x, ]$version & adjust == top2[x, ]$adjust & tune.args == top2[x, ]$tune.args)
          tmp[["top"]] <- top2[x, ]$version

          # p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution, we can assume the normality.
          # kontrolu na alespoň 3 replikace musím udělat úplně na začátku!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (nrow(tmp[["t"]]) >= 3) {
            tmp[["shapiro.p"]] <- shapiro.test(tmp[["t"]]$AUC)$p.value
          } else {
            print(" < 3 replikací !!!")
            tmp[["shapiro.p"]] <- 0
          }
          return(tmp)
        })

        top2v.normality <- sum(sapply(top2v, function(x) x$shapiro.p > 0.05)) == 2

        # # F test: variances of the two groups are equal? netřeba, pokud nejsou stejné v t.test-u se provede Welch místo klasického
        # vartest <- var.test(top2v[[1]][["t"]]$AUC, top2v[[2]][["t"]]$AUC)   # p < 0.05 (Variances are not equal)
        if (top2v.normality) {
          print("OK (normality)")
          ttest <- t.test(top2v[[1]][["t"]]$AUC, top2v[[2]][["t"]]$AUC, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
          tmp.top2[[cmp[1]]][[sp]][["test"]] <- ttest$p.value
          tmp.top2[[cmp[1]]][[sp]][["test.auc1"]] <- mean(top2v[[1]][["t"]]$AUC)
          tmp.top2[[cmp[1]]][[sp]][["test.auc2"]] <- mean(top2v[[2]][["t"]]$AUC)
          tmp.top2[[cmp[1]]][[sp]][["test.version1"]] <- top2v[[1]][["top"]]
          tmp.top2[[cmp[1]]][[sp]][["test.version2"]] <- top2v[[2]][["top"]]
        } else {
          print("auc nejsou normálně rozložené nebo < 3 replikací, jiný test...")
          tmp.top2[[cmp[1]]][[sp]][["test"]] <- NA
          tmp.top2[[cmp[1]]][[sp]][["test.auc1"]] <- NA
          tmp.top2[[cmp[1]]][[sp]][["test.auc2"]] <- NA
          tmp.top2[[cmp[1]]][[sp]][["test.version1"]] <- NA
          tmp.top2[[cmp[1]]][[sp]][["test.version2"]] <- NA
        }

        #
        # top1 vs null
        #
        print("top1null:")
        top1null <- tmptbl.null %>% filter(species == sp)

        tmp.null.t <- tbl.rep.null %>% filter(species == sp & version == top1null[1, ]$version & adjust == top1null[1, ]$adjust & tune.args == top1null[1, ]$tune.args)

        if (nrow(tmp.null.t) >= 3) {
          tmp.null <- shapiro.test(tmp.null.t$AUC)$p.value
        } else {
          print("null < 3 replikací !!!")
          tmp.null <- 0
        }

        top1null.normality <- sum(c(tmp.null > 0.05, top2v[[1]]$shapiro.p > 0.05)) == 2

        if (top1null.normality) {
          print("OK (normality)")
          ttest <- t.test(top2v[[1]][["t"]]$AUC, tmp.null.t$AUC, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
          tmp.top1null[[cmp[1]]][[sp]][["test"]] <- ttest$p.value
          tmp.top1null[[cmp[1]]][[sp]][["test.auc1"]] <- mean(top2v[[1]][["t"]]$AUC)
          tmp.top1null[[cmp[1]]][[sp]][["test.auc2"]] <- mean(tmp.null.t$AUC)
          tmp.top1null[[cmp[1]]][[sp]][["test.version1"]] <- top2v[[1]][["top"]]
          tmp.top1null[[cmp[1]]][[sp]][["test.version2"]] <- "un"
        } else {
          print("auc (null a 1) nejsou normálně rozložené nebo < 3 replikací, jiný test...")
          tmp.top1null[[cmp[1]]][[sp]][["test"]] <- NA
          tmp.top1null[[cmp[1]]][[sp]][["test.auc1"]] <- NA
          tmp.top1null[[cmp[1]]][[sp]][["test.auc2"]] <- NA
          tmp.top1null[[cmp[1]]][[sp]][["test.version1"]] <- NA
          tmp.top1null[[cmp[1]]][[sp]][["test.version2"]] <- NA
        }
      }
    }

    ### ### ###
    ### ### ### startA TEST
    ### ### ###

    tobs <- ggplot() +
      geom_bar(data = temp.g.c %>% ungroup() %>% group_by(species) %>% slice_max(AUCdiff, with_ties = FALSE), stat = "identity", mapping = aes(species, AUCdiff, fill = factor(version)), color = "yellow", size = 0.01) +
      geom_point(data = temp.g.c %>% ungroup(), mapping = aes(species, AUCdiff, fill = factor(version)), size = 1, stroke = 0.1, color = "yellow", shape = 21, alpha = 0.90) +
      scale_fill_manual(values = unique(unname(unlist(temp.g.c %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
      theme_light() +
      theme(
        text = element_text(size = 6),
        axis.text.y = element_markdown(),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
        panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
        strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
      )

    tt.null <- sapply(tmp.top1null[[cmp[1]]], function(x) x$test < 0.05)
    tt <- sapply(tmp.top2[[cmp[1]]], function(x) x$test < 0.05)

    # základní + změna řazení
    no <- temp.g.c %>%
      ungroup() %>%
      group_by(species) %>%
      slice_max(AUCdiff, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(tt.null = tt.null) %>%
      mutate(tt = tt) %>%
      mutate(title = paste0(ifelse(!is.na(tt.null) & tt.null & !is.na(tt) & tt, paste0("***", species, "***"), ifelse(!is.na(tt.null) & tt.null, paste0("**", species, "**"), ifelse(!is.na(tt) & tt, paste0("*", species, "*"), species))), " | ", formatC(AUC_null, digits = 2, format = "f"), "->", formatC(AUC, digits = 2, format = "f"))) %>%
      dplyr::select(title, occs.n, species, AUCdiff, AUC, AUC_null)

    title <- unname(unlist(no$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no$occs.n))))

    tobs + scale_x_discrete(labels = rev(title.occs.n), limits = rev) + xlab("species (ordered by: alphabet); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species.test.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    # occs.n
    no.temp <- no %>%
      arrange(occs.n)
    order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
    title <- unname(unlist(no.temp$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

    tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: sum of occupied squares); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species-orderByOccs.test.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    # AUCdiff
    no.temp <- no %>%
      arrange(AUCdiff)
    order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
    title <- unname(unlist(no.temp$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

    tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUCdiff); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species-orderByDiff.test.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    # AUC
    no.temp <- no %>%
      arrange(AUC)
    order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
    title <- unname(unlist(no.temp$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

    tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUC); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species-orderByAuc.test.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    # AUC_null
    no.temp <- no %>%
      arrange(AUC_null)
    order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
    title <- unname(unlist(no.temp$title))
    title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

    tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUC_null); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
    ggsave(paste0(path.PP, "trend-overall-best-species-orderByAucNull.test.", as.character(at), "---", cmp[1], ".png"), width = 1500, height = 2000, units = "px")

    ### ### ### endB
  }


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(occs.n, AUCdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.5) +
    geom_smooth(method = loess, size = 0.2) +
    scale_color_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      # legend.position = "none",
      text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    ) +
    facet_wrap(~version)
  ggsave(paste0(path.PP, "trend.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")


  ### ### ### start

  ## nesmím  nechávat thin verzi nové occs.n pořadí!! - přepsat jinou verzí - obecně všude?

  # a) vše - OVERALL
  ggplot(temp.g %>% ungroup(), aes(occs.n, AUCdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.05) +
    geom_smooth(method = loess, size = 0.2) +
    scale_color_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )
  ggsave(paste0(path.PP, "trend-overall.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  # a) vše - OVERALL - best
  ggplot(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(AUCdiff, with_ties = FALSE), aes(occs.n, AUCdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.5) +
    geom_smooth(method = loess, size = 0.2) +
    scale_color_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )
  ggsave(paste0(path.PP, "trend-overall-best.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")



  ### ### ### end





  # porovnání přispění verzí per species
  title <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(AUC, with_ties = FALSE) %>% ungroup() %>% mutate(title = paste0(species, "\n", formatC(AUC_null, digits = 2, format = "f"), " -> ", formatC(AUC, digits = 2, format = "f"))) %>% dplyr::select(title)))
  names(title) <- unname(unlist(temp.g %>% ungroup() %>% group_by(species) %>% slice_max(AUC, with_ties = FALSE) %>% ungroup() %>% dplyr::select(species)))
  ggplot(temp.g %>% ungroup(), aes(version, AUCdiff)) +
    geom_bar(stat = "identity", aes(fill = factor(version))) +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
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





  # boxploty - změna řazení
  temp.g.orig <- temp.g
  # vnucení řazení podle mediánu
  version_ordered <- with(droplevels(temp.g), reorder(version, AUCdiffMedian))
  temp.g <- droplevels(temp.g)
  temp.g$version <- factor(temp.g$version, levels = levels(version_ordered))
  # temp.g %<>% mutate(version = fct_reorder(version, AUCdiffMedian))


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(x = version, y = AUCdiff, fill = clr)) +
    # stat_summary(fun.data=boxplotCustom, geom="boxplot", lwd=0.1,notch=TRUE)  +
    geom_boxplot(size = 0.1, notch = TRUE, outlier.size = 0.1, outlier.stroke = 0.3) +
    geom_point(aes(y = AUCdiffQ_010), size = 1, shape = 3, stroke = 0.01) +
    geom_point(aes(y = AUCdiffQ_005), size = 1, shape = 4, stroke = 0.01) +
    theme_light() +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(clr)) %>% dplyr::select(clr))))) +
    theme(
      legend.position = "none", text = element_text(size = 4),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  ggsave(paste0(path.PP, "boxplot.test.", as.character(at), ".png"), width = 1500, height = 1000, units = "px")




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

#
# souhrny
#

pairs.main <- names(pairs.compare)
pairs.sso <- sapply(pairs.compare, function(x) paste(x, collapse = "|"))

for (v in c("val", "test")) {
  # párový přínost sso verzí
  combs.select <- c(pairs.main, pairs.sso)
  combs.select.res <- combs[["0"]][[v]] %>%
    filter(versionComb %in% combs.select) %>%
    arrange(versionComb) %>%
    dplyr::select(AUCdiffSum, mean, AUCdiffMedian, versionComb)
  # manual reorder
  combs.select.res <- combs.select.res[c(1:5, 8, 4, 7), ]
  write.csv(combs.select.res, paste0(path.PP, "pair-sso-version-cumsum.", v, ".csv"), row.names = FALSE)

  # TGOB + další řádné + sso
  tgob.main.res <- combs[["0"]][[v]] %>%
    filter(versionComb %in% c(pairs.main[1], paste(pairs.main, collapse = "|"))) %>%
    arrange(versionComb) %>%
    dplyr::select(AUCdiffSum, mean, AUCdiffMedian, versionComb)
  tgob.main.res %<>% add_row(combs[["0"]][[v]] %>% slice_head(n = 1) %>% dplyr::select(AUCdiffSum, mean, AUCdiffMedian, versionComb))
  write.csv(tgob.main.res, paste0(path.PP, "main-sso-version-cumsum.", v, ".csv"), row.names = FALSE)
}


if (verified.generate) {
  saveRDS(tmp.top2, paste0(path.eval, verified.top2))
  saveRDS(tmp.top1null, paste0(path.eval, verified.top1null))
}

end_time <- Sys.time()
print(end_time - start_time)