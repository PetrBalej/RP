start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "ENMeval", "sdmATM", "ggplot2", "grid", "gridExtra", "ggtext", "rmarkdown") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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

se <- function(x) sqrt(var(x) / length(x))

# select null if performs better
withNull <- FALSE
# remove thin to get only TGOB derived versions
withThin <- FALSE

# tgobXostatní
anyAlt <- FALSE

minReps <- 5

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
      path.tgob,
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

    # # dočasné: vzniká mi tu v kombinaci un/0 ve version2 NA které pak níže neprojdou přes na.omit(), přepsat
    rds %<>% mutate(version2 = ifelse(is.na(version2) & version == "un" & adjust == 0, "un", version2))
    #  rds %<>% dplyr::select(-version2)

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


#
# asi už nemůžu párovat (groupovat) podle adjustu a version - sjednocením vifcor to není možné? jen podle tune.args+species a nějak detekcí všech verzí z version2 ???!!!
#

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
    mutate(nRep = n_distinct(rep)) %>%
    filter(nRep >= minReps) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE)
  # vyloučené sloupce doplním znovu zpět
  summarised.not <- setdiff(names(tbl), names(tbl.avg)) # odstranit sloupec rep - nedává smysl
  tbl.not <- tbl %>%
    group_by(version, adjust, tune.args, species) %>%
    mutate(nRep = n_distinct(rep)) %>%
    filter(nRep >= minReps) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::select(all_of(summarised.not))

  tbl.avg %<>% add_column(tbl.not)

  tbl <- tbl.avg
  tbl %<>% ungroup() %>% mutate(id = row_number())
  saveRDS(tbl, paste0(path.eval, "tbl.rds"))
}

#
### version2 - rozdělení zgroupovaného sloupce vifcor-em na samostatné řádky
#

tbl.old <- tbl # záloha
tbl.longer <- tbl %>%
  ungroup() %>%
  mutate(id.old = id) %>%
  separate_longer_delim(version2, delim = "|") %>%
  separate_wider_delim(version2, delim = "_", names = c("versionD", "adjustD"), too_few = "align_start", cols_remove = FALSE) %>%
  rename(version.old = version) %>%
  rename(version = versionD) %>%
  mutate(id = row_number()) # raději přeindexuju

tbl <- tbl.longer


tbl.rep.old <- tbl.rep # záloha
tbl.rep.longer <- tbl.rep %>%
  ungroup() %>%
  separate_longer_delim(version2, delim = "|") %>%
  separate_wider_delim(version2, delim = "_", names = c("versionD", "adjustD"), too_few = "align_start", cols_remove = FALSE) %>%
  rename(version.old = version) %>%
  rename(version = versionD)

tbl.rep <- tbl.rep.longer

###

# vu <- unique(tbl$version)

# tbl %<>% add_row(tbl_5)

# saveRDS(tbl, paste0(path.eval, "tbl_5_6.rds"))
# tbl <- readRDS(modelsResults.avg)

# selection <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "un")
# # selection.f <- c("ssos.topA100", "ssos.topS10", "tgob", "topA100", "topS10", "un", "ssos") # soos záměrně na konci
# selection.f2 <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "un")
# selection.rename <- c("ssosTSAO", "ssosTGOB", "ssosTO", "ssosTS", "TGOB", "TO", "TS", "TSAO", "thin")


# selection <- c(
#   "TGOB",
#   "TS.w",
#   "TO.w",
#   "TGOB.sso.w",
#   "un"
# )
# selection <- c(
#   "TGOB",
#   "TGOB.bss",
#   "TO.w",
#   "TS.w",
#   "TGOB.sso.w",
#   "TGOB.sso.w.3",
#   "TGOB.sso.w.dnorm1.median",
#   "TGOB.sso.w.3.dnorm1.median",
#   "TGOB.sso.w.dnorm1.mean",
#   "TGOB.sso.w.3.dnorm1.mean"  ,
#   "TGOB.sso.w.dnorm2.median",
#   "TGOB.sso.w.3.dnorm2.median",
#   "TGOB.sso.w.dnorm2.mean",
#   "TGOB.sso.w.3.dnorm2.mean",
#   "un"
# )

selection <- c(
  "TGOB",
  "TGOB.sso",
  "TGOB.sso.w",
  # "TGOB.sso.w.3p",
  "TS.w",
  "TS.w.sso",

  "TO.w",
  "TO.w.sso",

  "TGOB.bss",

  "un"
)




# selection.f <- c("ssos.topA100", "ssos.topS10", "tgob", "topA100", "topS10", "un", "ssos") # soos záměrně na konci
selection.f2 <- selection
# selection.rename <- c(
#   "TGOB", "ssoTGOB",
#   "TSAO", "ssoTSAO",
#   "TO", "ssoTO",
#   "TS", "ssoTS",
#   "thin"
# )

# selection.rename <- c(
#   "TGOB",
#   "TGOB.sso",
#     "TGOB.sso.w2",
#   #"TGOB.sso.w.3p",
#   #"TGOB.bss",

#   ".thin" # tečka dočasně kvůli řazení v legendě
# )

selection.rename <- str_replace(selection, "un", ".thin")



selection.f2 <- paste0("^", selection.f2, "$")


names(selection.rename) <- selection.f2

pairs.compare <- list()
pairs.compare[[selection.rename[1]]] <- selection.rename[c(1, 2)]

# pairs.compare[[selection.rename[2]]] <- selection.rename[c(1, 3)]
# pairs.compare[[selection.rename[3]]] <- selection.rename[c(2, 3)]


selection.colors <- c(
  "#87CEFA",
  "#1E90FF",
  "#0066cc",

  "#7CFC00",
  "#008000",

  "#FA8072",
  "#FF0000",

  "#000000",

  "#C0C0C0"
)

# #
# # barvy generuju... - pokud si nedám vlastní výše
# #
# vu <- sort(unique(tbl$version))
# vu.colors <- palette(rainbow(length(vu)))
# selection.colors <- vu.colors
# selection <- selection.f2 <- selection.rename <- vu

# # znovu z dřívějška
# selection.f2 <- paste0("^", selection.f2, "$")
# names(selection.rename) <- selection.f2
# pairs.compare <- list()
# pairs.compare[[selection.rename[1]]] <- selection.rename[c(1, 2)]
# #
# # barvy generuju... (konec)
# #


# selection.colors <- paste0("#", selection.colors)
names(selection.colors) <- paste0("^", selection.rename, "$")

#
# dřívější verze s replikacemi pro t.test-y (rozdíl dvou nej verzí), neřeším null?
#

tbl.rep %<>% ungroup() %>% mutate(id = row_number())

# vyhodit kombinace s malým počtem replikací
tbl.rep.minReps <- tbl.rep %>%
  group_by(version, adjust, tune.args, species) %>%
  mutate(nRep = n_distinct(rep)) %>%
  filter(nRep < minReps) %>%
  ungroup() %>%
  dplyr::select(id)

tbl.rep.minReps <- unique(unlist(tbl.rep.minReps))
tbl.rep %<>% filter(id %notin% tbl.rep.minReps)

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

tbl.null.test.cor <- tbl.null.ids %>%
  group_by(species) %>%
  slice_max(cor, with_ties = FALSE) %>%
  dplyr::select(cor, species, id) %>%
  rename(cor_null = cor) %>%
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

  tbl %<>% left_join(tbl.null.test.cor, by = c("species"), suffix = c("", "__cor")) %>%
    mutate(cordiff = ifelse(cor > cor_null, cor - cor_null, 0)) %>%
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

  tbl %<>% left_join(tbl.null.test.cor, by = c("species"), suffix = c("", "__cor")) %>%
    mutate(cordiff = cor - cor_null) %>%
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

path.PP.val <- paste0(path.PP, "val/")
dir.create(path.PP.val, showWarnings = FALSE)
path.PP.val.test <- paste0(path.PP, "val.test/")
dir.create(path.PP.val.test, showWarnings = FALSE)
path.PP.test <- paste0(path.PP, "test/")
dir.create(path.PP.test, showWarnings = FALSE)

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
combs.n <- 3
if (anyAlt) {
  tbl.f %<>% mutate(version = ifelse(version == "TGOB", "TGOB", "anyAlt")) %>%
    mutate(clr = ifelse(clr == selection.colors[1], selection.colors[1], "#FF00FF")) %>%
    ungroup()
  selection.f2 <- c(selection.f2[1], "#FF00FF")
  combs.n <- 2
}

k6 <- comb_all(selection.f2, combs.n) # length(selection.f2)

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
  # níže zatím řeším jen korelovanost jak nejlepší val model odpovídá témuž test - neřeším ale nakolik nejlepší val trefuje nejlepší test - to ale dělám ve "val.test", které jsem zatím stopnul uprostřed změn... - dodělat!!!???










  #
  # @@@ nelze použít v současném systému verzí, nevím která je první vybraná ve version, musel bych nějak parsovat version2
  #


  # zasobnik <- list()
  # for (ver in unique(tbl.nn$version)) {
  #   tmp <- tbl.nn %>%
  #     ungroup() %>%
  #     filter(auc.val.avg >= at) %>%
  #     filter(version == ver)
  #   ct <- cor.test(tmp$auc.val.avg, tmp$AUC)
  #   zasobnik[[ver]][["cor"]] <- ct[["estimate"]][["cor"]]
  #   zasobnik[[ver]][["p.value"]] <- ct[["p.value"]]
  #   zasobnik[[ver]][["ci.min"]] <- ct[["conf.int"]][1]
  #   zasobnik[[ver]][["ci.max"]] <- ct[["conf.int"]][2]
  # }
  # zasobnik.t <- t(as_tibble(zasobnik))
  # rn <- row.names(zasobnik.t)
  # cn <- names(zasobnik[[1]])
  # zasobnik.t <- as.data.frame(zasobnik.t)
  # zasobnik.t <- as_tibble(sapply(zasobnik.t, as.numeric))
  # names(zasobnik.t) <- cn
  # zasobnik.t$version <- rn
  # zasobnik0[[as.character(at)]][["all"]] <- zasobnik.t %<>% arrange(desc(cor))
  # write.csv(zasobnik.t, paste0(path.PP, "korelace-", as.character(at), ".csv"), row.names = FALSE)
  #
  # # max verze
  # zasobnik <- list()
  # for (ver in unique(tbl.nn$version)) {
  #   tmp <- tbl.nn %>%
  #     ungroup() %>%
  #     group_by(species, version) %>%
  #     slice_max(auc.val.avg, with_ties = FALSE) %>%
  #     filter(auc.val.avg >= at) %>%
  #     filter(version == ver)
  #   ct <- cor.test(tmp$auc.val.avg, tmp$AUC)
  #   zasobnik[[ver]][["cor"]] <- ct[["estimate"]][["cor"]]
  #   zasobnik[[ver]][["p.value"]] <- ct[["p.value"]]
  #   zasobnik[[ver]][["ci.min"]] <- ct[["conf.int"]][1]
  #   zasobnik[[ver]][["ci.max"]] <- ct[["conf.int"]][2]
  # }
  # zasobnik.t <- t(as_tibble(zasobnik))
  # rn <- row.names(zasobnik.t)
  # cn <- names(zasobnik[[1]])
  # zasobnik.t <- as.data.frame(zasobnik.t)
  # zasobnik.t <- as_tibble(sapply(zasobnik.t, as.numeric))
  # names(zasobnik.t) <- cn
  # zasobnik.t$version <- rn
  # zasobnik0[[as.character(at)]][["max"]] <- zasobnik.t %<>% arrange(desc(cor))
  # write.csv(zasobnik.t, paste0(path.PP, "korelace-max-", as.character(at), ".csv"), row.names = FALSE)

  #
  # @@@ nelze použít v současném systému verzí, nevím která je první vybraná ve version, musel bych nějak parsovat version2
  #









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

  gc()
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

      # 111
      # top2v <- lapply(1:2, function(x) {
      #   tmp <- list()
      #   tmp[["t"]] <- tbl.rep.nn %>% filter(species == sp & version == top2[x, ]$version & adjust == top2[x, ]$adjust & tune.args == top2[x, ]$tune.args)
      #   tmp[["top"]] <- top2[x, ]$version
      #
      #   # p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution, we can assume the normality.
      #   # kontrolu na alespoň 3 replikace musím udělat úplně na začátku!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #   if (nrow(tmp[["t"]]) >= 333) {
      #     tmp[["shapiro.p"]] <- shapiro.test(tmp[["t"]]$auc.val.avg)$p.value
      #   } else {
      #     print(" < 3 replikací !!!")
      #     tmp[["shapiro.p"]] <- 0
      #   }
      #   return(tmp)
      # })

      top2v.normality <- FALSE # sum(sapply(top2v, function(x) x$shapiro.p > 0.05)) == 2

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

      if (nrow(tmp.null.t) >= 333) {
        tmp.null <- shapiro.test(tmp.null.t$auc.val.avg)$p.value
      } else {
        print("null < 3 replikací !!!")
        tmp.null <- 0
      }

      top1null.normality <- FALSE # sum(c(tmp.null > 0.05, top2v[[1]]$shapiro.p > 0.05)) == 2

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

  temp.g.next <- temp.g
  temp.g.clrn <- temp.g
  if (anyAlt) {
    ## stejný výsledek
    temp.g.next <- temp.g
    clrn <- temp.g %>%
      ungroup() %>%
      group_by(species) %>%
      slice_max(auc.val.avgdiff, with_ties = TRUE) %>%
      group_by(species) %>%
      summarise(clrn = length(id))
    clrn %<>% filter(clrn > 1)
    temp.g %<>% mutate(clr = ifelse(species %in% clrn$species, "#808080", clr))
    temp.g %<>% mutate(version = ifelse(species %in% clrn$species, "2+", version))
    temp.g.clrn <- temp.g
  }

  tobs <- ggplot() +
    geom_bar(data = temp.g %>% ungroup() %>% group_by(species) %>% slice_max(auc.val.avgdiff, with_ties = FALSE), stat = "identity", mapping = aes(species, auc.val.avgdiff, fill = factor(version)), color = "yellow", size = 0.01) +
    geom_point(data = temp.g %>% ungroup(), mapping = aes(species, auc.val.avgdiff, fill = factor(version)), size = 1, stroke = 0.1, color = "yellow", shape = 21, alpha = 0.90) +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      legend.text = element_text(size = 4),
      text = element_text(size = 6),
      axis.text.y = element_markdown(),
      # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )

  tt.null <- sapply(tmp.top1null[["all"]], function(x) x$val < 0.05)
  tt <- sapply(tmp.top2[["all"]], function(x) x$val < 0.05)

  # negativní korekce
  tt.null.minus <- sapply(tmp.top1null[["all"]], function(x) is.na(x$val.auc1 < x$val.auc2) | x$val.auc1 < x$val.auc2)
  # znegovat případné minusové opravy
  tt.null[tt.null == TRUE & tt.null.minus == TRUE] <- FALSE
  tt[tt == TRUE & tt.null.minus == TRUE] <- FALSE

  # omezením druhy pro AUC treshold
  tt.null <- tt.null[names(tt.null) %in% unique(temp.g$species)]
  tt <- tt[names(tt) %in% unique(temp.g$species)]

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
  ggsave(paste0(path.PP.val, "trend-overall-best-species.val.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")


  # occs.n
  no.temp <- no %>%
    arrange(occs.n)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: sum of occupied squares); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP.val, "trend-overall-best-species-orderByOccs.val.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  # auc.val.avgdiff
  no.temp <- no %>%
    arrange(auc.val.avgdiff)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avgdiff); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP.val, "trend-overall-best-species-orderByDiff.val.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  # auc.val.avg
  no.temp <- no %>%
    arrange(auc.val.avg)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avg); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP.val, "trend-overall-best-species-orderByAuc.val.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  # auc.val.avg_null
  no.temp <- no %>%
    arrange(auc.val.avg_null)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avg_null); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP.val, "trend-overall-best-species-orderByAucNull.val.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  ### ### ### endB

  #
  # páry
  #
  for (cmp in pairs.compare) {
    if (!anyAlt) {
      temp.g.c <- temp.g


      temp.g.c %<>% filter(version %in% cmp)

      gc()
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

          # 111
          #         top2v <- lapply(1:2, function(x) {
          #           tmp <- list()
          #           tmp[["t"]] <- tbl.rep.nn %>% filter(species == sp & version == top2[x, ]$version & adjust == top2[x, ]$adjust & tune.args == top2[x, ]$tune.args)
          #           tmp[["top"]] <- top2[x, ]$version
          #
          #           # p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution, we can assume the normality.
          #           # kontrolu na alespoň 3 replikace musím udělat úplně na začátku!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          #           if (nrow(tmp[["t"]]) >= 333) {
          #             tmp[["shapiro.p"]] <- shapiro.test(tmp[["t"]]$auc.val.avg)$p.value
          #           } else {
          #             print(" < 3 replikací !!!")
          #             tmp[["shapiro.p"]] <- 0
          #           }
          #           return(tmp)
          #         })

          top2v.normality <- FALSE # sum(sapply(top2v, function(x) x$shapiro.p > 0.05)) == 2

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

          if (nrow(tmp.null.t) >= 333) {
            tmp.null <- shapiro.test(tmp.null.t$auc.val.avg)$p.value
          } else {
            print("null < 3 replikací !!!")
            tmp.null <- 0
          }

          top1null.normality <- FALSE # sum(c(tmp.null > 0.05, top2v[[1]]$shapiro.p > 0.05)) == 2

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
          legend.text = element_text(size = 4),
          text = element_text(size = 6),
          axis.text.y = element_markdown(),
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
          panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
          strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
        )

      tt.null <- sapply(tmp.top1null[[cmp[1]]], function(x) x$val < 0.05)
      tt <- sapply(tmp.top2[[cmp[1]]], function(x) x$val < 0.05)

      # negativní korekce
      tt.null.minus <- sapply(tmp.top1null[[cmp[1]]], function(x) is.na(x$val.auc1 < x$val.auc2) | x$val.auc1 < x$val.auc2)
      # znegovat případné minusové opravy
      tt.null[tt.null == TRUE & tt.null.minus == TRUE] <- FALSE
      tt[tt == TRUE & tt.null.minus == TRUE] <- FALSE

      # omezením druhy pro AUC treshold
      tt.null <- tt.null[names(tt.null) %in% unique(temp.g.c$species)]
      tt <- tt[names(tt) %in% unique(temp.g.c$species)]


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
      ggsave(paste0(path.PP.val, "trend-overall-best-species.val.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")


      # occs.n
      no.temp <- no %>%
        arrange(occs.n)
      order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
      title <- unname(unlist(no.temp$title))
      title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

      tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: sum of occupied squares); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
      ggsave(paste0(path.PP.val, "trend-overall-best-species-orderByOccs.val.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      # auc.val.avgdiff
      no.temp <- no %>%
        arrange(auc.val.avgdiff)
      order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
      title <- unname(unlist(no.temp$title))
      title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

      tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avgdiff); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
      ggsave(paste0(path.PP.val, "trend-overall-best-species-orderByDiff.val.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      # auc.val.avg
      no.temp <- no %>%
        arrange(auc.val.avg)
      order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
      title <- unname(unlist(no.temp$title))
      title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

      tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avg); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
      ggsave(paste0(path.PP.val, "trend-overall-best-species-orderByAuc.val.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      # auc.val.avg_null
      no.temp <- no %>%
        arrange(auc.val.avg_null)
      order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
      title <- unname(unlist(no.temp$title))
      title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

      tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best auc.val.avg_null); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
      ggsave(paste0(path.PP.val, "trend-overall-best-species-orderByAucNull.val.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      ### ### ### endB
    }
  }


  temp.g <- temp.g.next


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(occs.n, auc.val.avgdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.1) +
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
  ggsave(paste0(path.PP.val, "trend.val.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")



  ### ### ### start
  temp.g <- temp.g.clrn


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
  ggsave(paste0(path.PP.val, "trend-overall.val.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  # a) vše - OVERALL - best
  temp.g.t <- temp.g %>%
    ungroup() %>%
    group_by(species) %>%
    slice_max(auc.val.avgdiff, with_ties = FALSE)
  temp.g.t.clr <- unique(unname(unlist(temp.g.t$clr)))
  temp.clr <- unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))
  temp.clr <- temp.clr[temp.clr %in% temp.g.t.clr] # ne vždy musí být ve výsledku všechny verze, nutno případně dofiltrovat, aby seděly barvy a pořadí

  ggplot(temp.g.t, aes(occs.n, auc.val.avgdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.5) +
    geom_smooth(method = loess, size = 0.2) +
    scale_color_manual(values = temp.clr) +
    theme_light() +
    theme(
      text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )
  ggsave(paste0(path.PP.val, "trend-overall-best.val.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  ### ### ### end

  temp.g <- temp.g.next







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
  ggsave(paste0(path.PP.val, "version-species.val.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")




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
  ggsave(paste0(path.PP.val, "boxplot.val.", as.character(at), ".png"), width = 1500, height = 1000, units = "px")



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
  write.csv(res, paste0(path.PP.val, "combs-val-", as.character(at), ".csv"), row.names = FALSE)


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
  ggsave(paste0(path.PP.val.test, "boxplot.val.test.", as.character(at), ".png"), width = 1500, height = 1000, units = "px")


  # a) vše
  ggplot(temp.g %>% ungroup(), aes(occs.n, auc.val.avgdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.01) +
    geom_smooth(method = loess, size = 0.2) +
    theme_light() +
    theme(
      legend.position = "none", text = element_text(size = 4),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    ) +
    facet_wrap(~version)
  ggsave(paste0(path.PP.val.test, "trend.val.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

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
  ggsave(paste0(path.PP.val.test, "version-species.val.test.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

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
  write.csv(res, paste0(path.PP.val.test, "combs-val.test-", as.character(at), ".csv"), row.names = FALSE)




  ################################################################################################################################################################################
  # kss AUCdiff
  ################################################################################################################################################################################

  #
  # se.plot
  #

  temp.g.tm <- list()
  temp.g.tm[["AUC"]] <- tbl.nn %>%
    ungroup() %>%
    group_by(species, version) %>%
    slice_max(AUCdiff, with_ties = FALSE) # %>% filter(AUC >= at) # nedělám, neznám odpovídající treshold pro cor...

  temp.g.tm[["cor"]] <- tbl.nn %>%
    ungroup() %>%
    group_by(species, version) %>%
    slice_max(cordiff, with_ties = FALSE)

  temp.g.tm.join <- list()

  temp.g.tm.join[["AUC"]] <- temp.g.tm[["AUC"]] %>%
    ungroup() %>%
    group_by(version) %>%
    summarise(aucMean = mean(AUC), aucSd = se(AUC), clr = first(clr))
  temp.g.tm.join[["cor"]] <- temp.g.tm[["cor"]] %>%
    ungroup() %>%
    group_by(version) %>%
    summarise(corMean = mean(cor), corSd = se(cor), clr = first(clr))

  temp.g.tm.joined <- temp.g.tm.join[["AUC"]] %>% left_join(temp.g.tm.join[["cor"]], by = "version")


  se.plot <- ggplot(temp.g.tm.joined, aes(x = aucMean, y = corMean, color = version)) +
    geom_point(shape = 16, alpha = 0.90) +
    scale_color_manual(values = temp.g.tm.joined$clr.x) +
    theme_light() +
    # geom_smooth(method=lm)+
    geom_errorbarh(aes(
      xmin = aucMean - aucSd,
      xmax = aucMean + aucSd
    ),
    height = 0.0, size = 0.2, alpha = 0.90
    ) +
    geom_errorbar(aes(
      ymin = corMean - corSd,
      ymax = corMean + corSd
    ),
    width = 0.0, size = 0.2, alpha = 0.90
    ) +
    theme(
      legend.text = element_text(size = 4),
      text = element_text(size = 4)
    )

  ggsave(paste0(path.PP.test, "se.plot.test.", as.character(at), ".png"), width = 1500, height = 800, units = "px")

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

  gc()
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

      # 111
      # top2v <- lapply(1:2, function(x) {
      #   tmp <- list()
      #   tmp[["t"]] <- tbl.rep.nn %>% filter(species == sp & version == top2[x, ]$version & adjust == top2[x, ]$adjust & tune.args == top2[x, ]$tune.args)
      #   tmp[["top"]] <- top2[x, ]$version
      #
      #   # p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution, we can assume the normality.
      #   # kontrolu na alespoň 3 replikace musím udělat úplně na začátku!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #   if (nrow(tmp[["t"]]) >= 333) {
      #     tmp[["shapiro.p"]] <- shapiro.test(tmp[["t"]]$AUC)$p.value
      #   } else {
      #     print(" < 3 replikací !!!")
      #     tmp[["shapiro.p"]] <- 0
      #   }
      #   return(tmp)
      # })

      top2v.normality <- FALSE # sum(sapply(top2v, function(x) x$shapiro.p > 0.05)) == 2

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

      if (nrow(tmp.null.t) >= 333) {
        tmp.null <- shapiro.test(tmp.null.t$AUC)$p.value
      } else {
        print("null < 3 replikací !!!")
        tmp.null <- 0
      }

      top1null.normality <- FALSE # sum(c(tmp.null > 0.05, top2v[[1]]$shapiro.p > 0.05)) == 2

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

  temp.g.next <- temp.g
  temp.g.clrn <- temp.g
  if (anyAlt) {
    temp.g.next <- temp.g
    clrn <- temp.g %>%
      ungroup() %>%
      group_by(species) %>%
      slice_max(AUCdiff, with_ties = TRUE) %>%
      group_by(species) %>%
      summarise(clrn = length(id))
    clrn %<>% filter(clrn > 1)
    temp.g %<>% mutate(clr = ifelse(species %in% clrn$species, "#808080", clr))
    temp.g %<>% mutate(version = ifelse(species %in% clrn$species, "2+", version))
    temp.g.clrn <- temp.g
  }


  tobs <- ggplot() +
    geom_bar(data = temp.g %>% ungroup() %>% group_by(species) %>% slice_max(AUCdiff, with_ties = FALSE), stat = "identity", mapping = aes(species, AUCdiff, fill = factor(version)), color = "yellow", size = 0.01) +
    geom_point(data = temp.g %>% ungroup(), mapping = aes(species, AUCdiff, fill = factor(version)), size = 1, stroke = 0.1, color = "yellow", shape = 21, alpha = 0.90) +
    scale_fill_manual(values = unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))) +
    theme_light() +
    theme(
      legend.text = element_text(size = 4),
      text = element_text(size = 6),
      axis.text.y = element_markdown(),
      # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )

  tt.null <- sapply(tmp.top1null[["all"]], function(x) x$test < 0.05)
  tt <- sapply(tmp.top2[["all"]], function(x) x$test < 0.05)


  # negativní korekce
  tt.null.minus <- sapply(tmp.top1null[["all"]], function(x) is.na(x$test.auc1 < x$test.auc2) | x$test.auc1 < x$test.auc2)
  # znegovat případné minusové opravy
  tt.null[tt.null == TRUE & tt.null.minus == TRUE] <- FALSE
  tt[tt == TRUE & tt.null.minus == TRUE] <- FALSE

  # omezením druhy pro AUC treshold
  tt.null <- tt.null[names(tt.null) %in% unique(temp.g$species)]
  tt <- tt[names(tt) %in% unique(temp.g$species)]



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
  ggsave(paste0(path.PP.test, "trend-overall-best-species.test.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  # occs.n
  no.temp <- no %>%
    arrange(occs.n)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: sum of occupied squares); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP.test, "trend-overall-best-species-orderByOccs.test.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  # AUCdiff
  no.temp <- no %>%
    arrange(AUCdiff)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUCdiff); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP.test, "trend-overall-best-species-orderByDiff.test.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  # AUC
  no.temp <- no %>%
    arrange(AUC)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUC); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP.test, "trend-overall-best-species-orderByAuc.test.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  # AUC_null
  no.temp <- no %>%
    arrange(AUC_null)
  order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
  title <- unname(unlist(no.temp$title))
  title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

  tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUC_null); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
  ggsave(paste0(path.PP.test, "trend-overall-best-species-orderByAucNull.test.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")

  ### ### ### endB

  #
  # páry
  #
  for (cmp in pairs.compare) {
    if (!anyAlt) {
      temp.g.c <- temp.g


      temp.g.c %<>% filter(version %in% cmp)

      gc()
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

          # 111
          # top2v <- lapply(1:2, function(x) {
          #   tmp <- list()
          #   tmp[["t"]] <- tbl.rep.nn %>% filter(species == sp & version == top2[x, ]$version & adjust == top2[x, ]$adjust & tune.args == top2[x, ]$tune.args)
          #   tmp[["top"]] <- top2[x, ]$version
          #
          #   # p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution, we can assume the normality.
          #   # kontrolu na alespoň 3 replikace musím udělat úplně na začátku!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          #   if (nrow(tmp[["t"]]) >= 333) {
          #     tmp[["shapiro.p"]] <- shapiro.test(tmp[["t"]]$AUC)$p.value
          #   } else {
          #     print(" < 3 replikací !!!")
          #     tmp[["shapiro.p"]] <- 0
          #   }
          #   return(tmp)
          # })

          top2v.normality <- FALSE # sum(sapply(top2v, function(x) x$shapiro.p > 0.05)) == 2

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

          if (nrow(tmp.null.t) >= 333) {
            tmp.null <- shapiro.test(tmp.null.t$AUC)$p.value
          } else {
            print("null < 3 replikací !!!")
            tmp.null <- 0
          }

          top1null.normality <- FALSE # sum(c(tmp.null > 0.05, top2v[[1]]$shapiro.p > 0.05)) == 2

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
          legend.text = element_text(size = 4),
          text = element_text(size = 6),
          axis.text.y = element_markdown(),
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
          panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
          strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
        )

      tt.null <- sapply(tmp.top1null[[cmp[1]]], function(x) x$test < 0.05)
      tt <- sapply(tmp.top2[[cmp[1]]], function(x) x$test < 0.05)

      # negativní korekce
      tt.null.minus <- sapply(tmp.top1null[[cmp[1]]], function(x) is.na(x$test.auc1 < x$test.auc2) | x$test.auc1 < x$test.auc2)
      # znegovat případné minusové opravy
      tt.null[tt.null == TRUE & tt.null.minus == TRUE] <- FALSE
      tt[tt == TRUE & tt.null.minus == TRUE] <- FALSE

      # omezením druhy pro AUC treshold
      tt.null <- tt.null[names(tt.null) %in% unique(temp.g.c$species)]
      tt <- tt[names(tt) %in% unique(temp.g.c$species)]

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
      ggsave(paste0(path.PP.test, "trend-overall-best-species.test.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      # occs.n
      no.temp <- no %>%
        arrange(occs.n)
      order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
      title <- unname(unlist(no.temp$title))
      title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

      tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: sum of occupied squares); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
      ggsave(paste0(path.PP.test, "trend-overall-best-species-orderByOccs.test.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      # AUCdiff
      no.temp <- no %>%
        arrange(AUCdiff)
      order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
      title <- unname(unlist(no.temp$title))
      title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

      tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUCdiff); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
      ggsave(paste0(path.PP.test, "trend-overall-best-species-orderByDiff.test.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      # AUC
      no.temp <- no %>%
        arrange(AUC)
      order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
      title <- unname(unlist(no.temp$title))
      title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

      tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUC); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
      ggsave(paste0(path.PP.test, "trend-overall-best-species-orderByAuc.test.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      # AUC_null
      no.temp <- no %>%
        arrange(AUC_null)
      order.new <- unname(unlist(no.temp %>% dplyr::select(species)))
      title <- unname(unlist(no.temp$title))
      title.occs.n <- paste0(title, " | ", sprintf("%04d", unname(unlist(no.temp$occs.n))))

      tobs + scale_x_discrete(labels = title.occs.n, limits = order.new) + xlab("species (ordered by: best AUC_null); AUCnull->AUCbest | sum of occupied squares") + coord_flip()
      ggsave(paste0(path.PP.test, "trend-overall-best-species-orderByAucNull.test.", as.character(at), "---", cmp[1], ".png"), width = 2000, height = 2000, units = "px")

      ### ### ### endB
    }
  }

  temp.g <- temp.g.next
  # a) vše
  ggplot(temp.g %>% ungroup(), aes(occs.n, AUCdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.1) +
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
  ggsave(paste0(path.PP.test, "trend.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")


  ### ### ### start

  temp.g <- temp.g.clrn
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
  ggsave(paste0(path.PP.test, "trend-overall.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")

  # a) vše - OVERALL - best
  temp.g.t <- temp.g %>%
    ungroup() %>%
    group_by(species) %>%
    slice_max(AUCdiff, with_ties = FALSE)
  temp.g.t.clr <- unique(unname(unlist(temp.g.t$clr)))
  temp.clr <- unique(unname(unlist(temp.g %>% group_by(version) %>% slice_head(n = 1) %>% ungroup() %>% arrange(tolower(version)) %>% dplyr::select(clr))))
  temp.clr <- temp.clr[temp.clr %in% temp.g.t.clr] # ne vždy musí být ve výsledku všechny verze, nutno případně dofiltrovat, aby seděly barvy a pořadí

  ggplot(temp.g.t, aes(occs.n, AUCdiff)) +
    geom_point(aes(colour = factor(version)), size = 0.5) +
    geom_smooth(method = loess, size = 0.2) +
    scale_color_manual(values = temp.clr) +
    theme_light() +
    theme(
      text = element_text(size = 6),
      panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.1),
      strip.text = element_text(margin = margin(t = 1, r = 1, b = 1, l = 1))
    )
  ggsave(paste0(path.PP.test, "trend-overall-best.test.", as.character(at), ".png"), width = 2000, height = 1500, units = "px")



  ### ### ### end


  temp.g <- temp.g.next



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
  ggsave(paste0(path.PP.test, "version-species.test.", as.character(at), ".png"), width = 2000, height = 2000, units = "px")





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
  ggsave(paste0(path.PP.test, "boxplot.test.", as.character(at), ".png"), width = 1500, height = 1000, units = "px")




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
  write.csv(res, paste0(path.PP.test, "combs-test-", as.character(at), ".csv"), row.names = FALSE)
}

if (verified.generate) {
  saveRDS(tmp.top2, paste0(path.eval, verified.top2))
  saveRDS(tmp.top1null, paste0(path.eval, verified.top1null))
}
gc()

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
# significant improvement (t.test)
#

tmp.top1null.out <- list()
tmp.top2.out <- list()
tmp.out <- list()
tmp.summary <- list()
# převedení verified na tabulky: top1null (nej vs null)
tmp.top1null.out <- list()
for (tn in names(tmp.top1null)) {
  #
  ### # top1null
  #
  tmp.top1null.t <- t(as_tibble(tmp.top1null[[tn]]))
  rn <- row.names(tmp.top1null.t)
  cn <- names(tmp.top1null[[tn]][[1]])
  tmp.top1null.t <- as_tibble(tmp.top1null.t)

  tmp.top1null.t <- as_tibble(sapply(tmp.top1null.t, as.character))
  names(tmp.top1null.t) <- cn
  tmp.top1null.t$species <- rn
  tmp.top1null.t <- as_tibble(tmp.top1null.t)

  # změna datových typů sloupců podle obsahu
  tmp.top1null.t %<>% mutate_all(funs(type.convert(as.character(.))))

  #
  ### # top2
  #
  tmp.top2.t <- t(as_tibble(tmp.top2[[tn]]))
  rn <- row.names(tmp.top2.t)
  cn <- names(tmp.top2[[tn]][[1]])
  tmp.top2.t <- as_tibble(tmp.top2.t)

  tmp.top2.t <- as_tibble(sapply(tmp.top2.t, as.character))
  names(tmp.top2.t) <- cn
  tmp.top2.t$species <- rn
  tmp.top2.t <- as_tibble(tmp.top2.t)

  # změna datových typů sloupců podle obsahu
  tmp.top2.t %<>% mutate_all(funs(type.convert(as.character(.))))

  # průběžné výsledky
  tmp.top1null.out[[tn]] <- tmp.top1null.t
  write.csv(tmp.top1null.t, paste0(path.PP, "verified-top1null-part-.", tn, ".csv"), row.names = FALSE)
  tmp.top2.out[[tn]] <- tmp.top2.t
  write.csv(tmp.top2.t, paste0(path.PP, "verified-top2-part-.", tn, ".csv"), row.names = FALSE)

  # výběr jen významných rozdílů + lepších než null
  tmp.joined <- tmp.top2.t %>% left_join(tmp.top1null.t, by = "species", suffix = c("", "_null"))
  tmp.joined %<>% mutate(val.lessThanNull = ifelse(val.auc1 < val.auc2_null, 1, 0)) %<>% mutate(test.lessThanNull = ifelse(test.auc1 < test.auc2_null, 1, 0))
  tmp.joined %<>% mutate(val.different = ifelse(val.lessThanNull == 0 & val < 0.05, 1, 0)) %<>% mutate(test.different = ifelse(test.lessThanNull == 0 & test < 0.05, 1, 0)) # !!! správně val_null test_null - ale je tam něco jiného!?!?!!!

  tmp.out[[tn]] <- tmp.joined
  write.csv(tmp.joined, paste0(path.PP, "verified-joined-part-.", tn, ".csv"), row.names = FALSE)

  # souhrny úspěšných
  tmp.summary[[tn]][["val"]] <- sum(unlist(tmp.joined$val.different) == 0, na.rm = TRUE)
  tmp.summary[[tn]][["test"]] <- sum(unlist(tmp.joined$test.different) == 0, na.rm = TRUE)
}

saveRDS(tmp.top1null.out, paste0(path.PP, "top1null.out.rds"))
saveRDS(tmp.top2.out, paste0(path.PP, "top2.out.rds"))
saveRDS(tmp.out, paste0(path.PP, "joined.out.rds"))
saveRDS(tmp.summary, paste0(path.PP, "summary.out.rds"))


#
# souhrny
#

pairs.main <- names(pairs.compare)
pairs.sso <- sapply(pairs.compare, function(x) paste(x, collapse = "|"))
combs.select.res.both <- list()

for (v in c("val", "test")) {
  print(v)
  # párový přínost sso verzí
  combs.select <- c(pairs.main, pairs.sso)
  combs.select.res <- combs[["0"]][[v]] %>%
    filter(versionComb %in% combs.select) %>%
    arrange(versionComb) %>%
    dplyr::select(AUCdiffSum, mean, AUCdiffMedian, versionComb)
  # manual reorder
  combs.select.res <- combs.select.res[c(1:5, 8, 6, 7), ]
  write.csv(combs.select.res, paste0(path.PP, "pair-sso-version-cumsum.", v, ".csv"), row.names = FALSE)

  # TGOB + další řádné + sso
  tgob.main.res <- combs[["0"]][[v]] %>%
    filter(versionComb %in% c(pairs.main[1], paste(pairs.main, collapse = "|"))) %>%
    arrange(versionComb) %>%
    dplyr::select(AUCdiffSum, mean, AUCdiffMedian, versionComb)
  tgob.main.res %<>% add_row(combs[["0"]][[v]] %>% slice_head(n = 1) %>% dplyr::select(AUCdiffSum, mean, AUCdiffMedian, versionComb))
  write.csv(tgob.main.res, paste0(path.PP, "main-sso-version-cumsum.", v, ".csv"), row.names = FALSE)

  #
  ### souhrnné tabulky
  #
  combs.select.res.both[[v]][["t1"]] <- combs.select.res %>%
    mutate(AUCdiffSum_contrib = ifelse(row_number() %% 2 == 0, AUCdiffSum - lag(AUCdiffSum), NA)) %>%
    mutate(mean_contrib = ifelse(row_number() %% 2 == 0, mean - lag(mean), NA)) %>%
    mutate(AUCdiffMedian_contrib = ifelse(row_number() %% 2 == 0, AUCdiffMedian - lag(AUCdiffMedian), NA))

  t1a <- combs.select.res.both[[v]][["t1"]] %>%
    dplyr::select(c(1, 2, 4)) %>%
    rename(version = "versionComb")
  t1b <- combs.select.res.both[[v]][["t1"]] %>%
    dplyr::select(4:6) %>%
    na.omit()
  t1a <- t1a[c(1, 3, 5, 7), ]
  t1a$s <- "|"
  tAll <- t1a %>% add_column(t1b)
  t1 <- tAll[, c(3, 1, 2, 4:7)]

  combs.select.res.both[[v]][["t2"]] <- tgob.main.res %>%
    mutate(AUCdiffSum_contrib = AUCdiffSum - lag(AUCdiffSum)) %>%
    mutate(mean_contrib = mean - lag(mean)) %>%
    mutate(AUCdiffMedian_contrib = AUCdiffMedian - lag(AUCdiffMedian))

  # statisticky významné zlepšení druhů
  tmp.different.test <- sapply(tmp.summary, function(x) x[[v]])

  combs.select.res.both[[v]][["t3"]] <- as_tibble(list("version" = names(tmp.different.test), speciesImproved = tmp.different.test)) %>% arrange(tolower(version))
  combs.select.res.both[[v]][["t3"]]$version <- c("all", t1b$versionComb) # přepsat poárovými, tak je to ve skutečnosti...
}




#
# reporty
#

# # včetně negativních!? odělat i odfiltrovanou verzi bez negativních hodnot diffů!!! - stačilo by nebrat negativní hodnoty?
# t1.temp <- combs.select.res.both[["test"]] %>%
#   mutate(AUCdiffSum_contrib = ifelse(row_number() %% 2 == 0, AUCdiffSum - lag(AUCdiffSum), NA)) %>%
#   mutate(mean_contrib = ifelse(row_number() %% 2 == 0, mean - lag(mean), NA)) %>%
#   mutate(AUCdiffMedian_contrib = ifelse(row_number() %% 2 == 0, AUCdiffMedian - lag(AUCdiffMedian), NA))

# t1a <- t1.temp %>%
#   dplyr::select(c(1, 2, 4)) %>%
#   rename(version = "versionComb")
# t1b <- t1.temp %>%
#   dplyr::select(4:6) %>%
#   na.omit()
# t1a <- t1a[c(1, 3, 5, 7), ]
# t1a$s <- "|"
# tAll <- t1a %>% add_column(t1b)
# t1 <- tAll[, c(3, 1, 2, 4:7)]

# t2 <- tgob.main.res %>%
#   mutate(AUCdiffSum_contrib = AUCdiffSum - lag(AUCdiffSum)) %>%
#   mutate(mean_contrib = mean - lag(mean)) %>%
#   mutate(AUCdiffMedian_contrib = AUCdiffMedian - lag(AUCdiffMedian))

# # statisticky významné zlepšení druhů
# tmp.different.test <- sapply(tmp.summary, function(x) x$test)
# t3 <- as_tibble(list("version" = names(tmp.different.test), speciesImproved = tmp.different.test)) %>% arrange(tolower(version))
# t3$version <- c("all", t1$versionComb) # přepsat poárovými, tak je to ve skutečnosti...

rmarkdown::render(paste0(path.wd, "report.Rmd"), "all", paste0(path.PP, "report.html"))


end_time <- Sys.time()
print(end_time - start_time)

# tmp.top2 <- list()
# tmp.top1null <- list()

# for (l1 in names(verified.top1nullA)) {
#   for (l2 in names(verified.top1nullA[[l1]])) {
#     tmp.top1null[[l1]][[l2]] <- append(verified.top1nullA[[l1]][[l2]], verified.top1nullB[[l1]][[l2]])
#     tmp.top2[[l1]][[l2]] <- append(verified.top2A[[l1]][[l2]], verified.top2B[[l1]][[l2]])
#   }
# }
# saveRDS(tmp.top2, paste0(path.eval, verified.top2))
# saveRDS(tmp.top1null, paste0(path.eval, verified.top1null))