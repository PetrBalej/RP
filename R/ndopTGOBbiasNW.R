cmd_arg <- commandArgs(trailingOnly = TRUE)
if (is.na(cmd_arg[1])) {
  print("*********************************************************************************************************************")
  print("nezadán parametr pro část druhů")
  print("*********************************************************************************************************************")
  cmd_arg <- "1_2"
} else {
  cmd_arg <- cmd_arg[1]
}

# aa.tbl <- (`t_90___ssos_Accipiter nisus_1` %>% mutate(aa = AUC_trainRandom - auc.train) %>% mutate(ac = (AUC_trainRandom+abs(cor_trainRandom))/2) %>%  dplyr::select(AUC, auc.train, AICc, cor.manual.estimate, cor.manual.p.value, AUC_trainRandom, cor_trainRandom,aa, ac))
# plot(`t_90___ssos_Accipiter nisus_1` %>% mutate(aa = AUC_trainRandom - auc.train) %>% mutate(ac = (AUC_trainRandom+abs(cor_trainRandom))/2) %>%  dplyr::select(AUC, auc.train, AICc, cor.manual.estimate, cor.manual.p.value, AUC_trainRandom, cor_trainRandom,aa, ac))
# plot(`t_90___ssos_Emberiza citrinella_1`  %>% mutate(aa = AUC_trainRandom - auc.train) %>% mutate(ac = (AUC_trainRandom+abs(cor_trainRandom))/2)  %>% filter(str_detect(version, "kernel2"))  %>%  dplyr::select(AUC, auc.train, AICc, cor.manual.estimate, cor.manual.p.value, AUC_trainRandom, cor_trainRandom,aa, ac))
# plot(`t_90___ssos_Emberiza citrinella_1`  %>% mutate(aa = AUC_trainRandom - auc.train) %>% mutate(ac = (AUC_trainRandom+abs(cor_trainRandom))/2)  %>% filter(str_detect(version, "kernel2|un"))  %>%  dplyr::select(AUC, auc.train, AICc, cor.manual.estimate, cor.manual.p.value, AUC_trainRandom, cor_trainRandom,aa, ac))





# cc.all <- c()
# for(cid in paste0("kernel", c(200:205, 210:215, 220:225, 230:235))){
#
#   cc <- cor(`t_90___ssos_Accipiter gentilis_1`  %>% mutate(aa = AUC_trainRandom - auc.train) %>% mutate(ac = (AUC_trainRandom+abs(cor_trainRandom))/2)  %>% filter(version == cid)  %>%  dplyr::select(auc.train, cor.manual.estimate))
#
#   print(cid)
#   print(cc)
#
#
#   cc.all[as.character(cid)] <- cc[1,2]
# }
#
# plot(cc.all)
# axis(1,at=1:length(cc.all), labels=names(cc.all))




print(cmd_arg)


ss <- as.numeric(unlist(strsplit(cmd_arg, split = "_")))


print(ss[1])
print(ss[2])
# Rscript /mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/R/ndopTGOB.R


# "C:\Users\svc.kagup\AppData\Local\Programs\R\R-4.2.1\bin\x64\Rscript.exe" "D:\PERSONAL_DATA\Balej\RP\RP\R\ndopTGOBbiasNW.R"
# "C:\Program Files\R\R-4.0.5\bin\x64\Rscript.exe" "D:\PERSONAL_DATA\Balej\RP\RP\R\ndopTGOBbiasNW.R" 1

# "C:\Program Files\R\R-4.2.1\bin\x64\Rscript.exe" "D:\PersonalWork\Balej\v2\RP\RP\R\ndopTGOBbiasNW.R" 1
# library(devtools)
# install_github("PetrBalej/tgbg", force = TRUE)
# install_github("PetrBalej/tgbg", force = TRUE, ref="dev", build=TRUE)
# install_local("D:/PersonalWork/Balej/tgbg-dev/tgbg-dev", force = TRUE, build=TRUE)
# install_local("/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/tgbg/tgbg", force = TRUE, build=TRUE)
# install_local("/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/maxnet-master", force = TRUE, build=TRUE)


# vlastní nutné modifikace existujících balíčků:
# install_github("PetrBalej/sdm", force = TRUE, ref="advancedTresholdMetrics", build=TRUE)
# install_github("PetrBalej/usdm", force = TRUE, ref="pb", build=TRUE)
# install_github("PetrBalej/ENMeval", force = TRUE, build=TRUE) # přímo v masteru (+ je vhodné ponechávat původní název? - různé plusy i minusy, už jsem řešil...)
# devtools::document("/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/tgbg/tgbg")

# *** Running ENMeval v2.0.4 with maxnet from maxnet package v0.1.4 *** -- na WS1
# "C:/Users/svc.kagup/AppData/Local/Programs/R/R-4.2.1/bin/x64/Rscript.exe" "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/R/ndopTGOB.R" 1
# que %<>% mutate(rep = ifelse(DRUH %in% c("Phoenicurus ochruros","Phoenicurus phoenicurus","Phylloscopus trochilus","Pica pica","Picus canus"), 0, rep))

start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "ENMeval", "ecospat", "blockCV", "tgbg", "sdmATM", "usdmPB", "classInt", "PRROC") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())), repos = "https://mirrors.nic.cz/R/") # určení repos vyžadované na WS1

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)


print(cmd_arg)

# sessionInfo()
# stop()

"%notin%" <- Negate("%in%")

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

# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- paste0(path.wd, "../../projects-data/") # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.prep <- paste0(path.wd, "../../dataPrep/")
source(paste0(path.wd, "shared.R"))
path.rgee <- paste0(path.wd, "../../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.ndop <- paste0(path.prep, "ndop/") # path.wd.prep.ndop
path.lsd <- paste0(path.prep, "lsd/") # path.wd.prep.ndop
path.tgob <- paste0(path.prep, "ndopTGOBimp999q/") # path.wd.prep.ndop
path.tgob.b <- paste0(path.prep, "ndopTGOBimp999b/") # path.wd.prep.ndop
path.tgob.tgbg <- paste0(path.prep, "ndopTGOBtgbg999p/") # path.wd.prep.ndop

# path.tgob.tgbg <- paste0(path.prep, "ndopTGOBtgbg999p-",cmd_arg, "/") # path.wd.prep.ndop
# dir.create(path.tgob.tgbg, showWarnings = FALSE)


############
# inputs
############


predictors <- stack(readRDS(paste0(path.prep, "predictors/20192022/preds_final_vifs.rds"))) # vifstep2 pcamap6
# predictors <- predictors/sum(as.vector(predictors), na.rm = TRUE) # na WS1 udělá RasterBrick
predictors <- stack(predictors / sum(as.vector(predictors), na.rm = TRUE))

# predictors3 <- disaggregate(predictors[[1]], fact=3)



# NDOP dofiltrované pro SSOS
# NDOP druhy presence
ndopP <- readRDS(paste0(path.ndop, "ndopP.rds")) %>% dplyr::select(POLE, DRUH, AUTOR)
ndopP.POLE <- readRDS(paste0(path.ndop, "ndopP.POLE.rds")) %>% dplyr::select(POLE, DRUH, AUTOR)

print("orig")
print(nrow(ndopP))
print(nrow(ndopP.POLE))


# LSD PA
lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds")) %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci
lsd.pa.min <- readRDS(paste0(path.lsd, "lsd.pa.min.rds")) #  %>% filter(TaxonNameLAT %in% c("Falco tinnunculus","Lanius collurio", "Podiceps cristatus", "Sylvia nisoria", "Acrocephalus arundinaceus", "Accipiter gentilis", "Upupa epops", "Ciconia nigra", "Emberiza citrinella", "Accipiter nisus"))
lsd.pa.min %<>% arrange((TaxonNameLAT)) # %>% filter(TaxonNameLAT %in% c("Accipiter gentilis", "Upupa epops", "Ciconia nigra", "Emberiza citrinella", "Accipiter nisus", "Falco tinnunculus"))
# lsd.pa.min <- lsd.pa.min[1:nrow(lsd.pa.min),]

# lsd.pa.min <- lsd.pa.min[ss[1]:ss[2],]


# odstranit LSD čtverce z prediktorů a NDOP
lsd.pa.centroids.unique <- lsd.pa.centroids %>%
  ungroup() %>%
  group_by(POLE) %>%
  slice_head(n = 1) %>%
  dplyr::select(geometry, POLE) %>%
  ungroup()

# pro kernel
predictors.lsdMask <- raster::mask(predictors, as_Spatial(lsd.pa.centroids.unique), inverse = TRUE)


# dle názvů kvadrátů - nespolehlivé, názvy ne vždy sedí!!!
ndopP.lsdMask2 <- ndopP %>% filter(POLE %notin% lsd.pa.centroids.unique$POLE)
ndopP.POLE.lsdMask2 <- ndopP.POLE %>% filter(POLE %notin% lsd.pa.centroids.unique$POLE)
print("LSD filter")
print(nrow(ndopP.lsdMask2))
print(nrow(ndopP.POLE.lsdMask2))

# prostorově
v.temp <- raster::extract(predictors.lsdMask[[1]], sf::st_coordinates(ndopP), cellnumbers = TRUE)
ndopP.lsdMask <- ndopP
ndopP.lsdMask$cellNumberValue <- v.temp[, 2]
occs.na <- nrow(ndopP.lsdMask %>% filter(is.na(cellNumberValue)))

if (occs.na > 0) {
  print(paste0(occs.na, " occurences removed, no cellnumbers (raster::extract)"))
  ndopP.lsdMask %<>% filter(!is.na(cellNumberValue))
}



v.temp <- raster::extract(predictors.lsdMask[[1]], sf::st_coordinates(ndopP.POLE %>% dplyr::select(geometry) %>% st_centroid()), cellnumbers = TRUE)

ndopP.POLE.lsdMask <- ndopP.POLE
ndopP.POLE.lsdMask$cellNumberValue <- v.temp[, 2]
occs.na <- nrow(ndopP.POLE.lsdMask %>% filter(is.na(cellNumberValue)))

if (occs.na > 0) {
  print(paste0(occs.na, " occurences removed, no cellnumbers (raster::extract)"))
  ndopP.POLE.lsdMask %<>% filter(!is.na(cellNumberValue))
}

print("prostorový filter")
print(nrow(ndopP.lsdMask))
print(nrow(ndopP.POLE.lsdMask))






# odstraním rovnou observery s nízkým počtem nálezů - TODO: normálně bych je měl jen označit NA - nemazat
passiveObservers <- as_tibble(ndopP.lsdMask) %>%
  mutate(uid = row_number()) %>%
  dplyr::select(-geometry) %>%
  ungroup() %>%
  group_by(AUTOR) %>%
  summarise(
    uid.n = n_distinct(uid)
  ) %>%
  filter(uid.n >= 10) %>%
  dplyr::select(AUTOR)
observers.unique <- unique(as.vector(unlist(passiveObservers)))

ndopP.lsdMask %<>% filter(AUTOR %in% observers.unique)
ndopP.POLE.lsdMask %<>% filter(AUTOR %in% observers.unique)

print("autorský filter")
print(nrow(ndopP.lsdMask))
print(nrow(ndopP.POLE.lsdMask))


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! per cca 100m jen jeden nález
# 100m zgrupování - jen body, polygony v ndopP.POLE nelze
r2 <- disaggregate(predictors.lsdMask[[1]], fact = 60) # 27 9
v.temp <- raster::extract(r2, sf::st_coordinates(ndopP.lsdMask), cellnumbers = TRUE)
ndopP.lsdMask$cellNumber2 <- v.temp[, 1]

ndopP.lsdMask %<>% ungroup() %>%
  group_by(AUTOR, cellNumber2, DRUH) %>% # původně bez species, pak to bylo ale jen z hlediska prostorové aktivity, ale bez počtu sesbíraných druhů
  slice_head(n = 1) %>%
  ungroup()


print("po redukci 100m")
print(nrow(ndopP.lsdMask))




ndopP.orig <- ndopP
ndopP.POLE.orig <- ndopP.POLE

ndopP <- ndopP.lsdMask
ndopP.POLE <- ndopP.POLE.lsdMask




# c("Falco tinnunculus","Lanius collurio", "Podiceps cristatus", "Sylvia nisoria", "Acrocephalus arundinaceus", "Accipiter gentilis", "Upupa epops", "Ciconia nigra", "Emberiza citrinella", "Accipiter nisus"))
# statistiky NDOP druhy
ndop.stat.res.sp <- readRDS(paste0(path.ndop, "ndop.stat.res.sp.rds"))


ndopP.POLE.unique <- ndopP.POLE %>%
  group_by(POLE) %>%
  slice_head(n = 1) # %>% dplyr::select(geometry )


#
# kdt <- tgbg::bg(ndopP$geometry, predictors[[1]], authors = ndopP$AUTOR )
#
# kdt.ppp <- kdt$ppp
#
# kdt.ppp.marks <- split(kdt.ppp)
#
#
#
# kdt.ppm <- ppm(kdt.ppp ~  marks , Poisson())
#
#
# kdt.ppm <- ppm(kdt.ppp ~  marks , Poisson()) # * polynom(x,y,3)
#
# ndopP.fs <- ndopP %>% filter(DRUH == "Accipiter gentilis")
#
# kdt.fs <- tgbg::bg(ndopP.fs$geometry, predictors[[1]], authors = ndopP.fs$AUTOR )
#
# kdt.fs.ppp <- kdt.fs$ppp
# kdt.fs.ppp <- unmark(kdt.fs.ppp)
#
# az <- as.im(kdt.ppp.marks$`abrahamek zdenek`)
#
# mr <- as.im( kdt.ppp.marks$`mulacek roman`)
#
#
# kdt.ppp <- unmark(kdt.ppp)
# all <- as.im(kdt.ppp)
#
# kdt.ppm <- ppm(kdt.fs.ppp ~  all , Poisson())
#
# kdt.pred <- predict.ppm(kdt.ppm)


############
# settings
############
# potřebuju být schopný dofiltrovat svou inovaci, null, Botellu, tradiční tgb jako varianty k porovnání
# a být schopen vybrat z různých adjust fittigů 1) prior dobrou variantu (nevím jak to dělá Botella?!) (+ porovnat i post nej validatci LSD)

# Botella et al. (2020):
#   UB – Uniform Bias - uniformly distributed background points (UB).
#   TGB – Target Group Bias - background points as the *sites* where there has been at least one presence among a Target-Group of species (sites je myšlena jedna konkrétní samplovaná lokalita) =“event”/“akce” – nutné znát, jinak jde o TGOB)
#   TGOB - Target-Group Occurrences Background - integrate all species *occurrences* from the Target-Group as background (jednotlivé nálezy, occurrences)
# PB:
#   SSOOG - shared species-observers occurrences group ~ Shared species-observers sampling effort



vn <- versionNames()

vf <- list("1" = 1:8, "0" = 9)


# scénáře výběru bias rasterů k prediktorům
scenario.selected.all <- c("TGOB", "TGOB.sso", "TO", "TO.sso", "TS", "TS.sso", "TO.w", "TO.w.sso", "TS.w", "TS.w.sso", "TS.AO", "TS.AO.sso", "TS.AO.w", "TS.AO.w.sso", "TO.AO", "TO.AO.sso", "TO.AO.w", "TO.AO.w.sso", "un")
scenario.selected.base <- scenario.selected.all[!str_detect(scenario.selected.all, "sso|w|AO")]
scenario.selected.base.w <- scenario.selected.all[str_detect(scenario.selected.all, "\\.w|TGOB")]


scenarios.selected.puv <- list(
  "1" = "TGOB",
  "2" = "TS.w",
  "3" = "TO.w",
  "4" = "TS.AO.w",
  "5" = "TO.AO.w",
  "1.sso" = "TGOB.sso",
  "2.sso" = "TS.w.sso",
  "3.sso" = "TO.w.sso",
  "4.sso" = "TS.AO.w.sso",
  "5.sso" = "TO.AO.w.sso",
  "6" = scenario.selected.base.w
)

scenarios.selected.puv2 <- list(
  "1" = "TGOB",
  "2" = "TS.w",
  "3" = "TO.w",
  "4" = "TGOB.sso.w",
  "5" = c("TGOB", "TS.w", "TO.w", "TGOB.sso.w")
)

scenarios.selected.puv8 <- list(
  "1" = "TGOB",
  "2" = "TS.w",
  "3" = "TO.w",
  "4" = "TGOB.sso.w"
)

scenarios.selected.puv9 <- list(
  "1" = "TGOB",
  "2" = "TS.w",
  "3" = "TO.w",
  "4" = "TGOB.sso.w",
  "5" = "TGOB.sso.w.dpois.mean",
  "6" = "TGOB.sso.w.dpois.median",
  "7" = "TGOB.sso.w.dnorm.mean",
  "8" = "TGOB.sso.w.dnorm.median"
)

scenarios.selected <- list(
  "1" = "TGOB",
  "2" = "TS.w",
  "3" = "TO.w",
  "4" = "TS.AO.w",
  "5" = "TO.AO.w",

  "6" = "TGOB.sso",
  "7" = "TS.sso",
  "8" = "TO.sso",
  "9" = "TS.AO.sso",
  "10" = "TO.AO.sso",

  "11" = "TS.w.sso",
  "12" = "TO.w.sso",
  "13" = "TS.AO.w.sso",
  "14" = "TO.AO.w.sso",
  #
  # "1.sso" = "TGOB.sso",
  # "2.sso" = "TS.w.sso",
  # "3.sso" = "TO.w.sso",
  # "4.sso" = "TS.AO.w.sso",
  # "5.sso" = "TO.AO.w.sso",
  #
  "1.sso.w" = "TGOB.sso.w",
  "1.sso.w.3" = "TGOB.sso.w.3",
  #
  "1.sso.w2" = "TGOB.sso.w2",
  "1.sso.w2.3" = "TGOB.sso.w2.3",
  #
  "1.sso.w.dnorm2.mean" = "TGOB.sso.w.dnorm2.mean",
  "1.sso.w.3.dnorm2.mean" = "TGOB.sso.w.3.dnorm2.mean",
  "1.sso.w.dnorm2.median" = "TGOB.sso.w.dnorm2.median",
  "1.sso.w.3.dnorm2.median" = "TGOB.sso.w.3.dnorm2.median",
  #
  "1.sso.w.dnorm3.mean" = "TGOB.sso.w.dnorm3.mean",
  "1.sso.w.3.dnorm3.mean" = "TGOB.sso.w.3.dnorm3.mean",
  #
  "1.sso.w2.dnorm3.mean" = "TGOB.sso.w2.dnorm3.mean",
  "1.sso.w2.3.dnorm3.mean" = "TGOB.sso.w2.3.dnorm3.mean",
  "1.sso.w2.dnorm3.median" = "TGOB.sso.w2.dnorm3.median",
  "1.sso.w2.3.dnorm3.median" = "TGOB.sso.w2.3.dnorm3.median"
)


ndop.fs <- list(
  "adjusts" = c(0.25, 0.5, 1, 2, 3, 4), # 0.1, 1, 2
  "tuneArgs" = list(fc = c("L", "Q", "LQ"), rm = c(1, 1.5, 2, 4, 6, 10)), # light, full: "L", "LQ", "LQP", "LQH", "LQHP" c(0.5, 1, 2, 3, 4)
  "bgRatio" = 1 / 5,
  "bg" = 5000,
  "speciesPerGroup" = 2, "speciesOccMin" = 30,
  "sq2rad" = c((1 / 6) / 4, 0.1 / 4), # kvadráty KFME 2rad, xy velikost ve stupních
  "sq2radDist" = c(1:5), # c(1:5)
  "replicates" = 10,
  "speciesPart" = cmd_arg, "version" = "v1",
  "bgRaster" = FALSE,
  "versionNames" = vn, "versionSmooting" = vf,
  "outputBgBr" = c("bg", "br"), # c("bg", "br")
  "scenario" = "bgAll", # "bg", "br", "brAll", "bgAll"
  "scenarios" = scenarios.selected
)


############
# execution
############

# odstranit LSD sites z rastru prediktorů - jen pro generování tgob/ssoob - možná až na finální, teď ne
# odstraním prezence (sites) LSD z NDOP pro zabráníení spatial (auto)correlation - LSD tak bude i prostorově nezávislý dataset




evalIndep <- function(m, lsd.temp, preds) {
  print("evalIndep!")
  first <- TRUE
  out.t <- NA
  for (layer in names(m@predictions)) { # v případě rozdělení per RM+FC je cyklus zbytečný (ale preventivně nechat), je tam jen jeden RasterLayer
    ev <- NA
    r.temp <- m@predictions[[layer]]

    # r.temp <- enm.maxnet@predict(m@models[[1]], preds, list(doClamp = FALSE, pred.type = "cloglog"))


    # @models[[1]] - také jednotlivé modely???

    ex.predicted <- extract(r.temp, st_coordinates(lsd.temp))
    #
    # dopočíst a přidat další metriky z performance() (rgee)!!!
    # tam ale není nic co zohledňuje tn tp (bez poměru k fp)
    #
    pr <- list()
    pr$auc.integral <- NA


    pr <- pr.curve(ex.predicted[lsd.temp$presence == 1], ex.predicted[lsd.temp$presence == 0])

    if (inherits(try({
      ev <- sdmATM::evaluates(lsd.temp$presence, ex.predicted)
    }), "try-error")) {
      ev <- NA
    }
    if (length(ev) == 0) {
      print("evalIndep! - failed+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
      next
    }
    ev.temp <- as.data.frame(ev@statistics[-3])
    ev.temp <- merge(ev.temp, t(as.data.frame(ev@statistics$COR)))
    ev.temp <- merge(ev.temp, ev@threshold_based[1, ]) # sp=se
    ev.temp[["tune.args"]] <- layer

    ev.temp[["prc.integral"]] <- pr$auc.integral




    if (first) {
      first <- FALSE
      out.t <- ev.temp
    } else {
      out.t %<>% add_row(ev.temp)
    }
  }

  return(out.t)
}


evalDep <- function(m, lsd.temp, preds, bias, p, nw) {

  # nw <- NA # zrusit docasne nefungujew
  print("evalDep!")


  first <- TRUE
  out.t <- NA
  for (layer in names(m@predictions)) { # v případě rozdělení per RM+FC je cyklus zbytečný (ale preventivně nechat), je tam jen jeden RasterLayer
    ev <- NA
    r.temp <- m@predictions[[layer]]

    # r.temp <- enm.maxnet@predict(m@models[[1]], preds, list(doClamp = FALSE, pred.type = "cloglog"))

    ex.predicted <- extract(r.temp, st_coordinates(lsd.temp))



    cor.manual <- list()
    cor.manual$estimate <- NA
    cor.manual$p.value <- NA

    cor.manual2 <- list()
    cor.manual2$estimate <- NA
    cor.manual2$p.value <- NA


    pr <- list()
    prw <- list()
    aucw <- list()
    pr$auc.integral <- NA
    prw$auc.integral <- NA
    aucw$auc.integral <- NA

    pr <- pr.curve(ex.predicted[lsd.temp$presence == 1], ex.predicted[lsd.temp$presence == 0])

    if (length(bias) > 1) {
      # korelace pro evaluaci
      v.temp <- raster::extract(bias, st_coordinates(p))

      v.temp.reverse <- v.temp * -1 + max(v.temp) + min(v.temp)

      prw <- pr.curve(ex.predicted[lsd.temp$presence == 1], ex.predicted[lsd.temp$presence == 0], nw, rep(0, length(ex.predicted[lsd.temp$presence == 0])))

      aucw <- roc.curve(ex.predicted[lsd.temp$presence == 1], ex.predicted[lsd.temp$presence == 0], nw, rep(0, length(ex.predicted[lsd.temp$presence == 0])))

      # prw <- pr.curve(ex.predicted[lsd.temp$presence == 1], ex.predicted[lsd.temp$presence == 0])
      #
      # aucw <- roc.curve(ex.predicted[lsd.temp$presence == 1], ex.predicted[lsd.temp$presence == 0])




      #
      # tutéžř korekci outlierů???
      #
      bp.min <- min(max(v.temp), quantile(v.temp, prob = 0.75) + 1.5 * IQR(v.temp))
      v.temp[v.temp >= bp.min] <- bp.min
      v.temp.reverse2 <- v.temp * -1 + max(v.temp) + min(v.temp)


      p.predicted <- extract(r.temp, st_coordinates(p))

      cor.manual <- cor.test(v.temp.reverse, p.predicted)
      cor.manual2 <- cor.test(v.temp.reverse2, p.predicted)
    }
    #
    # dopočíst a přidat další metriky z performance() (rgee)!!!
    # tam ale není nic co zohledňuje tn tp (bez poměru k fp)
    #

    if (inherits(try({
      ev <- sdmATM::evaluates(lsd.temp$presence, ex.predicted)
    }), "try-error")) {
      ev <- NA
    }
    if (length(ev) == 0) {
      #  print(lsd.temp$presence)
      # print(ex.predicted)
      print("nemám evaluaci!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      next
    }
    ev.temp <- as.data.frame(ev@statistics[-3])
    ev.temp <- merge(ev.temp, t(as.data.frame(ev@statistics$COR)))
    ev.temp <- merge(ev.temp, ev@threshold_based[1, ]) # sp=se
    ev.temp[["tune.args"]] <- layer
    ev.temp[["cor.manual.estimate"]] <- cor.manual$estimate
    ev.temp[["cor.manual.p.value"]] <- cor.manual$p.value
    ev.temp[["cor.manual2.estimate"]] <- cor.manual2$estimate
    ev.temp[["cor.manual2.p.value"]] <- cor.manual2$p.value


    ev.temp[["prc.integral_trainRandom"]] <- pr$auc.integral
    ev.temp[["prcw.integral_trainRandom"]] <- prw$auc.integral
    ev.temp[["aucw.integral_trainRandom"]] <- aucw$auc.integral

    if (first) {
      first <- FALSE
      out.t <- ev.temp
    } else {
      out.t %<>% add_row(ev.temp)
    }
  }

  return(out.t)
}




evalDepFolds <- function(m, lsd.temp, folds, preds) {
  first <- TRUE
  out.t <- NA




  for (layer in names(m@predictions)) { # v případě rozdělení per RM+FC je cyklus zbytečný (ale preventivně nechat), je tam jen jeden RasterLayer
    ev <- NA
    r.temp <- m@predictions[[layer]]

    # r.temp <- enm.maxnet@predict(m@models[[1]], preds, list(doClamp = FALSE, pred.type = "cloglog"))

    ex.predicted <- extract(r.temp, st_coordinates(lsd.temp))
    #
    # dopočíst a přidat další metriky z performance() (rgee)!!!
    # tam ale není nic co zohledňuje tn tp (bez poměru k fp)
    #


    for (fid in unique(folds)) {
      print(fid)



      lsd.temp.folds <- lsd.temp[which(folds == fid), ]
      ex.predicted.folds <- ex.predicted[which(folds == fid)]

      # print(lsd.temp.folds)
      # print(ex.predicted.folds)


      # return(list())


      if (inherits(try({
        ev <- sdmATM::evaluates(lsd.temp.folds$presence, ex.predicted.folds)
      }), "try-error")) {
        print("nemám evaluaci!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        ev <- NA
      }
      if (length(ev) == 0) {
        next
      }
      ev.temp <- as.data.frame(ev@statistics[-3])
      ev.temp <- merge(ev.temp, t(as.data.frame(ev@statistics$COR)))
      ev.temp <- merge(ev.temp, ev@threshold_based[1, ]) # sp=se
      # ev.temp[["tune.args"]] <- layer
      if (first) {
        first <- FALSE
        out.t <- ev.temp
      } else {
        out.t %<>% add_row(ev.temp)
      }
    }
  }

  out.t %<>% summarise_if(is.numeric, mean, na.rm = TRUE)
  out.t[["tune.args"]] <- layer


  return(out.t)
}

sumNormal <- function(x) {
  return(x / sum(as.vector(x), na.rm = TRUE))
}


# kešovat
fn <- paste0(path.tgob, "tgobV.rds")
# if (file.exists(fn)) {
#     print("načítám existující tgobV.rds")
#     # tgobV <- readRDS(fn)
#
#     # # dočasně, nepoužívám AO - odstranit, ať to nežere místo
#     # tgbg.names  <-  names(tgobV$p)
#     # tgobV$p %<>% dplyr::select(tgbg.names[!str_detect(tgbg.names, "\\.AO")])
#     # tgbg.names  <-  names(tgobV$p)
#     # tgobV$p %<>% dplyr::select(tgbg.names[!str_detect(tgbg.names, "sso_")])
#     # saveRDS(tgobV, paste0(path.tgob, "tgobV2.rds"))
#     tgobV <- readRDS("D:/PersonalWork/Balej/v2/RP/dataPrep/ndopTGOB/tgobV.rds")
# } else {
#     print("generuju nový tgobV.rds")
#    tgobV <- tgbg::tgobVersions(ndopP, predictors.lsdMask[[1]], species = "DRUH", observers = "AUTOR", species.select = lsd.pa.min$TaxonNameLAT, path = path.tgob.tgbg)


# tgobV <- tgbg::tgobVersions(ndopP, predictors.lsdMask[[1]], species = "DRUH", observers = "AUTOR", species.select = lsd.pa.min$TaxonNameLAT, path = path.tgob.tgbg, pp=ndopP.POLE.unique)


# # tgobV <- tgbg::tgobVersions(ndopP, predictors3, species = "DRUH", observers = "AUTOR", species.select = lsd.pa.min$TaxonNameLAT, path = path.tgob.tgbg)

#    saveRDS(tgobV, fn)
# }


#  # rozdělím per species tgobV - optimalizace
# cnames <- names(tgobV$p)
# for(sp in unique(lsd.pa.min$TaxonNameLAT)){ # unique(lsd.pa.min$TaxonNameLAT) unique(tgobV$p$DRUH)
#
#        cnames.c <- c(cnames[str_detect(cnames, paste0(sp, "$"))], "DRUH", "geometry", "nc_TO.w", "nc_TS.w", "nc_TO", "nc_TS") # paste0("nc_sso_", sp)
#
#        # !!sym(cnames.c))
# #Please use `all_of()` or `any_of()` instead. - namísto přímého vstupu text řetezců...
# saveRDS(tgobV$p  %>% dplyr::select(cnames.c), paste0(path.tgob.tgbg, sp,".rds"))
#
# }



# tgobV <- tgbg::tgobVersions(ndopP, predictors.lsdMask[[1]], species = "DRUH", observers = "AUTOR", species.select = lsd.pa.min$TaxonNameLAT, path = path.tgob.tgbg)
# Sys.sleep(30)
# stop()


prefix <- "nc_"

# # nepotřebuji při použití výše
# cn <- names(tgobV$p)
# cn_prefix <- cn[str_detect(cn, paste0("^", prefix))]

ndop.stat.res.sp.selected <- ndop.stat.res.sp %>%
  filter(DRUH %in% unique(lsd.pa.min$TaxonNameLAT)) %>%
  ###  filter(DRUH %in% c("Accipiter gentilis", "Emberiza citrinella")) %>%
  ###  filter(DRUH %in% c("Accipiter nisus")) %>%
  arrange(DRUH) ### arrange(POLE.n)
sp.diff <- setdiff(unique(lsd.pa.min$TaxonNameLAT), ndop.stat.res.sp.selected$DRUH)
if (length(sp.diff) > 0) {
  print("Problém s nejednotnou taxonomií nebo neprůnikem LSD druhů s NDOP druhy!")
}

# fronta

que <- function(replicates, ndop.stat.res.sp.selected, path.tgob) {
  out <- list("species" = NA, "replicates" = NA)
  que.fn <- "que.rds"
  if (file.exists(paste0(path.tgob, que.fn))) {
    print("načítám existující que.rds")
    que <- readRDS(paste0(path.tgob, que.fn))
  } else {
    print("generuju nový que.rds")
    que <- ndop.stat.res.sp.selected %>% arrange(DRUH) ### arrange(POLE.n)
    que$rep <- 0
    que$lock <- 0
    saveRDS(que, paste0(path.tgob, que.fn))
  }


  locked <- nrow(que %>% filter(lock == 1))
  que.r <- que %>%
    filter(rep < replicates) %>%
    arrange(DRUH) %>% ### arrange(POLE.n)
    slice_head(n = 1)


  if (locked > 0) {
    if (nrow(que.r) > 0) {
      print("table locked, waiting...")
      Sys.sleep(1)
      que(replicates, ndop.stat.res.sp.selected, path.tgob)
    } else {
      print("no tasks in que, finised")
      return(out)
    }
  } else {
    print("selecting next combination of species & replication, lock and save table")
    que$lock <- 1
    saveRDS(que, paste0(path.tgob, que.fn))


    if (nrow(que.r) > 0) {
      print("next combination of species & replication selected, unlock and save updated table, return species & replication")
      que.sr <- que %>%
        filter(rep < replicates) %>%
        arrange(DRUH) %>% ### arrange(POLE.n)
        slice_head(n = 1) %>%
        mutate(rep = rep + 1)
      print(que.sr)
      que <- dplyr::rows_update(que, que.sr)
      que$lock <- 0
      saveRDS(que, paste0(path.tgob, que.fn))
      return(list("species" = as.character(unname(unlist(que.r$DRUH))), "replicates" = as.numeric(unname(unlist(que.r$rep))) + 1))
    } else {
      print("no rows to select")
      return(out)
    }
  }
}


# druhy s malým počtem presencí se počítají mnohem rychleji, chci tyto per group upřednostnit (seřadit rovnoměrně v rámci skupin druhy od nejméně početných)
ndop.stat.res.sp.selected %<>% mutate(nthGroup = NA)
rowsTotal <- nrow(ndop.stat.res.sp.selected)
groups <- ceiling(rowsTotal / ndop.fs$speciesPerGroup)

first <- TRUE
for (nth in 1:groups) {
  fill <- seq(nth, rowsTotal, groups)
  temp.rows <- ndop.stat.res.sp.selected %>% filter(row_number() %in% fill)
  temp.rows %<>% mutate(nthGroup = nth)
  if (first) {
    first <- FALSE
    ndop.stat.res.sp.selected.regrouped <- temp.rows
  } else {
    ndop.stat.res.sp.selected.regrouped %<>% add_row(temp.rows)
  }
}

#
# run per species models and varying BG
#
e.mx.all <- list()
e.mx.all.new <- list()
tune.args <- ndop.fs$tuneArgs
ll <- c("longitude", "latitude")
druhy <- unique(as.vector(unlist(ndop.stat.res.sp.selected$DRUH)))



# "tradiční" TGOB - použít jako BG všechny unikátní per_pixel presence všech druhů (NDOP pokrývá cca 85% území s cca 7300 pixely (z 9700)), takže je to zde upočitatelné

tgob.trad.POLE <- ndopP.POLE %>%
  group_by(POLE) %>%
  slice_head(n = 1) %>%
  st_centroid() %>%
  dplyr::select(POLE)

tgob.trad <- tgob.trad.POLE %>% dplyr::select(-everything())

length(unique(tgob.trad.POLE$POLE))


# "tradiční" TGOB per pole a druh - pro thin
tgob.trad.sp <- ndopP.POLE %>%
  group_by(POLE, DRUH) %>%
  slice_head(n = 1) %>%
  st_centroid()

sp.group <- ndop.stat.res.sp.selected.regrouped %>% filter(nthGroup == ndop.fs$speciesPart)






# ###
# # per autor/bias raster
# ###
#
#
# fn <- paste0(path.tgob, "perAutorBr.rds")
# if (file.exists(fn)) {
#     print("načítám existující perAutorBr.rds")
#
#   collector.aut <- readRDS(fn)
# } else {
#     print("generuju nový perAutorBr.rds")
#
#
#   ndop.tgob <- readRDS("/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/dataPrep/ndopTGOBtgbg/Corvus frugilegus.rds") %>% dplyr::select(nc_AUTOR)
#
#   autors <- as_tibble(ndop.tgob) %>% dplyr::select(-geometry) %>% distinct(nc_AUTOR)
#
#   collector.aut <- list()
#   collector.cnt <- list()
#   for (au in as.vector(autors$nc_AUTOR)) {
#     p.temp <- NA
#     p.temp <- ndop.tgob %>% filter(nc_AUTOR == au)
#     print(au)
#
#     cnt <- nrow(p.temp)
#     print(cnt)
#     # collector[[au]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), predictors[[1]], sigma = ndop.fs$adjusts, weights = p.temp.inf.new, output = ndop.fs$outputBgBr, anisotropic = TRUE)
#     collector.aut[[au]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), predictors[[1]], sigma = ndop.fs$adjusts, output = "br", anisotropic = TRUE)
#
#     collector.cnt[[au]] <- cnt
#     # jen přidat váhy jednotlivým rasterům podle w3? (pronásobit je?)
#
#   }
#
#   saveRDS(collector.aut, fn)
#   saveRDS(collector.cnt, paste0(path.tgob, "perAutorBr.cnt.rds"))
# }

###







repeat {
  rs <- que(ndop.fs$replicates, ndop.stat.res.sp.selected, path.tgob)

  if (sum(is.na(rs)) > 0) {
    # obsahuje druh a počet replikací? - bez toho ukončím modelování
    break
  }

  rep <- rs$replicates
  druh <- rs$species
  print("RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR")
  print(rep)
  set.seed(rep)


  # https://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/inst/doc/tutorial_2.html
  # range <- cv_spatial_autocor(
  #     x = pa_data, # species data - jen PA????
  #     column = "occ", # column storing presence-absence records (0s and 1s)
  #     plot = FALSE
  # )

  # range$range

  # scv2 <- cv_nndm(
  #     x = pa_data,
  #     column = "occ",
  #     r = rasters,
  #     size = 360000, # range of spatial autocorrelation
  #     num_sample = 10000, # number of samples of prediction points
  #     sampling = "regular", # sampling methods; it can be random as well
  #     min_train = 0.1, # minimum portion to keep in each train fold
  #     plot = TRUE,
  #   presence_bg = FALSE  #nově pro p-bg
  # )
  # # see the number of folds in scv2 object
  # scv2$k

  # # Folds generated by cv_nndm function are used here (a training and testing fold for each record) to show how to use folds
  # # from this function (the cv_buffer is also similar to this approach) for evaluation species distribution models.
  # # Note that with cv_nndm using presence-absence data (and any other type of data except for presence-background data when presence_bg = TRUE is used),
  # # there is only one point in each testing fold, and therefore AUC cannot be calculated for each fold separately. Instead, the value of
  # # each point is first predicted to the testing point (of each fold), and then a unique AUC is calculated for the full set of predictions.

  # bloky (umožnil bych automatickou volbou velikosti bloku zvýšit korelaci mezi auc.val.avg a AUC)
  bCV <- blockCV::cv_spatial(st_as_sf(rasterToPoints(predictors[[1]], spatial = TRUE)) %>% dplyr::select(-everything()), size = 50000, deg_to_metre = 90000, k = 5, selection = "random", hexagon = TRUE, seed = rep, plot = FALSE)
  bCV.poly <- as_tibble(bCV[["blocks"]][["geometry"]])
  bCV.poly$fold <- bCV[["blocks"]][["folds"]]
  bCV.poly <- st_as_sf(bCV.poly)

  # výchozí random BG společný pro všechny # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # null.default <- generateRPall(predictors[[1]], nBackRatio = ndop.fs$bgRatio)

  ### nové null ze sample()?
  null.default <- tgbg::sample(predictors[[1]])

  #  for (druh in sp.group$DRUH) { # for DRUH
  first <- TRUE # pro každou replikaci a druh samostatný soubor
  collector <- list()
  collector.thin <- list()
  print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print(druh)
  print(sp.group)
  print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  e.mx.all <- list()
  e.mx.all.new <- list()
  out.t <- list()
  bg.col <- list()
  # LSD
  lsd.temp <- NA
  lsd.temp <- lsd.pa.centroids %>% filter(TaxonNameLAT == druh)
  #
  # základní verze backgroundů
  #

  # # TGOB ssos
  # ssos.sp <- cn_prefix[str_detect(cn_prefix, paste0(druh, "$"))]
  #
  # versions.base.puv <- c("TGOB", "TO", "TS", "TO.w", "TS.w")
  # versions.base <- c("TGOB", "TO.w", "TS.w")
  versions.base <- c("TGOB")
  # ao.sp <- ssos.sp[!str_detect(ssos.sp, "sso_")]
  # ao.sp.puv <- gsub(prefix, "", ao.sp)
  # ao.sp <- c()
  # ssos.w <- ssos.sp[str_detect(ssos.sp, "sso\\.w")]
  # ssos.w <- gsub(prefix, "", ssos.w)
  #
  # ssos.sp <- ssos.sp[str_detect(ssos.sp, "sso_")]
  # # ssos.sp <- gsub(prefix, "", ssos.sp)



  # QQQ přepsat nepotřebuju to tahat odsud? Kvůli autorům ano, ale to se řešilo jen pro kernely, tady už beru jen presence!!!


  tgobV <- list()

  # tgobV$p <- ndopP

  tgobV$p <- readRDS(paste0(path.tgob.tgbg, druh, ".rds"))

  tgobV$p %<>% rename(DRUH = DRUH.x)


  # tgobV$p %<>% dplyr::select( !matches("\\.w2z\\_|\\.w2y\\_|\\.w\\_|\\.w46\\_|\\.w99\\_|\\.w8\\_|\\.w2\\_|\\.w2x\\_|\\.w0\\_|\\.w0s\\_" ))

  # náhodně vyberu jeden bod per obs/px/sp
  # tgobV$p %<>% ungroup() %>%
  # group_by(nc_AUTOR, nc_cellNumber, DRUH) %>%
  # slice_sample(n = 1)




  tgobVnames <- names(tgobV$p)

  #
  for (vn in c(versions.base)) { # c(versions.base, ao.sp, ssos.w)//// "\\.w" |maxSE1_    tgobVnames[str_detect(tgobVnames, "\\.w23fa_|TT4_|TT2_|TT6_|TT8_|TT88_|TT9_|TTT")]
    print("----------------------------------------------------------------------------------")
    print("----------------------------------------------------------------------------------")
    print("----------------------------------------------------------------------------------")
    print(vn)

    if (vn == "TGOB" | vn == "TS" | vn == "TO" | vn == "TS.AO" | vn == "TO.AO") {
      vns <- vn
    } else {
      vns <- strsplit(vn, "_")[[1]][2] # odstraním případnou příponu s druhem -- posun na 2
    }
    print(vns)

    print("----------------------------------------------------------------------------------")
    print("----------------------------------------------------------------------------------")
    print("----------------------------------------------------------------------------------")

    p.temp <- NA


    print("odstraním focus species (pro generování bias rasterů) - vede k overprediction") # zrušeno (teď pro bg scénář)
    print(druh)
    print(nrow(tgobV$p %>% filter(DRUH == druh)))

    # základní verze
    if (vn == "TGOB") {
      p.temp <- tgobV$p %>%
        # filter(DRUH != druh) %>%
        dplyr::select(geometry)
      # zkresluju podle druhu
      p.temp.bss <- tgobV$p %>%
        filter(DRUH == druh)
    } else {
      if (str_detect(vn, "\\.w")) {

        # beru všechny Top, ale musím odstranit NA, bez vah
        p.temp <- tgobV$p %>%
          # filter(DRUH != druh) %>%
          filter(!is.na(!!sym(paste0(vn))))
      } else {
        if (vn == "TS" | vn == "TO" | vn == "TT1" | vn == "TT2" | vn == "TT3" | vn == "TT4" | vn == "TT5" | vn == "TT6" | vn == "TT7" | vn == "TT8" | vn == "TT88" | vn == "TT9" | vn == "TT0" | vn == "TT10" | vn == "TT11" | vn == "TT20" | vn == "TT21" | vn == "TT22" | vn == "TT23" | vn == "maxSE1") {
          p.temp <- tgobV$p %>%
            # filter(DRUH != druh) %>%
            filter(!!sym(paste0(prefix, vn)) == 1)
        } else {
          if (str_detect(vn, "TTT")) {
            p.temp <- tgobV$p %>%
              # filter(DRUH != druh) %>%
              filter(!!sym(paste0(prefix, vns, "_", druh)) == 1)
          }


          if (vn == "TS.AO" | vn == "TO.AO") {
            p.temp <- tgobV$p %>%
              # filter(DRUH != druh) %>%
              filter(!!sym(paste0(prefix, vn, "_", druh)) == 1)
          } else {
            p.temp <- tgobV$p %>%
              # filter(DRUH != druh) %>%
              filter(!!sym(paste0(vn)) == 1)
          }
        }
      }
    }

    if (nrow(p.temp) < 1) {
      print("neexistují presence základní verze!!!")
      next
    }

    if (str_detect(vn, "\\.w")) {
      p.temp.inf.new <- p.temp[[paste0(vn)]]
      if (str_detect(vn, "\\.w5i")) {
        p.temp.inf <- p.temp.inf.new > 1.0
        p.temp.inf.new[p.temp.inf] <- 1
      }


      collector[[vns]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), predictors.lsdMask[[1]], sigma = ndop.fs$adjusts, weights = p.temp.inf.new, output = ndop.fs$outputBgBr, anisotropic = TRUE)


      # saveRDS(p.temp.inf.new, paste0(path.tgob, "35ssos_collectorW_", vns, "_", druh, "_", rep, ".rds"))



      if (vns == "TS.w" | vns == "TO.w" | vns == "TS.AO.w" | vns == "TO.AO.w") {
        # extra sso zde
        print("___________________________________________________________________________________________________________________________________________________________________")
        print("extra sso zde")
        print(vns)
        print(vn)
        print(paste0(prefix, "sso_", druh))
        print(paste0(vns, ".sso"))
        print("___________________________________________________________________________________________________________________________________________________________________")

        p.temp.temp <- p.temp %>% filter(!!sym(paste0(prefix, "sso_", druh)) == 1)

        collector[[paste0(vns, ".sso")]] <- tgbg::bg(p.temp.temp %>% dplyr::select(geometry), predictors[[1]], sigma = ndop.fs$adjusts, weights = p.temp.temp[[paste0(vn)]], output = ndop.fs$outputBgBr, anisotropic = TRUE)
        # abs(p.temp.temp[[paste0(vn)]] - 1)
        # + zkombinovat i w druhu a w sso???? TS.w.sso.w?
      }
    } else {
      if (str_detect(vn, "TGOB")) {
        collector[["TGOB01"]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), n = 0.02, predictors.lsdMask[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)
        collector[["TGOB02"]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), n = 0.03, predictors.lsdMask[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)

        collector[["TGOB03"]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), n = 0.05, predictors.lsdMask[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)

        collector[["TGOB04"]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), n = 0.1, predictors.lsdMask[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)

        collector[["TGOB05"]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), n = 0.2, predictors.lsdMask[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)

        collector[["TGOB06"]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), n = 0.3, predictors.lsdMask[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)
      } else {
        collector[[vns]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), predictors.lsdMask[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)
      }




      #             if (vn == "TGOB") {
      #                 # zkresluju podle druhu
      #                 temp.bss.revert <- tgbg::bg(p.temp.bss %>% dplyr::select(geometry), predictors[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)
      #
      # temp.bss.revert.temp <- list()
      # for(adjR in names(temp.bss.revert[["br"]])){
      #
      # temp.bss.revert.temp[["bg"]][[adjR]] <- tgbg::sample(temp.bss.revert[["br"]][[adjR]], prob=TRUE)
      # temp.bss.revert.temp[["br"]][[adjR]] <- temp.bss.revert[["br"]][[adjR]]
      # }
      # 		   collector[[paste0(vns, ".bss")]] <- temp.bss.revert.temp
      #             }
    }


    # #
    # # sso
    # #
    # if (vn == "TGOB" | vn == "TS" | vn == "TO" | vn == "TS.AO" | vn == "TO.AO") {
    # #### # základní čisté nevážené sso druhově nespecifické
    # if (vn == "TGOB") {
    # p.temp <- tgobV$p %>%
    # # filter(DRUH != druh) %>%

    # filter(!!sym(paste0(prefix, "sso_", druh)) == 1) %>%
    # dplyr::select(geometry)
    # }

    # if (vn == "TS" | vn == "TO") {
    # p.temp <- tgobV$p %>%
    # # filter(DRUH != druh) %>%
    # filter(!!sym(paste0(prefix, vn)) == 1) %>%
    # filter(!!sym(paste0(prefix, "sso_", druh)) == 1) %>%
    # dplyr::select(geometry)
    # }

    # if (vn == "TS.AO" | vn == "TO.AO") {
    # p.temp <- tgobV$p %>%
    # # filter(DRUH != druh) %>%
    # filter(!!sym(paste0(prefix, vn, "_", druh)) == 1) %>%
    # filter(!!sym(paste0(prefix, "sso_", druh)) == 1) %>%
    # dplyr::select(geometry)
    # }

    # if (nrow(p.temp) < 1) {
    # print("neexistují presence základní verze!!!")
    # next
    # }

    # collector[[paste0(vns, ".sso")]] <- tgbg::bg(p.temp, predictors[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)
    # }



    # if (!str_detect(vn, "sso\\.w")) {
    # # ssos podverze ze základních, nedělám ale z už vážených sso, tam je ten výběr proveden už při vzniku
    # p.temp <- NA
    # vns <- paste0(vns, ".", strsplit(ssos.sp, "_")[[1]][2]) # přidám příponu sso
    # print("sso ------------")
    # print(vns)

    # if (vn == "TGOB") {
    #     p.temp <- tgobV$p %>%
    #         filter(!!sym(ssos.sp) == 1) %>%
    #         dplyr::select(geometry)
    # } else {
    #     if (str_detect(vn, "\\.w")) {
    #         # beru všechny Top, ale musím odstranit NA, bez vah
    #         p.temp <- tgobV$p %>%
    #             filter(!is.na(!!sym(paste0(prefix, vn)))) %>%
    #             filter(!!sym(ssos.sp) == 1) %>%
    #             dplyr::select(geometry, !!sym(paste0(prefix, vn)))
    #     } else {
    #         p.temp <- tgobV$p %>%
    #             filter(!!sym(paste0(prefix, vn)) == 1) %>%
    #             filter(!!sym(ssos.sp) == 1) %>%
    #             dplyr::select(geometry)
    #     }
    # }
    # if (nrow(p.temp) < 1) {
    #     print("neexistují presence základní verze!!!")
    #     next
    # }

    # if (str_detect(vn, "\\.w")) {
    #     collector[[vns]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), predictors[[1]], sigma = ndop.fs$adjusts, weights = p.temp[[paste0(prefix, vn)]], output = ndop.fs$outputBgBr, anisotropic = TRUE)
    # } else {
    #     collector[[vns]] <- tgbg::bg(p.temp, predictors[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)
    # }

    # }
  }


  gc()





























  pres <- ndopP %>% filter(DRUH == druh)

  pres.unique <- ndopP.POLE %>%
    filter(DRUH == druh) %>%
    group_by(POLE) %>%
    slice_head(n = 1) %>%
    st_centroid() %>%
    dplyr::select(-everything())

  pres.unique.poly <- ndopP.POLE %>%
    filter(DRUH == druh) %>%
    group_by(POLE) %>%
    slice_head(n = 1) %>%
    dplyr::select(-everything())

  # presence do modelu a k thinningu
  df.temp <- as.data.frame(st_coordinates(pres.unique))
  names(df.temp) <- ll

  # pro thinning - rozsahy
  sq2rad.dist <- ndop.fs$sq2rad[1] + 0.00001 # delší strana kvadrátu, přičtení drobné vzdálenosti
  sq2rad.dist.range <- sq2rad.dist * ndop.fs$sq2radDist

  # # buffer kolem presencí
  # # náhodně může vznikat více typů geometrií, pak to spadne do sfc_GEOMETRY, to nechci (neumí s tím pak pracovat dále některé funkce), vybírám pouze (MULTI)POLYGON
  # # buffer je 4*strana 2rad, od středu kvadrátu tedy cca 12 km, reálně mi jde o to získat buffer cca 3 čtverců okolo sousedících
  # pres.unique.buffer <- st_as_sf(st_union(st_buffer(pres.unique, 10000))) %>% filter(st_geometry_type(x) %in% c("MULTIPOLYGON", "POLYGON"))
  # predictors.ssos.buffer <- mask(predictors[[1]], pres.unique.buffer)


  tgob.trad.tm <- tgob.trad.sp %>% filter(DRUH == druh)
  tgob.trad.tm.df <- as.data.frame(st_coordinates(tgob.trad.tm))
  colnames(tgob.trad.tm.df) <- c("x", "y")

  # # # buffer
  # id <- paste(c("area", "buffer"), collapse = "_")
  # collector.temp <- list()
  # collector.temp[["bg"]][["0"]] <- st_as_sf(rasterToPoints(kss[["buffer"]], spatial = TRUE)) %>% dplyr::select(-everything()) # vizuální kontrola!!!
  # collector[[id]] <- collector.temp

  #
  # un je fix sám o sobě, přidám
  #
  # id <- paste(c("cz", "un"), collapse = "_")
  id <- "un"
  collector[[id]][["bg"]][["0"]] <- null.default

  saveRDS(collector, paste0(path.tgob, "55ssos_collector_", druh, "_", rep, ".rds"))



  # thinning  presencí pro UN
  for (thinDist in sq2rad.dist.range) {
    occs.thinned <- ecospat.occ.desaggregation(xy = tgob.trad.tm.df, min.dist = thinDist, by = NULL)
    thinDist <- as.character(thinDist)
    collector.thin[[id]][[thinDist]] <- occs.thinned %>% st_as_sf(coords = c("x", "y"), crs = 4326)
  }





  # #
  # # 333 přiložení nových bodů z bias rasterů nově ohodnocených
  # #
  # print("lastChance")
  #
  # nw <- readRDS(paste0(path.tgob.tgbg, "i_",druh, ".rds"))
  #
  # #### musím normalizovat per pixel/druh/autor!!! Nemůžů počítat více nálezů téhož druhu pro autora vícekrát!!!
  # p.temp333base <- tgobV$p %>%  group_by(nc_cellNumber, nc_AUTOR, DRUH) %>% slice_head(n=1) %>% ungroup()
  #
  #
  # # jen autoři s nálezy
  # p.temp333a <- p.temp333base %>% filter(!is.na(cellNumber.all))
  #
  # # jen fs
  # p.temp333s <- p.temp333base %>%  filter(DRUH == druh)
  #
  # # ohodnocení autorů per obsazené pixely
  # p.temp333s.perPixelAutor <- p.temp333s %>%  group_by(nc_cellNumber) %>% mutate(perPixelAutor = n())  %>% ungroup() %>%  group_by(nc_AUTOR) %>% mutate(perAutorWeights = sum(1 / perPixelAutor)) %>% slice_head(n=1)
  #
  # p.temp333s.cnu <- unique(p.temp333s$nc_cellNumber)
  # #
  # p.temp333s.au <- unique(p.temp333s$nc_AUTOR)
  #
  #
  # collector.stack <- list()
  # # jdu po jednotlivých autorech
  # initStack <- TRUE
  # for(au in p.temp333s.au){
  #
  #
  #   au.stats  <- p.temp333s.perPixelAutor %>%  filter(nc_AUTOR == au)
  #
  #
  #   au.stats.nw  <-  nw[["weights"]] %>%  filter(nc_AUTOR == au)
  #
  #   print("au")
  #   print(au)
  #   print(au.stats.nw$nc_TGOB.sso.w333c)
  #
  #   adjusts <-  names(collector.aut[[au]][["br"]])
  #
  #   for(adjust in adjusts){
  #
  #     if(initStack){
  #       collector.stack[[adjust]] <- stack()
  #
  #     }
  #
  #
  #
  #     mm <- raster::setMinMax(collector.aut[[au]][["br"]][[adjust]])
  #
  #    if(!is.na(raster::minValue(mm)) & !is.na(raster::maxValue(mm)) & !is.na(au.stats.nw$nc_TGOB.sso.w333c) ) {
  #
  #
  #      collector.stack[[adjust]] <- stack(collector.stack[[adjust]], collector.aut[[au]][["br"]][[adjust]] * au.stats.nw$nc_TGOB.sso.w333c  ) # au.stats$perAutorWeights
  #
  #    }else{
  #
  #      print("vylučuju vylučuju vylučuju vylučuju vylučuju vylučuju")
  #    }
  #
  #   }
  #
  #   initStack <- FALSE
  # }
  # print("au1")
  #
  # collector.r.mean <- list()
  #
  # for(adjust in adjusts){
  #
  #   br.temp <- calc(collector.stack[[adjust]], mean)
  #   br.temp <- raster::setMinMax(br.temp)
  #
  #   # normalize to 0-1
  #   r.min <- raster::minValue(br.temp)
  #   r.max <- raster::maxValue(br.temp)
  #   br.temp <- ((br.temp - r.min) / (r.max - r.min))
  #   br.temp  <- raster::setMinMax(br.temp)
  #
  #   collector.r.mean[[adjust]] <- br.temp
  #   print("au2")
  #
  #   collector[["lastChance"]][["br"]][[adjust]]    <- collector.r.mean[[adjust]]
  #   collector[["lastChance"]][["bg"]][[adjust]]    <- tgbg::sample(collector.r.mean[[adjust]], n = 0.2, prob = TRUE)
  #
  # }
  #
  saveRDS(collector, paste0(path.tgob, "55ssos_collector.all_", druh, "_", rep, ".rds"))




  # příprava rasterů pro scénář: seskupení do jedné společné verze per adjust
  # ndop.fs$scenario.selected
  if (ndop.fs$scenario == "brAll") {
    collector.temp <- list()
    scenario.names <- names(ndop.fs$scenarios)
    for (scenario.name in scenario.names) {
      scenario <- unname(unlist(ndop.fs$scenarios[scenario.name]))
      collector.selected <- collector[scenario]

      collector.selected.r <- lapply(collector.selected, function(x) x$br)

      adjust.names <- names(collector.selected.r[[1]])
      version.names <- names(collector.selected.r)
      brs <- list()
      for (vn in version.names) {
        for (an in adjust.names) {
          brs[[an]][[vn]] <- collector.selected.r[[vn]][[an]]
        }
      }

      # počet verzí bias rasterů
      scenario.n <- length(scenario)

      brs.stack <- sapply(brs, function(x) {
        rs <- raster::stack(x)
        names(rs) <- names(x)
        return(rs)
      })

      brs.stack.vif <- list()
      if (scenario.n > 1) {
        for (ran in names(brs.stack)) {
          rs.temp <- brs.stack[[ran]]
          # vifstep
          vs <- vifstep(rs.temp, th = 3)
          # vifcor - jen pro dofiltr, asi není nutné
          vs.vc <- vifcor(rs.temp[[vs@results$Variables]], th = 0.7)
          rs.temp_vifs <- rs.temp[[vs.vc@results$Variables]]

          brs.stack.vif[[ran]] <- rs.temp_vifs
        }
      } else {
        brs.stack.vif <- brs.stack
      }

      collector.temp[[scenario.name]][["bg"]] <- as.list(sapply(adjust.names, function(x) NA))
      collector.temp[[scenario.name]][["br"]] <- brs.stack.vif
    }

    collector.temp[["un"]] <- collector[["un"]]

    collector <- collector.temp
  }











  # odstraním TO
  cnr <- names(collector)
  #  collector <-  collector[!str_detect(cnr, "TO")]
  # collector <-  collector[!str_detect(cnr, "w30")]
  # collector <-  collector[!str_detect(cnr, "w0")]

  # collector <-  collector[!str_detect(cnr, "w4b")]
  # collector <-  collector[!str_detect(cnr, "w8")]
  # collector <-  collector[!str_detect(cnr, "w9")]
  # collector <-  collector[!str_detect(cnr, "w5b")]
  # collector <-  collector[!str_detect(cnr, "w5")]
  # collector <-  collector[!str_detect(cnr, "w21")]
  # collector <-  collector[!str_detect(cnr, "w26")]
  # collector <-  collector[!str_detect(cnr, "w27")]
  # collector <-  collector[!str_detect(cnr, "w4c")]
  # collector <-  collector[!str_detect(cnr, "w40")]
  # collector <-  collector[!str_detect(cnr, "w24")]
  # collector <-  collector[!str_detect(cnr, "TS.sso")]
  # collector <-  collector[!str_detect(cnr, "TS.w.sso")]

  # příprava bg pro scénář: seskupení korelovaných rasterů do jedné společné verze generování bg per adjust
  # ndop.fs$scenario.selected
  if (ndop.fs$scenario == "bgAllneeeeeeeeeeeeeeeeeeeee") {
    collector.subset <- list()
    collector.rl <- sapply(collector, function(x) x$br)
    cn.names <- c()
    for (collector.name in names(collector.rl)) {
      for (collector.name.adjust in names(collector.rl[[collector.name]])) {
        cn.names <- c(cn.names, paste0(collector.name, "_", collector.name.adjust))
      }
    }


    collector.rs <- stack(unlist(collector.rl[1:(length(collector.rl) - 1)]))
    names(collector.rs) <- cn.names

    # předběžně uložit
    fnv <- paste0(path.tgob, "v_vs.vc_", druh, ".rds")
    if (file.exists(fnv)) {
      vs.vc <- readRDS(fnv)
    } else {
      if (rep > 1) {
        print("čekám 1 ------------------------------------------------")
        # 1
        Sys.sleep(10)

        # počkat na vytvoření souboru s vifcor, ať mám jeden jednotný
        if (file.exists(fnv)) {
          vs.vc <- readRDS(fnv)
        } else {
          # 2
          print("čekám 2 ------------------------------------------------")
          Sys.sleep(60)
          if (file.exists(fnv)) {
            vs.vc <- readRDS(fnv)
          } else {
            # 3
            print("čekám 3 ------------------------------------------------")
            Sys.sleep(60 * 5)
            vs.vc <- readRDS(fnv)
          }
        }
      } else {
        vs.vc <- usdmPB::vifcor(collector.rs, th = 0.95)
        saveRDS(vs.vc, fnv)
      }
    }



    select.bg <- vs.vc@results$Variables
    pairs <- vs.vc@excludedPairs
    add.singles <- vs.vc@results$Variables[vs.vc@results$Variables %notin% names(vs.vc@excludedPairs)]

    add.singles.lenght <- length(add.singles)
    append.singles <- list()
    if (add.singles.lenght > 0) {
      append.singles <- rep(NA, add.singles.lenght)
      names(append.singles) <- add.singles

      append.singles <- as.list(append.singles)
    }

    pairs.prep <- append(pairs, append.singles)

    pairs.prep <- pairs.prep[sort(names(pairs.prep))]

    pairs.prep.names <- c()
    for (pairs.prep.name in names(pairs.prep)) {
      pairs.prep.names <- c(pairs.prep.names, paste(na.omit(c(pairs.prep.name, pairs.prep[[pairs.prep.name]])), collapse = "|"))
    }

    pairs.prep.names <- sort(pairs.prep.names)

    cntr <- 1
    for (ss in names(pairs.prep)) {
      str.ss <- strsplit(ss, split = "_")[[1]]

      collector.subset[[str.ss[1]]][["bg"]][[str.ss[2]]] <- collector[[str.ss[1]]][["bg"]][[str.ss[2]]]
      collector.subset[[str.ss[1]]][["name"]][[str.ss[2]]] <- pairs.prep.names[cntr]
      cntr <- cntr + 1
    }
    collector.subset[["un"]] <- collector[["un"]]


    #  saveRDS(collector, paste0(path.tgob, "50ssos_collector.puv_", druh, "_", rep, ".rds"))




    collector <- collector.subset
  }






  # #
  # # 333 přiložení nových bodů z bias rasterů nově ohodnocených
  # #
  #
  #
  # print("lastChance")
  #
  #
  #
  # # neodstranil jsem zde presence z bg!!!!!!!
  #
  # # kernelIDs <- c(paste0("00",0:5), paste0("01",0:5), paste0("02",0:5), paste0("03",0:5),  paste0("10",0:5), paste0("11",0:5), paste0("12",0:5), paste0("13",0:5),   paste0("20",0:5), paste0("21",0:5), paste0("22",0:5), paste0("23",0:5)     )
  #
  # kernelIDs <- c(paste0("00",0:5), paste0("01",0:5), paste0("02",0:5), paste0("03",0:5),  paste0("10",0:5), paste0("11",0:5), paste0("12",0:5), paste0("13",0:5),   paste0("20",0:5), paste0("21",0:5), paste0("22",0:5), paste0("23",0:5)     )
  #
  #
  #
  # for(kdri in kernelIDs){
  #
  #
  #
  #
  # nw <- readRDS(paste0(path.tgob.tgbg, "kdr", kdri, "_",druh, ".rds"))
  # newBg <- list()
  # for(adj in as.character(ndop.fs$adjusts)){
  #
  #   br.temp <- raster::resample(raster::raster(nw[[adj]]), predictors[[1]], method = "bilinear")
  #   crs <- raster::crs(predictors[[1]])
  #   raster::crs(br.temp) <- crs
  #   br.temp <- raster::mask(raster::crop(br.temp, raster::extent(predictors[[1]])), predictors[[1]])
  #   br.temp[br.temp == -Inf] <- 0
  #   br.temp <- raster::setMinMax(br.temp)
  #
  #   # normalize to 0-1
  #   r.min <- raster::minValue(br.temp)
  #   r.max <- raster::maxValue(br.temp)
  #   br.temp <- ((br.temp - r.min) / (r.max - r.min))
  #   br.temp <- raster::setMinMax(br.temp)
  #
  #   # br.temp <- raster::mask(br.temp, as_Spatial(pres.unique))
  #
  #
  #   newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.2, prob = TRUE, crs = crs)
  #
  #   newBg[["name"]][[adj]] <- paste0("kernel", kdri)
  # }
  #
  #
  #
  # collector[[paste0("kernel", kdri)]] <- newBg
  #
  # }



  #  kernelIDs <- c("211", "999")
  # kernelIDs <- c(paste0("20",0:5), paste0("21",0:5), paste0("22",0:5), paste0("23",0:5)  , "990"  , "997" ,"998"  , "999" )



  tempD <- unlist(strsplit(list.files(path.tgob.tgbg, pattern = "^kdr"), "_"))
  kernelIDs <- unique(gsub("[kdr]", "", tempD[str_detect(tempD, "kdr")]))


  # kernelIDs <-  c("000000", "000001", "000010" ,"000050" ,"000100", "000200", "000500", "001000" ,"002000", "004000" ,"005000" ,"008000" ,"009000", "011000" ,"012000" ,"015000", "017000", "020000", "025000" ,"030000", "040000" ,"052000","100000" ,"300000" ,"990" ,   "997"   , "998"  ,  "999" )




  # kernelIDs <- c(paste0("22",0:5)  , paste0("23",0:5) , "999"   )
  # # kernelIDs <- 231


  # kernelIDs <- c(999, paste0("00",0:5), paste0("01",0:5), paste0("02",0:5), paste0("03",0:5),  paste0("10",0:5), paste0("11",0:5), paste0("12",0:5), paste0("13",0:5),   paste0("20",0:5), paste0("21",0:5), paste0("22",0:5), paste0("23",0:5)     )

  # kernelIDs <-  c("211")

  for (kdri in kernelIDs) {
    nw <- readRDS(paste0(path.tgob.tgbg, "kdr", kdri, "_", druh, ".rds"))
    newBg <- list()
    for (adj in as.character(ndop.fs$adjusts)) {
      br.temp <- raster::resample(raster::raster(sumNormal(nw[[adj]])), predictors.lsdMask[[1]], method = "bilinear") # sumNormal
      crs <- raster::crs(predictors.lsdMask[[1]])
      raster::crs(br.temp) <- crs
      br.temp <- raster::mask(raster::crop(br.temp, raster::extent(predictors.lsdMask[[1]])), predictors.lsdMask[[1]])
      # br.temp[br.temp == -Inf] <- 0
      # br.temp[br.temp < 0] <- 0
      br.temp <- raster::setMinMax(br.temp)




      #     xxx.t <-   terrain(br.temp, v="slope", neighbors=8, na.rm=TRUE)
      #     xxx.t4 <-   terrain(br.temp, opt="slope", neighbors=4)
      #     xxx.tr <-   terrain(br.temp, opt="TRI")
      #
      #     br.temp2 <-  br.temp ^ (1/10)
      #
      # # normalize to 0-1
      # r.min <- raster::minValue(br.temp2)
      # r.max <- raster::maxValue(br.temp2)
      # br.temp2 <- ((br.temp2 - r.min) / (r.max - r.min))
      # br.temp2 <- raster::setMinMax(br.temp2)
      #
      #     br.temp2 <- br.temp2/sum(as.vector(br.temp2), na.rm = TRUE)
      #
      #     xxx.t <-   terrain(br.temp, v="slope")
      #     xxx <- zonal(x=br.temp, fun=mean)





      br.temp.o <- br.temp
      br.temp <- br.temp / sum(as.vector(br.temp), na.rm = TRUE)
      br.temp <- raster::setMinMax(br.temp)



      r.min <- raster::minValue(br.temp)
      r.max <- raster::maxValue(br.temp)
      br.temp <- ((br.temp - r.min) / (r.max - r.min))
      br.temp <- raster::setMinMax(br.temp)


      #     br.temp.i <- raster::mask(br.temp, as_Spatial(pres.unique))
      #       v.tempX  <- na.omit( as.vector(br.temp.i) )
      #
      # bp.min <- min(max(v.tempX), quantile(v.tempX, prob=0.75) + 1.5 * IQR(v.tempX))
      # bp.max <- max(min(v.tempX), quantile(v.tempX, prob=0.25) + 1.5 * IQR(v.tempX))
      #
      # br.temp[br.temp >= bp.min] <- bp.min
      # br.temp[br.temp <= bp.max] <- bp.max
      #
      # br.temp <- raster::setMinMax(br.temp)

      newBg[["br"]][[adj]] <- br.temp

      # 0.05 0.1 ... 0.5
      if (str_detect(paste0("kdr", kdri), "kdr01")) {
        print(1)
        newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.02, prob = TRUE, crs = crs, replace = TRUE)
      }
      if (str_detect(paste0("kdr", kdri), "kdr02")) {
        print(2)
        newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.03, prob = TRUE, crs = crs, replace = TRUE)
      }
      if (str_detect(paste0("kdr", kdri), "kdr03")) {
        print(3)
        newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.05, prob = TRUE, crs = crs, replace = TRUE)
      }
      if (str_detect(paste0("kdr", kdri), "kdr04")) {
        print(4)
        newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.1, prob = TRUE, crs = crs, replace = TRUE)
      }
      if (str_detect(paste0("kdr", kdri), "kdr05")) {
        print(5)
        newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.2, prob = TRUE, crs = crs, replace = TRUE)
      }
      if (str_detect(paste0("kdr", kdri), "kdr06")) {
        print(6)
        newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.3, prob = TRUE, crs = crs, replace = TRUE)
      }

      # # změnu počtu backgroundu musím aplikovat i na TGOB!!! - musím mít více jeho verzí - to ale pomůže i TGOB...
      # if(!str_detect(paste0("kdr", kdri), "kdr0")){
      #   print(6)
      #   newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.2, prob = TRUE, crs = crs, replace=FALSE)
      # }




      # teoreticky bch pak z nasamplovaných bg bodů jejich součtem po průniku z pravděpodobností mohl sečís jestli se rovnají 1 nebo ne? Nebo i samplovat jen do 1 a pak to v tgbd:sample ukončit?
      # dál už nenormalizovat na 0-1?


      # br.temp <- br.temp/sum(as.vector(br.temp), na.rm = TRUE)
      # br.temp <- raster::setMinMax(br.temp)
      # newBg[["brOrig"]][[adj]] <- br.temp
      #
      # # # normalize to 0-1
      # # r.min <- raster::minValue(br.temp)
      # # r.max <- raster::maxValue(br.temp)
      # # br.temp <- ((br.temp - r.min) / (r.max - r.min))
      # # br.temp <- raster::setMinMax(br.temp)
      #
      #
      # br.temp.i  <- br.temp * -1
      # #  br.temp.i  <- br.temp^(2)/sum(as.vector(br.temp^(2)), na.rm=TRUE) * -1
      #
      # br.temp.i <- raster::mask(br.temp.i, as_Spatial(pres.unique), updatevalue=0.0)
      #
      # br.temp <- raster::mask(br.temp, as_Spatial(pres.unique), inverse=TRUE, updatevalue=0.0)
      #
      # newBg[["brOrig"]][[adj]] <- br.temp + raster::mask((newBg[["brOrig"]][[adj]]^(0.5) ) / sum(as.vector((newBg[["brOrig"]][[adj]]^(0.5) )),na.rm=TRUE), as_Spatial(pres.unique), updatevalue=0.0) * -1
      #
      # br.temp <- br.temp + br.temp.i
      #
      # # normalize to 0-1
      # r.min <- raster::minValue(br.temp)
      # r.max <- raster::maxValue(br.temp)
      # br.temp <- ((br.temp - r.min) / (r.max - r.min))
      # br.temp <- raster::setMinMax(br.temp)
      #
      # br.temp <- br.temp/sum(as.vector(br.temp), na.rm = TRUE)
      #
      #
      # br.temp <- raster::setMinMax(br.temp)
      # newBg[["br"]][[adj]] <- br.temp
      #
      # # normalize to 0-1 ---- jen pro bg sampling aby nebyla negativní pravděpodobnost,reálně stejně nepoužívám?
      # r.min <- raster::minValue(br.temp)
      # r.max <- raster::maxValue(br.temp)
      # br.temp <- ((br.temp - r.min) / (r.max - r.min))
      # br.temp <- raster::setMinMax(br.temp)





      # ####
      # #### pokus o stejný počet presencí podle intensity bg v jejich místě
      # ####
      #
      # newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.2, prob = TRUE, crs = crs)
      #
      #
      # testtt <- newBg[["bg"]][[adj]]
      #
      # testtt <-   st_as_sf( as_tibble(testtt) %>% mutate(bg = row_number()))
      #
      #
      # pres.unique.poly <-   st_as_sf( as_tibble(pres.unique.poly) %>% mutate(p = row_number()))
      #
      # res.unique.poly.centroid <-  pres.unique.poly %>% st_centroid()
      #
      # testtt.i <-  st_intersection( pres.unique.poly, testtt)
      #
      # res.unique.poly.centroid.notin <-   res.unique.poly.centroid  %>% filter(p %notin% unique(unlist(testtt.i$p) ))
      #
      # testtt.i.prep <- testtt.i %>% dplyr::select(-bg)
      #
      # testtt.i.prep  %<>% add_row(res.unique.poly.centroid.notin)
      #
      # newBg[["p"]][[adj]] <-  testtt.i.prep %>% dplyr::select(geometry)


      #
      # testtt.i.count <-   as_tibble(testtt.i) %>% group_by(p) %>% mutate(count = n()) %>% slice_head(n=1) %>% dplyr::select(-geometry)
      #
      # res.unique.poly.centroid.count <- res.unique.poly.centroid %>% left_join(testtt.i.count, by="p") %>% mutate(count = ifelse(is.na(count),1, count) )
      #
      # for(pid in res.unique.poly.centroid.count$p){
      #  tr <-  res.unique.poly.centroid.count %>% filter(p == pid)
      #  trc <- tr$count
      #  for(cid in 1:trc){
      #  if(cid > 1){
      # print(cid)
      #    res.unique.poly.centroid.count %<>% add_row(tr)
      #  }
      #  }
      # }
      # newBg[["p"]][[adj]] <- res.unique.poly.centroid.count  %>% dplyr::select(geometry)



      ####
      #### pokus o inverzní počet presencí podle intensity bg v jejich místě
      ####

      v.temp <- na.omit(raster::extract(br.temp.o, st_coordinates(pres.unique)))

      v.temp.bg <- na.omit(raster::extract(br.temp.o, st_coordinates(collector[["un"]][["bg"]][["0"]] %>% add_row(pres.unique) %>% distinct())))

      if (length(v.temp) != nrow(pres.unique)) {
        print("presence v NA rasterů!!!")
        stop()
      }

      # max outliery - dám na nové maximum
      bp.min <- min(max(v.temp), quantile(v.temp, prob = 0.75) + 1.5 * IQR(v.temp))
      bp.max <- max(min(v.temp), quantile(v.temp, prob = 0.25) + 1.5 * IQR(v.temp))

      v.temp[v.temp >= bp.min] <- bp.min
      v.temp[v.temp <= bp.max] <- bp.max


      # extrémní + a - hodnoty sjednotím - oboje je problematické.
      # + hodně lidí druh videlo na daném ístě - oversampling
      # - jeden člověk druh viděl, ale spousta ho tam nenašlo - podezřelý/náhodný/neopakovatelný/(těžce zopakovatelný)/(extrémní přesamplování jedním člověkem, časem najde všechno...) nález
      # v.temp <- abs(v.temp)



      bp.min.bg <- min(max(v.temp.bg), quantile(v.temp.bg, prob = 0.75) + 1.5 * IQR(v.temp.bg))
      v.temp.bg[v.temp.bg >= bp.min.bg] <- bp.min.bg



      # v.temp.puv <- v.temp
      # v.temp <- (v.temp )^(1/3)
      v.temp.reverse <- v.temp * -1 + max(v.temp) + min(v.temp)
      v.temp.bg.reverse <- v.temp.bg * -1 + max(v.temp.bg) + min(v.temp.bg)



      # v.temp.diff <-  max(v.temp) - min(v.temp)
      # # kolik jednotek bude rozdíl? - zkusit odhadnout ze sumy biasu jednotně pro všechny?
      # unit.diff <- 5
      #
      #
      #
      #  v.temp.reverse.inverse <-  ceiling((v.temp.reverse *  unit.diff / max(v.temp)))
      #  pres.unique.temp <-  as_tibble(pres.unique) %>% mutate(uid = row_number())
      #  pres.unique.temp$count <- v.temp.reverse.inverse
      #
      #  pres.unique.inverse <-   st_as_sf(pres.unique.temp )
      #
      #
      #  for(pid in  pres.unique.inverse$uid){
      #   tr <-  pres.unique.inverse %>% filter(uid == pid)
      #   trc <- tr$count
      #   for(cid in 1:trc){
      #   if(cid > 1){
      #     # print(cid)
      #     pres.unique.inverse %<>% add_row(tr)
      #   }
      #   }
      #  }
      #  newBg[["p"]][[adj]] <- pres.unique.inverse %>% dplyr::select(geometry)



      v.temp <- v.temp.reverse
      r.min <- min(v.temp)
      r.max <- max(v.temp)
      v.temp <- ((v.temp - r.min) / (r.max - r.min))


      v.temp.bg <- v.temp.bg.reverse
      r.min.bg <- min(v.temp.bg)
      r.max.bg <- max(v.temp.bg)
      v.temp.bg <- ((v.temp.bg - r.min.bg) / (r.max.bg - r.min.bg))

      # nechci nulové váhy - reálně to bude stejně jen jedna hodnota...

      v.temp[v.temp == 0] <- sort(unique(v.temp))[2] / 2
      newBg[["pW"]][[adj]] <- v.temp

      v.temp.bg[v.temp.bg == 0] <- sort(unique(v.temp.bg))[2] / 2
      newBg[["pW.bg"]][[adj]] <- (v.temp.bg + 0.5) * 100


      newBg[["name"]][[adj]] <- paste0("kernel", kdri)
    }

    # plot(`34ssos.new_Accipiter gentilis_1`[["Accipiter gentilis"]][["kernel233"]][["2"]][["LQ"]][["2"]])
    #  plot(pres.unique, add=TRUE, pch=4)
    #  plot(lsd.pa.centroids %>% filter(TaxonNameLAT == "Accipiter gentilis") %>% dplyr::select(presence), add=TRUE, col=as.factor(presence+1), pch=22)

    collector[[paste0("kernel", kdri)]] <- newBg
  }

  # t.ag %>% group_by(adjust, version) %>% mutate(AICc.sum = sum(AICc)) %>% mutate(AICc.mean = mean(AICc)) %>% mutate(AICc.sd = sd(AICc)) %>% slice_head(n=1) %>% dplyr::select(adjust, version, AICc.sum, AICc.mean, AICc.sd) %>% arrange(AICc.sum)

  #  t.ag %>% group_by(adjust, version) %>% mutate(AICc.sum = sum(cor.manual)) %>% mutate(AICc.mean = mean(cor.manual)) %>% mutate(AICc.sd = sd(cor.manual)) %>% slice_head(n=1) %>% dplyr::select(adjust, version, AICc.sum, AICc.mean, AICc.sd) %>% arrange(desc(AICc.mean))



  #
  #   newBg <- list()
  #   nw <- readRDS(paste0(path.tgob.b, "34ssos_",druh, "_1.rds"))
  #   for(adj in as.character(ndop.fs$adjusts)){
  #
  #
  #     br.temp <- nw[[druh]][["TGOB"]][[adj]][["L"]][["1"]]@predictions@layers[[1]]
  #       r.min <- raster::minValue(br.temp)
  #       r.max <- raster::maxValue(br.temp)
  #       br.temp <- ((br.temp - r.min) / (r.max - r.min))
  #       br.temp <- raster::setMinMax(br.temp)
  #
  #       br.temp <- raster::mask(br.temp, as_Spatial(pres.unique), inverse=TRUE, updatevalue=0.0)
  #
  #     newBg[["bg"]][[adj]] <- tgbg::sample(br.temp, 0.2, prob = TRUE, crs = crs)
  #
  #     newBg[["name"]][[adj]] <- paste0("kernel", 999)
  #   }
  #
  #   collector[[paste0("kernel", 999)]] <- newBg
  #




  #
  #   ###
  #   ### simulace randomizace
  #   ###
  #
  #   print("simulace randomizace")
  #
  #
  #   sampled.bg.collector <- list()
  #   for(i in 1:10){
  #     r.sampled <- sampleRandom(predictors[[1]], size=2000, cells=TRUE)
  #
  #
  #     sso.temp <- ndopP
  #
  #
  #     v.temp <- raster::extract(predictors[[1]], sf::st_coordinates(sso.temp), cellnumbers = TRUE)
  #     sso.temp$cellNumber3 <-  v.temp[, 1]
  #
  #     presences.selected <- unique( unlist(as_tibble(sso.temp) %>% ungroup()%>% filter(DRUH == druh )   %>% dplyr::select(-geometry)   %>% dplyr::select( cellNumber3)) )
  #     r.sampled.noPresences <- r.sampled[r.sampled[,1] %notin% presences.selected,]
  #
  #
  #     observers.selected <- unique(  unlist(as_tibble(sso.temp)  %>% dplyr::select(-geometry) %>% ungroup() %>% filter(cellNumber3 %in%  r.sampled[,1] ) %>% filter(DRUH == druh ) %>% dplyr::select(AUTOR)))
  #
  #     # presences.selected <- unique( unlist(as_tibble(sso.temp) %>% ungroup()%>% filter(DRUH == druh ) %>% filter(AUTOR %in%  observers.selected)   %>% dplyr::select(-geometry)   %>% dplyr::select( cellNumber3)) )
  #     # r.sampled.noPresences <- r.sampled[r.sampled[,1] %notin% presences.selected,]
  #
  #
  #
  #     # pokud prázdné, sampluju random!!!
  #
  #
  #     if(length(observers.selected) > 0 ){
  #       observers.selected.records <- sso.temp %>% ungroup() %>% filter(AUTOR %in%  observers.selected )
  #
  #
  #       sampled <- tgbg::bg(observers.selected.records %>% dplyr::select(geometry), predictors[[1]], sigma = ndop.fs$adjusts, anisotropic = TRUE)
  #
  #
  #       for(a in ndop.fs$adjusts){
  #
  #
  #         v.temp <- raster::extract(predictors[[1]], sf::st_coordinates(sampled$bg[[as.character(a)]]$geometry), cellnumbers = TRUE)
  #
  #         if(is.null(sampled.bg.collector$bg[[as.character(a)]])){
  #           sampled.bg.collector$bg[[as.character(a)]] <-c()
  #         }
  #
  #
  #
  #         # sampled.bg.collector$bg[[as.character(a)]]  <-   rbind(sampled.bg.collector$bg[[as.character(a)]], sampled$bg[[as.character(a)]][v.temp[, 1] %in% r.sampled.noPresences[,1],])
  #
  # sampled.bg.collector$bg[[as.character(a)]]  <-   rbind(sampled.bg.collector$bg[[as.character(a)]], sampled$bg[[as.character(a)]][v.temp[, 1] %notin% presences.selected,])
  #
  # sampled.bg.collector$name[[as.character(a)]]  <- paste0("rsampled_", as.character(a))
  #
  #
  #       }
  #
  #
  #
  #     }else{
  #       # random
  #       print("random, není presence xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
  #
  #
  #       sampled <- tgbg::sample(predictors[[1]])
  #
  #
  #       for(a in ndop.fs$adjusts){
  #
  #         if(is.null(sampled.bg.collector$bg[[as.character(a)]])){
  #           sampled.bg.collector$bg[[as.character(a)]] <-c()
  #         }
  #
  #         sampled.bg.collector$bg[[as.character(a)]]  <-   rbind(sampled.bg.collector$bg[[as.character(a)]], sampled$bg[[as.character(a)]])
  #
  #         sampled.bg.collector$name[[as.character(a)]]  <- paste0("rsampled_", as.character(a))
  #
  #
  #       }
  #
  #
  #
  #     }
  #
  #
  #
  #
  #
  #
  #
  #   }
  #
  # collector[["rsampled"]] <- sampled.bg.collector


  saveRDS(collector, paste0(path.tgob, "50ssos_collector.subset_", druh, "_", rep, ".rds"))







  gc()
  # run vsech variant BG s ENMeval

  for (id in c("TGOB01", "TGOB02", "TGOB03", "TGOB04", "TGOB05", "TGOB06", "un", paste0("kernel", kernelIDs))) { #  names(collector)  c("TGOB", "un", "kernel") , paste0("kernel",kernelIDs) c("TGOB", "un")  "kernel211", "kernel999"  , paste0("kernel",kernelIDs)
    # id.names <- unlist(strsplit(id, "_"))
    print("")
    print("")
    print("")
    print("***********************************************************************************************************************************")
    print(id)
    print("***********************************************************************************************************************************")



    # if ("un" == id) {
    # nms <-  names(collector[["kernel231"]][["bg"]])
    #
    # }else{
    #   nms <-  names(collector[[id]][["bg"]])
    # }
    #
    for (adjust in names(collector[[id]][["bg"]])) {
      if (str_detect(id, "kernel|TGOB")) {
        print("random bg ke kernelům a tgob !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§§§")
        print(id)
        # # # collector[[id]][["bg"]][[adjust]]  <- collector[["un"]][["bg"]][["0"]] %>% add_row(pres.unique)  %>% distinct() # # # xxx
        #   collector[[id]][["bg"]][[adjust]]  <- collector[["un"]][["bg"]][["0"]] %>% add_row(collector[[id]][["bg"]][[adjust]])
        print(nrow(collector[[id]][["bg"]][[adjust]] %>% add_row(pres.unique))) # %>% distinct())

        collector[[id]][["bg"]][[adjust]] <- collector[[id]][["bg"]][[adjust]] %>% add_row(pres.unique) #  %>% distinct() # # # xxx


        print(nrow(collector[[id]][["bg"]][[adjust]]))
        #  collector[[id]][["bg"]][[adjust]]  <- collector[["un"]][["bg"]][["0"]]

        collector[[id]][["allW"]][[adjust]] <- c(collector[[id]][["pW"]][[adjust]], rep(100, nrow(collector[[id]][["bg"]][[adjust]])))
        #  collector[[id]][["allW"]][[adjust]]  <- c(collector[[id]][["pW"]][[adjust]], collector[[id]][["pW.bg"]][[adjust]] )
      } else {
        collector[[id]][["pW"]][[adjust]] <- NA
        collector[[id]][["allW"]][[adjust]] <- NA
      }


      print(adjust)

      if (ndop.fs$scenario == "bg" | ndop.fs$scenario == "bgAll") {
        bg.temp <- as.data.frame(st_coordinates(collector[[id]][["bg"]][[adjust]]))
        block.bg <- collector[[id]][["bg"]][[adjust]] %>% st_join(bCV.poly)
        predictors.temp <- predictors
        print("původní prediktory")

        #### !!!!!!!!!!!! jen kvůli manipulaci presencí s random bg
        block.bg.null <- collector[["un"]][["bg"]][["0"]] %>% st_join(bCV.poly)

        # block.bg.null$fold
        predictors.temp.orig <- predictors
      }



      if (ndop.fs$scenario == "br" | ndop.fs$scenario == "brAll") {
        # dávám všude stený random, přidávám bias raster(y) k prediktorům
        bg.temp <- as.data.frame(st_coordinates(collector[["un"]][["bg"]][["0"]]))
        block.bg <- collector[["un"]][["bg"]][["0"]] %>% st_join(bCV.poly)
        block.bg.null <- collector[["un"]][["bg"]][["0"]] %>% st_join(bCV.poly)


        if ("un" == id) {
          # k random nepřidávám bias rastery
          predictors.temp <- predictors
        } else {
          predictors.temp <- stack(predictors, collector[[id]][["br"]][[adjust]])


          if (str_detect(id, "TGOB")) {
            predictors.temp.orig <- stack(predictors, collector[[id]][["br"]][[adjust]])
          } else {
            # vyhodnocení nad nezmanipulovaným bia rasterem
            predictors.temp.orig <- stack(predictors, collector[[id]][["brOrig"]][[adjust]])
          }

          print("přidány bias raster(y)")
          print(length(names(predictors.temp)))
          print(names(predictors.temp))
        }


        #         nw <- readRDS(paste0(path.tgob.b, "34ssos_",druh, "_1.rds"))
        #
        #         br.temp <- nw[[druh]][["TGOB"]][[adjust]][["L"]][["1"]]@predictions@layers[[1]]
        #         r.min <- raster::minValue(br.temp)
        #         r.max <- raster::maxValue(br.temp)
        #         br.temp <- ((br.temp - r.min) / (r.max - r.min))
        #         br.temp <- raster::setMinMax(br.temp)
        #
        # predictors.temp <- stack(predictors,br.temp)

        # nw <- readRDS(paste0(path.tgob.tgbg, "kdr999_",druh, ".rds"))
        # br.temp <- raster::resample(raster::raster(nw[[adjust]]), predictors[[1]], method = "bilinear")
        # crs <- raster::crs(predictors[[1]])
        # raster::crs(br.temp) <- crs
        # br.temp <- raster::mask(raster::crop(br.temp, raster::extent(predictors[[1]])), predictors[[1]])
        # br.temp[br.temp == -Inf] <- 0
        # br.temp <- raster::setMinMax(br.temp)
        # # normalize to 0-1
        # r.min <- raster::minValue(br.temp)
        # r.max <- raster::maxValue(br.temp)
        # br.temp <- ((br.temp - r.min) / (r.max - r.min))
        # br.temp <- raster::setMinMax(br.temp)
        # pstack <- stack(br.temp )
        #
        # br.temp.tgob <- br.temp
        #
        # nw <- readRDS(paste0(path.tgob.tgbg, "kdr211_",druh, ".rds"))
        # br.temp <- raster::resample(raster::raster(nw[[adjust]]), predictors[[1]], method = "bilinear")
        # crs <- raster::crs(predictors[[1]])
        # raster::crs(br.temp) <- crs
        # br.temp <- raster::mask(raster::crop(br.temp, raster::extent(predictors[[1]])), predictors[[1]])
        # br.temp[br.temp == -Inf] <- 0
        # br.temp <- raster::setMinMax(br.temp)
        # # normalize to 0-1
        # r.min <- raster::minValue(br.temp)
        # r.max <- raster::maxValue(br.temp)
        # br.temp <- ((br.temp - r.min) / (r.max - r.min))
        # br.temp <- raster::setMinMax(br.temp)
        # pstack <- stack(pstack, br.temp )
        # #
        # # # nw <- readRDS(paste0(path.tgob.tgbg, "kdr221_",druh, ".rds"))
        # # # br.temp <- raster::resample(raster::raster(nw[[adjust]]), predictors[[1]], method = "bilinear")
        # # # crs <- raster::crs(predictors[[1]])
        # # # raster::crs(br.temp) <- crs
        # # # br.temp <- raster::mask(raster::crop(br.temp, raster::extent(predictors[[1]])), predictors[[1]])
        # # # br.temp[br.temp == -Inf] <- 0
        # # # br.temp <- raster::setMinMax(br.temp)
        # # # # normalize to 0-1
        # # # r.min <- raster::minValue(br.temp)
        # # # r.max <- raster::maxValue(br.temp)
        # # # br.temp <- ((br.temp - r.min) / (r.max - r.min))
        # # # br.temp <- raster::setMinMax(br.temp)
        # # # pstack <- stack(pstack, br.temp )
        # # #
        # predictors.temp  <- pstack
      }

      version2 <- NA
      if (ndop.fs$scenario == "bgAll") {
        if (is.null(collector[[id]][["name"]][[adjust]])) {
          version2 <- id
        } else {
          version2 <- collector[[id]][["name"]][[adjust]]
        }
      }




      names(bg.temp) <- ll
      block.p <- pres.unique %>% st_join(bCV.poly)

      # if("TGOB" == id | "un" == id){
      #   block.p <- pres.unique %>% st_join(bCV.poly)
      #
      #   # if("un" == id){
      #   #   block.p <- pres.unique %>% st_join(bCV.poly)
      #   #
      #   #   df.temp <- as.data.frame(st_coordinates(collector[["kernel231"]][["p"]][[adjust]]))
      #   #   names(df.temp) <- ll
      #   # }
      # }else{
      #   print("měním presence - multi...")
      #   print(nrow(pres.unique))
      #   print(nrow(collector[[id]][["p"]][[adjust]]))
      #
      #
      #   ##nee - vstupuje do  evalDepFolds - tam, chci unikátní presence??? špatně bloky???  (pokud to neykoriguje sdm balíček)
      #   ## block.p <- collector[[id]][["p"]][[adjust]] %>% st_join(bCV.poly)
      #   block.p <- pres.unique %>% st_join(bCV.poly)
      #
      #   df.temp <- as.data.frame(st_coordinates(collector[[id]][["p"]][[adjust]]))
      #   names(df.temp) <- ll
      # }


      # # # xxx
      nw <- NA
      # if (str_detect(id, "kernel")) {
      #   print("NW do glmnet  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§§§")
      #   nw <-collector[[id]][["allW"]][[adjust]]
      #
      # }



      print("základní:")

      # FC
      for (FC in tune.args$fc) {
        # RM
        for (RM in tune.args$rm) {
          tune.args.ind <- list("fc" = FC, "rm" = RM)
          ee.temp <- NA
          if (inherits(try({
            ee.temp <- ENMevaluate(
              # user.grp = list("occs.grp" = block.p$fold, "bg.grp" = block.bg$fold),
              occs = df.temp,
              envs = predictors.temp,
              bg = bg.temp,
              algorithm = "maxnet", partitions = "none", # user
              # partition.settings = list("kfolds" = 3),
              tune.args = tune.args.ind,
              other.settings = list("addsamplestobackground" = TRUE, "newWeights" = nw, "other.args" = list("addsamplestobackground" = TRUE, "newWeights" = nw))
            )
            e.mx.all[[druh]][[id]][[adjust]][[FC]][[as.character(RM)]] <- ee.temp

            # při individuálních výpočtech možno takto, jinak je tam více layerů, ne jen jeden!!!
            e.mx.all.new[[druh]][[id]][[adjust]][[FC]][[as.character(RM)]] <- ee.temp@predictions[[1]]

            # e.mx.all.new[[druh]][[id]][[adjust]][[FC]][[as.character(RM)]] <- enm.maxnet@predict(ee.temp@models[[1]], predictors.temp.orig, list(doClamp = TRUE, pred.type = "cloglog"))# doClamp = FALSE - bez toho hází error?!?


            # export results + add indep test
            ei <- evalIndep(ee.temp, lsd.temp, predictors.temp.orig)

            # stop()
            #
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # print("*****************************")
            # stop()

            # random, dependent  evaluation

            depAUC.p <- pres.unique
            depAUC.p["presence"] <- 1

            depAUC.a <- collector[["un"]][["bg"]][["0"]]
            depAUC.a["presence"] <- 0


            depAUC.p %<>% add_row(depAUC.a)

            if ("un" == id) {
              ei2 <- evalDep(ee.temp, depAUC.p, predictors.temp.orig, NA, NA, NA)
              print("ei2 ei2 ei2")
              print(ei2)
              print(names(ei2))
            } else {

              # !!! zruseno
              # ei2 <- evalDep(ee.temp, depAUC.p, predictors.temp.orig, collector[[id]][["br"]][[adjust]], pres.unique, collector[[id]][["pW"]][[adjust]])
              ei2 <- evalDep(ee.temp, depAUC.p, predictors.temp.orig, NA, NA, NA)
              # random, dependent  evaluation
            }

            ee.temp@results[["species"]] <- druh
            ee.temp@results[["version"]] <- id
            ee.temp@results[["version2"]] <- version2
            ee.temp@results[["adjust"]] <- adjust

            ee.temp@results[["entropy"]] <- ee.temp@models[[1]][["entropy"]]


            ee.temp@results[["bg.n"]] <- nrow(bg.temp)



            # random, dependent  evaluation - folds ------------------------------


            ei3 <- evalDepFolds(ee.temp, depAUC.p, c(block.p$fold, block.bg.null$fold), predictors.temp.orig)


            if (length(ei) > 0 & length(ei2) > 0 & length(ei3) > 0) {
              temp.t <- ee.temp@results %>%
                left_join(ei, by = "tune.args") %>%
                left_join(ei2, by = "tune.args", suffix = c("", "_trainRandom")) %>%
                left_join(ei3, by = "tune.args", suffix = c("", "_trainRandomFold"))
            } else {
              print("evaluace indep daty nefunkční!!!")
              next
            }



            temp.t$occs.n <- nrow(ee.temp@occs)

            if (first) {
              first <- FALSE
              out.t <- temp.t
            } else {
              out.t %<>% add_row(temp.t)
            }
          }), "try-error")) {
            e.mx.all[[druh]][[id]][[adjust]][[FC]][[as.character(RM)]] <- NA
          }

          # varianta s thinningem presencí přidaných do UN
          if ("un" == id) {
            # next  # zatím přeskakuji
            print("thin:")
            for (thinDist in names(collector.thin[[id]])) {
              df.temp.thin <- as.data.frame(st_coordinates(collector.thin[[id]][[thinDist]]))
              names(df.temp.thin) <- ll
              # podstrčit thinning presence místo původních
              print("měním presenční dataset pro thinnovací verze")

              block.pt <- collector.thin[[id]][[thinDist]] %>% st_join(bCV.poly)
              ee.temp <- NA
              if (inherits(try({
                ee.temp <- ENMevaluate(
                  # user.grp = list("occs.grp" = block.pt$fold, "bg.grp" = block.bg$fold),
                  occs = df.temp.thin,
                  envs = predictors.temp,
                  bg = bg.temp,
                  algorithm = "maxnet", partitions = "none", # user
                  # partition.settings = list("kfolds" = 3),
                  tune.args = tune.args.ind,
                  other.settings = list("addsamplestobackground" = TRUE, "other.args" = list("addsamplestobackground" = TRUE))
                )
                e.mx.all[[druh]][[id]][[thinDist]][[FC]][[as.character(RM)]] <- ee.temp

                e.mx.all.new[[druh]][[id]][[thinDist]][[FC]][[as.character(RM)]] <- enm.maxnet@predict(ee.temp@models[[1]], predictors.temp.orig, list(doClamp = FALSE, pred.type = "cloglog"))


                # export results + add indep test (thin)
                ei <- evalIndep(ee.temp, lsd.temp, predictors.temp.orig)



                # random, dependent  evaluation - tady asi zbytečné, ale musím zachovat strukturu atributů

                depAUC.p <- pres.unique
                depAUC.p["presence"] <- 1

                depAUC.a <- collector[["un"]][["bg"]][["0"]]
                depAUC.a["presence"] <- 0


                depAUC.p %<>% add_row(depAUC.a)

                ei2 <- evalDep(ee.temp, depAUC.p, predictors.temp.orig, NA, NA)
                # random, dependent  evaluation





                ee.temp@results[["species"]] <- druh
                ee.temp@results[["version"]] <- id
                ee.temp@results[["version2"]] <- id
                ee.temp@results[["adjust"]] <- thinDist

                ee.temp@results[["bg.n"]] <- nrow(bg.temp)



                # random, dependent  evaluation - folds ------------------------------

                ei3 <- evalDepFolds(ee.temp, depAUC.p, c(block.p$fold, block.bg.null$fold), predictors.temp.orig)



                if (length(ei) > 0 & length(ei2) > 0 & length(ei3) > 0) {
                  temp.t <- ee.temp@results %>%
                    left_join(ei, by = "tune.args") %>%
                    left_join(ei2, by = "tune.args", suffix = c("", "_trainRandom")) %>%
                    left_join(ei3, by = "tune.args", suffix = c("", "_trainRandomFold"))
                } else {
                  print("evaluace indep daty nefunkční!!!")
                  next
                }


                temp.t$occs.n <- nrow(ee.temp@occs)

                if (first) {
                  first <- FALSE
                  out.t <- temp.t
                } else {
                  out.t %<>% add_row(temp.t)
                }
              }), "try-error")) {
                e.mx.all[[druh]][[id]][[thinDist]][[FC]][[as.character(RM)]] <- NA
              }
            }
          }
        }
      }
    }
  }
  ### saveRDS(e.mx.all, paste0(path.tgob, "34ssos_", druh, "_", rep, ".rds"))
  saveRDS(e.mx.all.new, paste0(path.tgob, "34ssos.new_", druh, "_", rep, ".rds"))
  saveRDS(out.t, paste0(path.tgob, "t_90___ssos_", druh, "_", rep, ".rds"))
  gc()
  # } # for DRUH
}
end_time <- Sys.time()
print(end_time - start_time)

print(sp.group)