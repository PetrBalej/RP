cmd_arg <- commandArgs(trailingOnly = TRUE)
if (is.na(cmd_arg[1])) {
    print("*********************************************************************************************************************")
    print("nezadán parametr pro část druhů")
    print("*********************************************************************************************************************")
    cmd_arg <- 1
} else {
    cmd_arg <- cmd_arg[1]
}

print(cmd_arg)
# "C:\Program Files\R\R-4.2.1\bin\x64\Rscript.exe" "D:\PersonalWork\Balej\v2\RP\RP\R\ndopTGOB.R" 1
# library(devtools)
# install_github("PetrBalej/tgbg", force = TRUE)
# install_github("PetrBalej/tgbg", force = TRUE, ref="dev", build=TRUE)
# install_local("/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/tgbg/tgbg", force = TRUE, build=TRUE)
# vlastní nutné modifikace existujících balíčků:
# install_github("PetrBalej/sdm", force = TRUE, ref="advancedTresholdMetrics", build=TRUE)
# install_github("PetrBalej/usdm", force = TRUE, ref="pb", build=TRUE)
# install_github("PetrBalej/ENMeval/", force = TRUE, build=TRUE) # přímo v masteru (+ je vhodné ponechávat původní název? - různé plusy i minusy, už jsem řešil...)
# devtools::document("/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/tgbg/tgbg")

# *** Running ENMeval v2.0.4 with maxnet from maxnet package v0.1.4 *** -- na WS1
# "C:/Users/svc.kagup/AppData/Local/Programs/R/R-4.2.1/bin/x64/Rscript.exe" "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/R/ndopTGOB.R" 1
start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "ENMeval", "ecospat", "blockCV", "tgbg", "sdmATM", "usdmPB") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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
path.tgob <- paste0(path.prep, "ndopTGOB/") # path.wd.prep.ndop



############
# inputs
############

grid <- readRDS(paste0(path.prep, "sitmap_2rad/sitmap_2rad-czechia.rds")) # 4326 sitmap_2rad.czechia
predictors <- stack(readRDS(paste0(path.prep, "predictors/20192022/preds_final_vifs.rds"))) # vifstep2 pcamap6

# NDOP dofiltrované pro SSOS
# NDOP druhy presence
ndopP <- readRDS(paste0(path.ndop, "ndopP.rds")) %>% dplyr::select(POLE, DRUH, AUTOR)
ndopP.POLE <- readRDS(paste0(path.ndop, "ndopP.POLE.rds")) %>% dplyr::select(POLE, DRUH, AUTOR)

# LSD PA
lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds")) %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci
lsd.pa.min <- readRDS(paste0(path.lsd, "lsd.pa.min.rds"))

# statistiky NDOP druhy
ndop.stat.res.sp <- readRDS(paste0(path.ndop, "ndop.stat.res.sp.rds"))

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

# plocha ČR
czechia.sq <- nrow(grid)

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
    "1.sso.w.dnorm1.mean" = "TGOB.sso.w.dnorm1.mean",
    "1.sso.w.3.dnorm1.mean" = "TGOB.sso.w.3.dnorm1.mean",
    "1.sso.w.dnorm1.median" = "TGOB.sso.w.dnorm1.median",
    "1.sso.w.3.dnorm1.median" = "TGOB.sso.w.3.dnorm1.median",
    #
    "1.sso.w.dnorm2.mean" = "TGOB.sso.w.dnorm2.mean",
    "1.sso.w.3.dnorm2.mean" = "TGOB.sso.w.3.dnorm2.mean",
    "1.sso.w.dnorm2.median" = "TGOB.sso.w.dnorm2.median",
    "1.sso.w.3.dnorm2.median" = "TGOB.sso.w.3.dnorm2.median"
)


ndop.fs <- list(
    "adjusts" = c(0.1, 1, 2, 3),
    "tuneArgs" = list(fc = c("L", "LQ", "LQH"), rm = c(0.5, 1, 2, 3, 4)), # , "H", "LQH", "LQHP"
    "bgRatio" = 1 / 10,
    "bg" = 1000,
    "speciesPerGroup" = 2, "speciesOccMin" = 30,
    "sq2rad" = c((1 / 6) / 4, 0.1 / 4), # kvadráty KFME 2rad, xy velikost ve stupních
    "sq2radDist" = c(1:5),
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

"%notin%" <- Negate("%in%")


evalIndep <- function(m, lsd.temp) {
    first <- TRUE
    out.t <- NA
    for (layer in names(m@predictions)) { # v případě rozdělení per RM+FC je cyklus zbytečný (ale preventivně nechat), je tam jen jeden RasterLayer
        ev <- NA
        r.temp <- m@predictions[[layer]]
        ex.predicted <- extract(r.temp, st_coordinates(lsd.temp))
        #
        # dopočíst a přidat další metriky z performance() (rgee)!!!
        # tam ale není nic co zohledňuje tn tp (bez poměru k fp)
        #

        if (inherits(try({
            ev <- sdmATM::evaluates(lsd.temp$presence, ex.predicted)
        }), "try-error")) {
            ev <- NA
        }
        if (is.na(ev)) {
            next
        }
        ev.temp <- as.data.frame(ev@statistics[-3])
        ev.temp <- merge(ev.temp, t(as.data.frame(ev@statistics$COR)))
        ev.temp <- merge(ev.temp, ev@threshold_based[1, ]) # sp=se
        ev.temp[["tune.args"]] <- layer
        if (first) {
            first <- FALSE
            out.t <- ev.temp
        } else {
            out.t %<>% add_row(ev.temp)
        }
    }

    return(out.t)
}


# kešovat
fn <- paste0(path.tgob, "tgobV.rds")
if (file.exists(fn)) {
    print("načítám existující tgobV.rds")
    # tgobV <- readRDS(fn)

    # # dočasně, nepoužívám AO - odstranit, ať to nežere místo
    # tgbg.names  <-  names(tgobV$p)
    # tgobV$p %<>% dplyr::select(tgbg.names[!str_detect(tgbg.names, "\\.AO")])
    # tgbg.names  <-  names(tgobV$p)
    # tgobV$p %<>% dplyr::select(tgbg.names[!str_detect(tgbg.names, "sso_")])
    # saveRDS(tgobV, paste0(path.tgob, "tgobV2.rds"))
    tgobV <- readRDS("D:/PersonalWork/Balej/v2/RP/dataPrep/ndopTGOB/tgobV2.rds")
} else {
    print("generuju nový tgobV.rds")
    tgobV <- tgbg::tgobVersions(ndopP, predictors[[1]], species = "DRUH", observers = "AUTOR", species.select = lsd.pa.min$TaxonNameLAT)
    saveRDS(tgobV, fn)
}

prefix <- "nc_"
cn <- names(tgobV$p)
cn_prefix <- cn[str_detect(cn, paste0("^", prefix))]

ndop.stat.res.sp.selected <- ndop.stat.res.sp %>%
    filter(DRUH %in% unique(lsd.pa.min$TaxonNameLAT)) %>%
    arrange(POLE.n)
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
        que <- ndop.stat.res.sp.selected %>% arrange(POLE.n)
        que$rep <- 0
        que$lock <- 0
        saveRDS(que, paste0(path.tgob, que.fn))
    }


    locked <- nrow(que %>% filter(lock == 1))
    que.r <- que %>%
        filter(rep < replicates) %>%
        arrange(POLE.n) %>%
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
                arrange(POLE.n) %>%
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
    null.default <- sample(predictors[[1]])

    #  for (druh in sp.group$DRUH) { # for DRUH
    first <- TRUE # pro každou replikaci a druh samostatný soubor
    collector <- list()
    collector.thin <- list()
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(druh)
    print(sp.group)
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    e.mx.all <- list()
    out.t <- list()
    bg.col <- list()
    # LSD
    lsd.temp <- NA
    lsd.temp <- lsd.pa.centroids %>% filter(TaxonNameLAT == druh)
    #
    # základní verze backgroundů
    #

    # TGOB ssos
    ssos.sp <- cn_prefix[str_detect(cn_prefix, paste0(druh, "$"))]

    versions.base.puv <- c("TGOB", "TO", "TS", "TO.w", "TS.w")
    versions.base <- c("TGOB", "TO.w", "TS.w")
    ao.sp <- ssos.sp[!str_detect(ssos.sp, "sso_")]
    ao.sp.puv <- gsub(prefix, "", ao.sp)
    ao.sp <- c()
    ssos.w <- ssos.sp[str_detect(ssos.sp, "sso\\.w")]
    ssos.w <- gsub(prefix, "", ssos.w)

    ssos.sp <- ssos.sp[str_detect(ssos.sp, "sso_")]
    # ssos.sp <- gsub(prefix, "", ssos.sp)

    for (vn in c(versions.base, ao.sp, ssos.w)) {
        print("----------------------------------------------------------------------------------")
        print(vn)
        vns <- strsplit(vn, "_")[[1]][1] # odstraním případnou příponu s druhem
        print(vns)
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
                filter(DRUH == druh) %>%
                dplyr::select(geometry)
        } else {
            if (str_detect(vn, "\\.w")) {
                # beru všechny Top, ale musím odstranit NA, bez vah
                p.temp <- tgobV$p %>%
                    # filter(DRUH != druh) %>%
                    filter(!is.na(!!sym(paste0(prefix, vn)))) %>%
                    dplyr::select(geometry, !!sym(paste0(prefix, vn)))
            } else {
                p.temp <- tgobV$p %>%
                    # filter(DRUH != druh) %>%
                    filter(!!sym(paste0(prefix, vn)) == 1) %>%
                    dplyr::select(geometry)
            }
        }

        if (nrow(p.temp) < 1) {
            print("neexistují presence základní verze!!!")
            next
        }

        if (str_detect(vn, "\\.w")) {
            collector[[vns]] <- tgbg::bg(p.temp %>% dplyr::select(geometry), predictors[[1]], sigma = ndop.fs$adjusts, weights = p.temp[[paste0(prefix, vn)]], output = ndop.fs$outputBgBr, anisotropic = TRUE)
        } else {
            collector[[vns]] <- tgbg::bg(p.temp, predictors[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)

            if (vn == "TGOB") {
                # zkresluju podle druhu
                collector[[paste0(vns, ".bss")]] <- tgbg::bg(p.temp.bss, predictors[[1]], sigma = ndop.fs$adjusts, output = ndop.fs$outputBgBr, anisotropic = TRUE)
            }
        }

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


    # thinning  presencí pro UN
    for (thinDist in sq2rad.dist.range) {
        occs.thinned <- ecospat.occ.desaggregation(xy = tgob.trad.tm.df, min.dist = thinDist, by = NULL)
        thinDist <- as.character(thinDist)
        collector.thin[[id]][[thinDist]] <- occs.thinned %>% st_as_sf(coords = c("x", "y"), crs = 4326)
    }

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


    # příprava bg pro scénář: seskupení korelovaných rasterů do jedné společné verze generování bg per adjust
    # ndop.fs$scenario.selected
    if (ndop.fs$scenario == "bgAll") {
        collector.subset <- list()
        collector.rl <- sapply(collector, function(x) x$br)
        cn.names <- c()
        for (collector.name in names(collector.rl)) {
            for (collector.name.adjust in names(collector.rl[[collector.name]])) {
                cn.names <- c(cn.names, paste0(collector.name, "_", collector.name.adjust))
            }
        }

        collector.rs <- stack(collector.rl)
        names(collector.rs) <- cn.names

        vs.vc <- usdmPB::vifcor(collector.rs, th = 0.7)

        # předběžně uložit
        saveRDS(vs.vc, paste0(path.tgob, "v_vs.vc_", druh, "_", rep, ".rds"))

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
        for (ss in select.bg) {
            str.ss <- strsplit(ss, split = "_")[[1]]

            collector.subset[[str.ss[1]]][["bg"]][[str.ss[2]]] <- collector[[str.ss[1]]][["bg"]][[str.ss[2]]]
            collector.subset[[str.ss[1]]][["name"]][[str.ss[2]]] <- pairs.prep.names[cntr]
            cntr <- cntr + 1
        }
        collector.subset[["un"]] <- collector[["un"]]

        collector <- collector.subset
    }

    gc()
    # run vsech variant BG s ENMeval
    for (id in names(collector)) {
        # id.names <- unlist(strsplit(id, "_"))

        print(id)
        for (adjust in names(collector[[id]][["bg"]])) {
            print(adjust)

            if (ndop.fs$scenario == "bg" | ndop.fs$scenario == "bgAll") {
                bg.temp <- as.data.frame(st_coordinates(collector[[id]][["bg"]][[adjust]]))
                block.bg <- collector[[id]][["bg"]][[adjust]] %>% st_join(bCV.poly)
                predictors.temp <- predictors
                print("původní prediktory")
            }
            if (ndop.fs$scenario == "br" | ndop.fs$scenario == "brAll") {
                # dávám všude stený random, přidávám bias raster(y) k prediktorům
                bg.temp <- as.data.frame(st_coordinates(collector[["un"]][["bg"]][["0"]]))
                block.bg <- collector[["un"]][["bg"]][["0"]] %>% st_join(bCV.poly)
                if ("un" == id) {
                    # k random nepřidávám bias rastery
                    predictors.temp <- predictors
                } else {
                    predictors.temp <- stack(predictors, collector[[id]][["br"]][[adjust]])
                    print("přidány bias raster(y)")
                    print(length(names(predictors.temp)))
                    print(names(predictors.temp))
                }
            }

            version2 <- NA
            if (ndop.fs$scenario == "bgAll") {
                version2 <- collector[[id]][["name"]][[adjust]]
            }

            names(bg.temp) <- ll

            block.p <- pres.unique %>% st_join(bCV.poly)

            print("základní:")

            # FC
            for (FC in tune.args$fc) {
                # RM
                for (RM in tune.args$rm) {
                    tune.args.ind <- list("fc" = FC, "rm" = RM)
                    ee.temp <- NA
                    if (inherits(try({
                        ee.temp <- ENMevaluate(
                            user.grp = list("occs.grp" = block.p$fold, "bg.grp" = block.bg$fold),
                            occs = df.temp,
                            envs = predictors.temp,
                            bg = bg.temp,
                            algorithm = "maxnet", partitions = "user",
                            # partition.settings = list("kfolds" = 3),
                            tune.args = tune.args.ind,
                            other.settings = list("addsamplestobackground" = FALSE, "other.args" = list("addsamplestobackground" = FALSE))
                        )
                        e.mx.all[[druh]][[id]][[adjust]][[FC]][[as.character(RM)]] <- ee.temp
                        # export results + add indep test
                        ei <- evalIndep(ee.temp, lsd.temp)

                        ee.temp@results[["species"]] <- druh
                        ee.temp@results[["version"]] <- id
                        ee.temp@results[["version2"]] <- version2
                        ee.temp@results[["adjust"]] <- adjust

                        temp.t <- ee.temp@results %>% left_join(ei, by = "tune.args")

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
                                    user.grp = list("occs.grp" = block.pt$fold, "bg.grp" = block.bg$fold),
                                    occs = df.temp.thin,
                                    envs = predictors.temp,
                                    bg = bg.temp,
                                    algorithm = "maxnet", partitions = "user",
                                    # partition.settings = list("kfolds" = 3),
                                    tune.args = tune.args.ind,
                                    other.settings = list("addsamplestobackground" = FALSE, "other.args" = list("addsamplestobackground" = FALSE))
                                )
                                e.mx.all[[druh]][[id]][[thinDist]][[FC]][[as.character(RM)]] <- ee.temp

                                # export results + add indep test (thin)
                                ei <- evalIndep(ee.temp, lsd.temp)

                                ee.temp@results[["species"]] <- druh
                                ee.temp@results[["version"]] <- id
                                ee.temp@results[["version2"]] <- id
                                ee.temp@results[["adjust"]] <- thinDist

                                temp.t <- ee.temp@results %>% left_join(ei, by = "tune.args")

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

    saveRDS(e.mx.all, paste0(path.tgob, "10ssos_", druh, "_", rep, ".rds"))
    saveRDS(out.t, paste0(path.tgob, "t_10ssos_", druh, "_", rep, ".rds"))
    gc()
    # } # for DRUH
}
end_time <- Sys.time()
print(end_time - start_time)

print(sp.group)