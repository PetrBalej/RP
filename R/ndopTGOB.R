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
# install_github("PetrBalej/ENMeval", force = TRUE)
# *** Running ENMeval v2.0.4 with maxnet from maxnet package v0.1.4 *** -- na WS1
# "C:/Users/svc.kagup/AppData/Local/Programs/R/R-4.2.1/bin/x64/Rscript.exe" "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/R/ndopTGOB.R" 1
start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "ENMeval", "ecospat", "blockCV") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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

ndop.topAUTOR100 <- readRDS(paste0(path.ndop, "ndop.topAUTOR100.rds")) %>% dplyr::select(POLE, DRUH, AUTOR)
ndop.topAUTOR10 <- readRDS(paste0(path.ndop, "ndop.topAUTOR10.rds")) %>% dplyr::select(POLE, DRUH, AUTOR)

ndop.topDRUH50 <- readRDS(paste0(path.ndop, "ndop.topDRUH50.rds")) %>% dplyr::select(POLE, DRUH, AUTOR)
ndop.topDRUH10 <- readRDS(paste0(path.ndop, "ndop.topDRUH10.rds")) %>% dplyr::select(POLE, DRUH, AUTOR)


# LSD PA
# lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds"))  %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci
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

ndop.fs <- list(
    "adjusts" = c(0.1, 0.5, 1, 2, 3, 4),
    "tuneArgs" = list(fc = c("L", "LQ", "LQH"), rm = c(0.5, 1, 2, 3, 4)), # , "H", "LQH", "LQHP"
    "bgRatio" = 1 / 10,
    "speciesPerGroup" = 2, "speciesOccMin" = 30,
    "sq2rad" = c((1 / 6) / 4, 0.1 / 4), # kvadráty KFME 2rad, xy velikost ve stupních
    "sq2radDist" = c(1:5),
    "replicates" = 3,
    "speciesPart" = cmd_arg, "version" = "v1",
    "bgRaster" = FALSE,
    "versionNames" = vn, "versionSmooting" = vf
)
# dopočty
ndop.fs$bg <- round(cellStats(!is.na(predictors[[1]]), sum) * ndop.fs$bgRatio)
# !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!!
ndop.fs$bgRatio <- ndop.fs$bg # dočasně přepisuju - normálně dělat obě varianty!!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!! !!!!


############
# execution
############

# odstranit LSD sites z rastru prediktorů - jen pro generování tgob/ssoob - možná až na finální, teď ne
# odstraním prezence (sites) LSD z NDOP pro zabráníení spatial (auto)correlation - LSD tak bude i prostorově nezávislý dataset

"%notin%" <- Negate("%in%")


ndop.stat.res.sp.selected <- ndop.stat.res.sp %>%
    filter(DRUH %in% unique(lsd.pa.min$TaxonNameLAT)) %>%
    arrange(POLE.n)
sp.diff <- setdiff(unique(lsd.pa.min$TaxonNameLAT), ndop.stat.res.sp.selected$DRUH)
if (length(sp.diff) > 0) {
    print("Problém s nejednotnou taxonomií nebo neprůnikem LSD druhů s NDOP druhy!")
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


generateRPall <- function(RasterLayer_OR_sf_POINT, nBackRatio = 1 / 10, prob = FALSE) {
    # počítám s tím, že dostanu už unikátní body per pixel!!!
    # pak můžu vždy poměrově vybrat část k "ploše"

    # https://github.com/danlwarren/ENMTools/blob/004a4a1e182127900a5f62bc015770479bcd0415/R/check.bg.R#L138-L144
    # raster::sampleRandom ani sample(.., replace=FALSE) nejsou vhodné/nefungují pro menší N (clipnuté rastery skillem)
    out <- list()
    if (is(RasterLayer_OR_sf_POINT, "RasterLayer")) {
        background.points <- as.data.frame(rasterToPoints(RasterLayer_OR_sf_POINT))
    } else if (is(RasterLayer_OR_sf_POINT, "sf") && st_geometry_type(RasterLayer_OR_sf_POINT, by_geometry = FALSE) == "POINT" && prob == FALSE) {
        background.points <- as.data.frame(st_coordinates(RasterLayer_OR_sf_POINT))
        colnames(background.points) <- c("x", "y")
    } else {
        print("Vstup generateRPall musí být RasterLayer nebo sf-POINT!")
        print("sf-POINT pouze s random=TRUE")
        return(out)
    }

    if (nBackRatio > 1) {
        # neudávám poměr ale přesný počet
        total <- nrow(background.points)
        if (nBackRatio > total) {
            # nemůžu nasamplovat více
            proportion <- total
        } else {
            proportion <- nBackRatio
        }
    } else {
        # udávám poměr z rasteru nebo bodů
        total <- nrow(background.points)
        proportion <- round(total * nBackRatio)
    }



    probs <- NULL
    if (prob == TRUE) {
        probs <- background.points[, 3]
    }
    inds <- sample(1:total,
        size = proportion,
        prob = probs,
        replace = FALSE
    )
    bg.temp <- background.points[inds, 1:2]
    return(bg.temp %>% st_as_sf(coords = c("x", "y"), crs = 4326))
}

smoothingRP <- function(raster.template, adjusts, occs.sf, nBackRatio = 1 / 10, bgRaster = FALSE) {
    if (!is(raster.template, "RasterLayer")) {
        print("Vstup raster.template musí být RasterLayer!")
        return(list())
    }

    rcrs <- crs(raster.template)
    ext <- extent(raster.template)
    ow <- owin(xrange = c(ext@xmin, ext@xmax), yrange = c(ext@ymin, ext@ymax))

    res.ndop.coords <- st_coordinates(occs.sf)
    res.ndop.coords.ppp <- ppp(res.ndop.coords[, 1], res.ndop.coords[, 2], window = ow)

    # raster_stack.bias.col <- list()
    bg <- list()
    bg.raster <- list()

    # spíše do nastavení?
    sq2rad <- c((1 / 6) / 4, 0.1 / 4)
    adjusts.df <- as.data.frame(cbind("adj_i" = adjusts, "x" = sq2rad[1] * adjusts, "y" = sq2rad[2] * adjusts))

    for (adj in adjusts) {
        fa <- as.vector(unlist(adjusts.df %>% filter(adj_i == adj) %>% dplyr::select(x, y)))
        raster_stack.bias <- resample(raster(density.ppp(res.ndop.coords.ppp, sigma = fa)), raster.template, method = "bilinear")
        crs(raster_stack.bias) <- rcrs
        raster_stack.bias <- mask(crop(raster_stack.bias, extent(raster.template)), raster.template)
        raster_stack.bias <- setMinMax(raster_stack.bias)

        # normalizace
        r.min <- minValue(raster_stack.bias)
        r.max <- maxValue(raster_stack.bias)
        raster_stack.bias <- ((raster_stack.bias - r.min) / (r.max - r.min))
        raster_stack.bias <- setMinMax(raster_stack.bias)

        bg.raster[[as.character(adj)]] <- NA
        if (bgRaster == TRUE) {
            bg.raster[[as.character(adj)]] <- raster_stack.bias
        }

        bg[[as.character(adj)]] <- generateRPall(raster_stack.bias, nBackRatio = nBackRatio, prob = TRUE)
    }

    return(list("bg" = bg, "bg.raster" = bg.raster))
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



# společný základ BG pro všechny druhy
bgSources <- list(
    "tgob" = ndopP,
    # "ssos" = tgob,
    # "ssos.2" = ssos,
    # "topA.10" = ndop.topAUTOR10, # pokrývají jen část ČR - nerovnoměrně, přestože mají značnou část nálezů
    "topA100" = ndop.topAUTOR100,
    "topS10" = ndop.topDRUH10
    # "topS.50" = ndop.topDRUH50
)

# fix area
bgSources.fix <- list()
for (bgSourcesName in names(bgSources)) {
    temp.POLE <- as_tibble(bgSources[[bgSourcesName]]) %>%
        dplyr::select(-geometry) %>%
        group_by(POLE) %>%
        slice_head(n = 1) %>%
        dplyr::select(POLE)
    bgSources.fix[[bgSourcesName]] <- tgob.trad.POLE %>%
        ungroup() %>%
        filter(POLE %in% unique(unname(unlist(temp.POLE)))) %>%
        dplyr::select(-everything())
}




for (rep in 1:ndop.fs$replicates) {
    print("RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR")
    print(rep)
    set.seed(rep)

    # bloky
    bCV <- blockCV::cv_spatial(st_as_sf(rasterToPoints(predictors[[1]], spatial = TRUE)) %>% dplyr::select(-everything()), size = 50000, deg_to_metre = 90000, k = 5, selection = "random", hexagon = TRUE, seed = rep, plot = FALSE)
    bCV.poly <- as_tibble(bCV[["blocks"]][["geometry"]])
    bCV.poly$fold <- bCV[["blocks"]][["folds"]]
    bCV.poly <- st_as_sf(bCV.poly)

    # výchozí random BG společný pro všechny
    null.default <- generateRPall(predictors[[1]], nBackRatio = ndop.fs$bgRatio)

    first <- TRUE
    for (druh in sp.group$DRUH) { #  as.vector(sp.group$DRUH) speciesParts[[ndop.fs$speciesPart]] c("Turdus viscivorus")
        collector <- list()
        collector.thin <- list()
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print(druh)
        print(sp.group)
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        e.mx.all <- list()
        bg.col <- list()

        pres <- ndopP %>% filter(DRUH == druh)

        pres.au <- unique(as.vector(unlist(pres$AUTOR)))
        # "skill TGOB" - vyberu pro TGOB pouze nálezy těch autorů, kteří už daný druh pozorovali
        # Ostatní jsou nejistí, respektive jejich sampling effort je irelevantní, protože není jisté, zda-li:
        # - druh vůbec poznají,
        # - jsou ochotní ho zaznamenávat nebo
        # - vůbec nemapují v oblasti výskytu druhu

        tgob <- ndopP %>% filter(AUTOR %in% pres.au)

        # SSOS - nad top50 druhy

        # tgob.top50 <- bgSources$topS.50 %>% filter(AUTOR %in% pres.au)
        tgob.top10 <- bgSources$topS10 %>% filter(AUTOR %in% pres.au)
        tgob.top100 <- bgSources$topA100 %>% filter(AUTOR %in% pres.au)


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

        # buffer kolem presencí
        # náhodně může vznikat více typů geometrií, pak to spadne do sfc_GEOMETRY, to nechci (neumí s tím pak pracovat dále některé funkce), vybírám pouze (MULTI)POLYGON
        # buffer je 4*strana 2rad, od středu kvadrátu tedy cca 12 km, reálně mi jde o to získat buffer cca 3 čtverců okolo sousedících
        pres.unique.buffer <- st_as_sf(st_union(st_buffer(pres.unique, 10000))) %>% filter(st_geometry_type(x) %in% c("MULTIPOLYGON", "POLYGON"))
        predictors.ssos.buffer <- mask(predictors[[1]], pres.unique.buffer)


        tgob.trad.tm <- tgob.trad.sp %>% filter(DRUH == druh)
        tgob.trad.tm.df <- as.data.frame(st_coordinates(tgob.trad.tm))
        colnames(tgob.trad.tm.df) <- c("x", "y")

        #
        # ssos - omezení autorů s "nepřirozeným"" poměrem sběru daného druhu
        #

        ### celkové
        # ssos.DRUH <- as_tibble(tgob) %>%
        #     dplyr::select(-geometry) %>%
        #     group_by(DRUH) %>%
        #     summarise(POLE.n = n_distinct(POLE))
        # ssos.DRUH.sum.all <- sum(ssos.DRUH$POLE.n)
        # ssos.DRUH.sum <- as.numeric(ssos.DRUH %>% filter(DRUH == druh) %>% dplyr::select(POLE.n))
        # # "přirozený" poměr
        # ssos.DRUH.ratio <- (ssos.DRUH.sum / ssos.DRUH.sum.all) * 0.1 # snížím hranici o 50%, odstraní extrémy (lidi co druh značí jen občas)

        ### autorské
        ssos.DRUH.AUTOR <- as_tibble(tgob) %>%
            dplyr::select(-geometry) %>%
            group_by(DRUH, AUTOR) %>%
            summarise(POLE.n = n_distinct(POLE))
        ssos.DRUH.AUTOR.sum.all <- ssos.DRUH.AUTOR %>%
            ungroup() %>%
            group_by(AUTOR) %>%
            summarise(POLE.n.sum.all = sum(POLE.n))
        ssos.DRUH.AUTOR.sum <- ssos.DRUH.AUTOR %>%
            ungroup() %>%
            filter(DRUH == druh) %>%
            group_by(AUTOR) %>%
            summarise(POLE.n.sum = sum(POLE.n))

        ssos.DRUH.AUTOR.ratio <- ssos.DRUH.AUTOR.sum %>%
            left_join(ssos.DRUH.AUTOR.sum.all, by = "AUTOR") %>%
            mutate(ratio = POLE.n.sum / POLE.n.sum.all)
        # %>%  mutate(ratioLess = ifelse(ratio < ssos.DRUH.ratio, 1, 0))

        ssos.DRUH.ratio.quantile <- unname(quantile(ssos.DRUH.AUTOR.ratio$ratio, probs = c(0.01, 0.10, 0.25))) #- asi nechci automaticky odstranit 1/4, jen lidi s nižším než 50% z SSOS poměru zastoupení druhu oproti ostatním druhům

        ssos2_001 <- ssos.DRUH.AUTOR.ratio %>% filter(ratio > ssos.DRUH.ratio.quantile[1])
        pres.au <- unique(ssos2_001$AUTOR)
        ssos2_001.P <- ndopP %>% filter(AUTOR %in% pres.au)
        # ssos2_001.P.top50 <- bgSources$topS.50 %>% filter(AUTOR %in% pres.au)
        ssos2_001.P.top10 <- bgSources$topS10 %>% filter(AUTOR %in% pres.au)
        ssos2_001.P.top100 <- bgSources$topA100 %>% filter(AUTOR %in% pres.au)

        ssos2_010 <- ssos.DRUH.AUTOR.ratio %>% filter(ratio > ssos.DRUH.ratio.quantile[2])
        pres.au <- unique(ssos2_010$AUTOR)
        ssos2_010.P <- ndopP %>% filter(AUTOR %in% pres.au)
        # ssos2_010.P.top50 <- bgSources$topS.50 %>% filter(AUTOR %in% pres.au)
        ssos2_010.P.top10 <- bgSources$topS10 %>% filter(AUTOR %in% pres.au)
        ssos2_010.P.top100 <- bgSources$topA100 %>% filter(AUTOR %in% pres.au)

        ssos.DRUH.AUTOR.ratio.puv <- ssos.DRUH.AUTOR.ratio


        bgSources.ssos.temp1 <- list(
            "ssos" = tgob,
            # "ssos.0.topS50" = tgob.top50,
            "ssos.topS10" = tgob.top10,
            "ssos.topA100" = tgob.top100,
            "ssos2" = ssos2_001.P,
            # "ssos.1.010" = ssos2_010.P,
            # "ssos.1.025" = ssos2_025.P,
            # "ssos.1.001.topS50" = ssos2_001.P.top50,
            # "ssos.1.010.topS50" = ssos2_010.P.top50,
            "ssos2.topS10" = ssos2_001.P.top10,
            # "ssos.1.010.topS10" = ssos2_010.P.top10,
            "ssos2.topA100" = ssos2_001.P.top100
            # "ssos.1.010.topA100" = ssos2_010.P.top100
            # "ssos.1.025.topS50" = ssos2_025.P.top
        )


        # ### autorské bez jedničkových
        # ### autorské

        # ssos.DRUH.AUTOR.ratio %<>% filter(POLE.n.sum > 1)

        # ssos.DRUH.ratio.quantile <- unname(quantile(ssos.DRUH.AUTOR.ratio$ratio, probs = c(0.01, 0.10, 0.25))) #- asi nechci automaticky odstranit 1/4, jen lidi s nižším než 50% z SSOS poměru zastoupení druhu oproti ostatním druhům

        # ssos2_001 <- ssos.DRUH.AUTOR.ratio %>% filter(ratio > ssos.DRUH.ratio.quantile[1])
        # pres.au <- unique(ssos2_001$AUTOR)
        # ssos2_001.P <- ndopP %>% filter(AUTOR %in% pres.au)
        # ssos2_001.P.top50 <- bgSources$topS.50 %>% filter(AUTOR %in% pres.au)
        # ssos2_001.P.top10 <- bgSources$topS.10 %>% filter(AUTOR %in% pres.au)

        # ssos2_010 <- ssos.DRUH.AUTOR.ratio %>% filter(ratio > ssos.DRUH.ratio.quantile[2])
        # pres.au <- unique(ssos2_010$AUTOR)
        # ssos2_010.P <- ndopP %>% filter(AUTOR %in% pres.au)
        # ssos2_010.P.top50 <- bgSources$topS.50 %>% filter(AUTOR %in% pres.au)
        # ssos2_010.P.top10 <- bgSources$topS.10 %>% filter(AUTOR %in% pres.au)


        # bgSources.ssos.temp2 <- list(
        #     "ssos.2.001" = ssos2_001.P,
        #     "ssos.2.010" = ssos2_010.P,
        #     # "ssos.1.025" = ssos2_025.P,
        #     "ssos.2.001.topS50" = ssos2_001.P.top50,
        #     "ssos.2.010.topS50" = ssos2_010.P.top50,
        #     "ssos.2.001.topS10" = ssos2_001.P.top10,
        #     "ssos.2.010.topS10" = ssos2_010.P.top10
        #     # "ssos.1.025.topS50" = ssos2_025.P.top
        # )

        # bgSources.ssos <- append(bgSources.ssos.temp1, bgSources.ssos.temp2)

        bgSources.ssos <- bgSources.ssos.temp1


        # ze všech fixy
        bgSources.ssos.fix <- list()

        for (bgSourcesName in names(bgSources.ssos)) {
            temp.POLE <- as_tibble(bgSources.ssos[[bgSourcesName]]) %>%
                dplyr::select(-geometry) %>%
                group_by(POLE) %>%
                slice_head(n = 1) %>%
                dplyr::select(POLE)

            bgSources.ssos.fix[[bgSourcesName]] <- tgob.trad.POLE %>%
                ungroup() %>%
                filter(POLE %in% unique(unname(unlist(temp.POLE)))) %>%
                dplyr::select(-everything())
        }
        # spojím se statickými

        bgSources.all <- append(bgSources, bgSources.ssos)
        bgSources.all.fix <- append(bgSources.fix, bgSources.ssos.fix)


        # # rastery z fixů
        # predictors.v <- list()
        # for (bgSourcesName in names(bgSources.all.fix)) {
        #     predictors.v[[bgSourcesName]] <- mask(predictors[[1]], bgSources.all.fix[[bgSourcesName]])
        # }


        # kss <- predictors.v
        kss <- list()
        kss$cz <- predictors[[1]]

        # kss <- list(
        #   "cz" = predictors[[1]]
        #   # ,"buffer" = predictors.ssos.buffer
        # )


        print("start sampling variant BG")
        #
        # KSS
        #

        for (kssName in names(kss)) {
            if (str_detect(kssName, "top")) {
                next
            }

            print("raster kss: --------")
            print(kssName)
            for (bgSourcesName in names(bgSources.all)) {
                # if (str_detect(kssName, "ssos.0") & (str_detect(bgSourcesName, "ssos.1") | str_detect(bgSourcesName, "ssos.2"))) {
                #     next
                # }
                # if (str_detect(kssName, "ssos.1") & (str_detect(bgSourcesName, "ssos.0") | str_detect(bgSourcesName, "ssos.2"))) {
                #     next
                # }
                # if (str_detect(kssName, "ssos.2") & (str_detect(bgSourcesName, "ssos.0") | str_detect(bgSourcesName, "ssos.1"))) {
                #     next
                # }

                print("bg source:")
                print(bgSourcesName)
                # id <- paste(c(kssName, bgSourcesName), collapse = "_")
                id <- bgSourcesName
                collector.temp <- smoothingRP(kss[[kssName]], ndop.fs$adjusts, bgSources.all[[bgSourcesName]], nBackRatio = ndop.fs$bgRatio, bgRaster = ndop.fs$bgRaster)
                collector[[id]] <- collector.temp
            }
        }
        # # # buffer
        # id <- paste(c("kss", "buffer"), collapse = "_")
        # collector.temp <- smoothingRP(kss[["buffer"]], c(0.1, 1, 4), pres, nBackRatio = ndop.fs$bgRatio, bgRaster = ndop.fs$bgRaster) # vizuální kontrola!!!
        # collector[[id]] <- collector.temp

        # #
        # # RSS - ne z ssos (může být hodně málo... - stačí area_ssos.x)
        # #

        # for (bgSourcesName in names(bgSources.all.fix)) {
        #     if (bgSourcesName %in% c("ssos", "ssos.2", "ssos.top")) {
        #         next
        #     }
        #     id <- paste(c("rss", bgSourcesName), collapse = "_")
        #     collector.temp <- list()
        #     collector.temp[["bg"]][["0"]] <- generateRPall(bgSources.all.fix[[bgSourcesName]], nBackRatio = ndop.fs$bgRatio)
        #     collector[[id]] <- collector.temp
        # }
        # # # buffer
        # id <- paste(c("rss", "buffer"), collapse = "_")
        # collector.temp <- list()
        # collector.temp[["bg"]][["0"]] <- generateRPall(kss[["buffer"]], nBackRatio = ndop.fs$bgRatio) # vizuální kontrola!!!
        # collector[[id]] <- collector.temp


        #
        # area
        #

        # for (bgSourcesName in names(bgSources.all.fix)) {
        #     id <- paste(c("area", bgSourcesName), collapse = "_")
        #     collector.temp <- list()
        #     collector.temp[["bg"]][["0"]] <- bgSources.all.fix[[bgSourcesName]]
        #     collector[[id]] <- collector.temp
        # }

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


        # ze všech fixů (skoro) udělat masky pro clip
        #
        # přidání clip verze tam kde dává smysl
        #

        # ### čím clipuju
        # for (idc in names(collector)) {
        #     id.namesc <- unlist(strsplit(idc, "_"))
        #
        #     # nemá smysl clipovat jiné než celorepublikové varianty (ostatní už jsou samy o sobě nějakým způsobem clipnuté)
        #     # neclipuju area+ssos.x
        #     if ("area" == id.namesc[1] & "buffer" != id.namesc[2] & "tgob" != id.namesc[2] & !str_detect(id.namesc[2], "top") & is.na(id.namesc[3])) {
        #         for (adjustc in names(collector[[idc]][["bg"]])) {
        #             print("clip čím *****************************")
        #             print(idc)
        #             ssos.mask <- mask(predictors[[1]], collector[[idc]][["bg"]][[adjustc]])
        #
        #             ### co clipuju
        #             for (id in names(collector)) {
        #                 id.names <- unlist(strsplit(id, "_"))
        #                 if ((id.namesc[2] == id.names[2]) | (str_detect(id.namesc[2], "ssos") & str_detect(id.names[2], "ssos"))) {
        #                     # neclipuju sebe sama a ssos vzájemně
        #                     next
        #                 }
        #
        #                 if ("kss" == id.names[1] & "buffer" != id.names[2] & is.na(id.names[3])) {
        #                     for (adjust in names(collector[[id]][["bg"]])) {
        #                         print("clip čeho ******")
        #                         print(id)
        #                         ex.temp <- extract(ssos.mask, st_coordinates(collector[[id]][["bg"]][[adjust]]))
        #                         collector[[paste0(id, "_", id.namesc[2])]][["bg"]][[adjust]] <- collector[[id]][["bg"]][[adjust]][!is.na(ex.temp), ]
        #                     }
        #                 }
        #             }
        #         }
        #     }
        # }

        gc()
        # run vsech variant BG s ENMeval
        for (id in names(collector)) {
            # id.names <- unlist(strsplit(id, "_"))

            print(id)
            for (adjust in names(collector[[id]][["bg"]])) {
                print(adjust)

                bg.temp <- as.data.frame(st_coordinates(collector[[id]][["bg"]][[adjust]]))
                names(bg.temp) <- ll

                block.bg <- collector[[id]][["bg"]][[adjust]] %>% st_join(bCV.poly)
                block.p <- pres.unique %>% st_join(bCV.poly)

                print("základní:")

                # FC
                for (FC in tune.args$fc) {
                    # RM
                    for (RM in tune.args$rm) {
                        tune.args.ind <- list("fc" = FC, "rm" = RM)

                        if (inherits(try({
                            e.mx.all[[druh]][[id]][[adjust]][[FC]][[as.character(RM)]] <- ENMevaluate(
                                user.grp = list("occs.grp" = block.p$fold, "bg.grp" = block.bg$fold),
                                occs = df.temp,
                                envs = predictors,
                                bg = bg.temp,
                                algorithm = "maxnet", partitions = "user",
                                # partition.settings = list("kfolds" = 3),
                                tune.args = tune.args.ind,
                                other.settings = list("addsamplestobackground" = FALSE, "other.args" = list("addsamplestobackground" = FALSE))
                            )
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

                                if (inherits(try({
                                    e.mx.all[[druh]][[id]][[thinDist]][[FC]][[as.character(RM)]] <- ENMevaluate(
                                        user.grp = list("occs.grp" = block.pt$fold, "bg.grp" = block.bg$fold),
                                        occs = df.temp.thin,
                                        envs = predictors,
                                        bg = bg.temp,
                                        algorithm = "maxnet", partitions = "user",
                                        # partition.settings = list("kfolds" = 3),
                                        tune.args = tune.args.ind,
                                        other.settings = list("addsamplestobackground" = FALSE, "other.args" = list("addsamplestobackground" = FALSE))
                                    )
                                }), "try-error")) {
                                    e.mx.all[[druh]][[id]][[thinDist]][[FC]][[as.character(RM)]] <- NA
                                }
                            }
                        }
                    }
                }
            }
        }

        saveRDS(e.mx.all, paste0(path.tgob, "3ssos_", druh, "_", ndop.fs$speciesPart, "_", rep, ".rds"))
        gc()
    }
}
end_time <- Sys.time()
print(end_time - start_time)

print(sp.group)