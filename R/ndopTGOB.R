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

library(tidyverse)
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

path.wd <- paste0(gcfl(), "/../")


# *** Running ENMeval v2.0.4 with maxnet from maxnet package v0.1.4 *** -- na WS1
# "C:/Users/svc.kagup/AppData/Local/Programs/R/R-4.2.1/bin/x64/Rscript.exe" "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/R/ndopTGOB.R" 1
start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "geosphere", "ENMeval", "ecospat") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- paste0(path.wd, "/../projects-data/") # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.rgee <- paste0(path.wd, "/../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "dataPrep/ndopTGOB/")
path.wd.prep.ndop <- paste0(path.wd, "dataPrep/ndop/")
source(paste0(path.wd, "shared.R"))

############
# inputs
############

sitmap_2rad.czechia <- readRDS(paste0(path.wd, "dataPrep/sitmap_2rad/sitmap_2rad-czechia.rds")) # 4326

pcamap6 <- readRDS(paste0(path.wd, "dataPrep/observerSkill/pcamap6.rds"))
vifstep2 <- stack(paste0(path.wd, "dataPrep/observerSkill/vifstep2.grd")) # .rds varianta nefunkční
pcamap6 <- vifstep2 # chci původní rastery, použít PCA je blbost, poslední osy postupně dávají nesmysly
# souřadnice na SF - původní počty pro PPP
ndop.ff.au.more.sf <- readRDS(paste0(path.wd.prep.ndop, "ndop.ff.au.more.sf.rds"))
ndop.ff.au.more.sf.POLE <- readRDS(paste0(path.wd.prep.ndop, "ndop.ff.au.more.sf.POLE.rds"))
ndop.ff.au.more.sf %<>% arrange(ID_NALEZ)
ndop.ff.au.more.sf.POLE %<>% arrange(ID_NALEZ)
# sum(ndop.ff.au.more.sf$ID_NALEZ == ndop.ff.au.more.sf.POLE$ID_NALEZ)
ndop.ff.au.more.sf$POLE <- ndop.ff.au.more.sf.POLE$POLE


# per pixel - unikátní hodnoty druhů na pixel a autora
ndop.ff.au.more.sf.POLE.pp <- readRDS(paste0(path.wd.prep.ndop, "ndop.ff.au.more.sf.POLE.pp.rds"))
# jen body (centroidy pixelů)
ndop.ff.au.more.sf.POLE.pp.c <- readRDS(paste0(path.wd.prep.ndop, "ndop.ff.au.more.sf.POLE.pp.c.rds"))

# vybrané druhy (musím znovu přepočítat po odebrání LSD)
ndop.ff.au.more.sf.POLE.pp.sp <- readRDS(paste0(path.wd.prep.ndop, "ndop.ff.au.more.sf.POLE.pp.sp.rds"))

### LSD pro filtraci druhů které lze ověřit
lsd.pa.min <- readRDS(paste0(path.wd, "dataPrep/lsd/lsd.pa.min.rds"))
lsd.POLE <- st_read(paste0(path.wd, "dataPrep/lsd/overview-selected-2rad.shp"))

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

ndop.fs <- list(
    "adjusts" = c(0.1, 0.5, 1, 2, 3, 4),
    "tuneArgs" = list(fc = c("L", "LQ"), rm = c(1, 2, 3, 5, 10)),
    "bg" = 5000, "speciesPerGroup" = 7, "speciesOccMin" = 30,
    "sq2rad" = c((1 / 6) / 4, 0.1 / 4), # kvadráty KFME 2rad, xy velikost ve stupních
    "sq2radDist" = c(1:5),
    "speciesPart" = cmd_arg, "version" = "v1",
    "versionNames" = vn, "versionSmooting" = vf
)

############
# execution
############

# odstranit LSD sites z rastru prediktorů - jen pro generování tgob/ssoob
lsd.POLE.c <- st_centroid(lsd.POLE)
pcamap6.lsdClip <- stack(mask(pcamap6, lsd.POLE.c, inverse = TRUE))

"%notin%" <- Negate("%in%")
# odstraním prezence (sites) LSD z NDOP pro zabráníení spatial (auto)correlation - LSD tak bude i prostorově nezávislý dataset
ndop.ff.au.more.sf.POLE.pp.c %<>% ungroup() %>% filter(POLE %notin% lsd.POLE$POLE)


# musím znovu spočítat počty presencí po odebrání LSD sites
ndop.sp.selected <- ndop.ff.au.more.sf.POLE.pp.sp %>%
    ungroup() %>%
    filter(POLE %notin% lsd.POLE$POLE) %>%
    st_drop_geometry() %>%
    group_by(DRUH) %>%
    count(DRUH) %>%
    arrange(desc(n))

ndop.sp.selected %<>% filter(n >= ndop.fs$speciesOccMin) %>% filter(DRUH %in% unique(ndop.ff.au.more.sf.POLE.pp.c$DRUH))

### vybrat pouze druhy z LSD
# sjednocení synonym LSD
lsd.pa.min %<>% rename(species = TaxonNameLAT)
lsd.pa.min.syn <- synonyms_unite(lsd.pa.min)
lsd.pa.min.syn.species <- unique(lsd.pa.min.syn$species)
# alespoň 10 presencí a 10 absencí - doplnit do lsd.R a tam vygenerovat

# sjednocení synonym NDOP
ndop.sp.selected %<>% rename(species = DRUH)
ndop.sp.selected <- synonyms_unite(ndop.sp.selected)

ndop.sp.selected %<>% filter(species %in% lsd.pa.min.syn.species)
ndop.sp.selected %<>% rename(DRUH = species)

# druhy s malým počtem presencí se počítají mnohem rychleji, chci tyto per group upřednostnit (seřadit rovnoměrně v rámci skupin druhy od nejméně početných)
ndop.sp.selected %<>% arrange(n)
ndop.sp.selected %<>% ungroup()
ndop.sp.selected %<>% mutate(nthGroup = NA)
rowsTotal <- nrow(ndop.sp.selected)
groups <- ceiling(rowsTotal / ndop.fs$speciesPerGroup)

first <- TRUE
for (nth in 1:groups) {
    fill <- seq(nth, rowsTotal, groups)
    temp.rows <- ndop.sp.selected %>% filter(row_number() %in% fill)
    temp.rows %<>% mutate(nthGroup = nth)
    if (first) {
        first <- FALSE
        ndop.sp.selected.regrouped <- temp.rows
    } else {
        ndop.sp.selected.regrouped %<>% add_row(temp.rows)
    }
}

# vyřezané testovací LSD z presencí NDOP-u
ex.temp <- extract(pcamap6.lsdClip[[1]], st_coordinates(ndop.ff.au.more.sf))
ndop.ff.au.more.sf.extract <- ndop.ff.au.more.sf %>% filter(ID_NALEZ %in% ndop.ff.au.more.sf$ID_NALEZ[!is.na(ex.temp)])

generateRPall <- function(RasterLayer_OR_sf_POINT, nback = 5000, random = TRUE) {
    # https://github.com/danlwarren/ENMTools/blob/004a4a1e182127900a5f62bc015770479bcd0415/R/check.bg.R#L138-L144
    # raster::sampleRandom ani sample(.., replace=FALSE) nejsou vhodné/nefungují pro menší N (clipnuté rastery skillem)
    out <- list()
    if (is(RasterLayer_OR_sf_POINT, "RasterLayer")) {
        background.points <- as.data.frame(rasterToPoints(RasterLayer_OR_sf_POINT))
    } else if (is(RasterLayer_OR_sf_POINT, "sf") && st_geometry_type(RasterLayer_OR_sf_POINT, by_geometry = FALSE) == "POINT" && random == TRUE) {
        background.points <- as.data.frame(st_coordinates(RasterLayer_OR_sf_POINT))
        colnames(background.points) <- c("x", "y")
    } else {
        print("Vstup generateRPall musí být RasterLayer nebo sf-POINT!")
        print("sf-POINT pouze s random=TRUE")
        return(out)
    }

    full <- FALSE
    if (nrow(background.points) < nback) {
        full <- TRUE
    }

    if (random == FALSE) {
        # a1)
        inds <- sample(1:nrow(background.points),
            size = nback,
            prob = background.points[, 3],
            replace = TRUE
        )
        bg.temp <- background.points[inds, 1:2]
        out["dupl"] <- bg.temp %>% st_as_sf(coords = c("x", "y"), crs = 4326)
        out["dupld"] <- bg.temp %>%
            distinct() %>%
            st_as_sf(coords = c("x", "y"), crs = 4326)
        if (full) {
            bg.temp <- background.points[, 1:2]
        } else {
            inds <- sample(1:nrow(background.points),
                size = nback,
                prob = background.points[, 3],
                replace = FALSE
            )
            bg.temp <- background.points[inds, 1:2]
        }

        # a2)
        out["uniq"] <- bg.temp %>%
            st_as_sf(coords = c("x", "y"), crs = 4326)
    }

    if (random == TRUE) {
        # b1)
        inds <- sample(1:nrow(background.points),
            size = nback,
            prob = NULL,
            replace = TRUE
        )
        bg.temp <- background.points[inds, 1:2]
        out["dupl"] <- bg.temp %>% st_as_sf(coords = c("x", "y"), crs = 4326)
        out["dupld"] <- bg.temp %>%
            distinct() %>%
            st_as_sf(coords = c("x", "y"), crs = 4326)
        if (full) {
            bg.temp <- background.points[, 1:2]
        } else {
            inds <- sample(1:nrow(background.points),
                size = nback,
                prob = NULL,
                replace = FALSE
            )
            bg.temp <- background.points[inds, 1:2]
        }
        # b2)
        out["uniq"] <- bg.temp %>%
            st_as_sf(coords = c("x", "y"), crs = 4326)
    }
    return(out)
}

smoothingRP <- function(raster.template, adjusts, occs.sf, nback = 5000) {
    rcrs <- crs(raster.template)
    ext <- extent(raster.template)
    ow <- owin(xrange = c(ext@xmin, ext@xmax), yrange = c(ext@ymin, ext@ymax))

    res.ndop.coords <- st_coordinates(occs.sf)
    res.ndop.coords.ppp <- ppp(res.ndop.coords[, 1], res.ndop.coords[, 2], window = ow)

    raster_stack.bias.col <- list()
    bg <- list()
    ppp.deg.distance <- bw.scott(res.ndop.coords.ppp)
    prague <- c(14.4, 50.0)
    # vzdálenosti v jednotkovém ppp na jedné lat a lon
    d.lon <- distHaversine(c(prague[1], prague[2] - ppp.deg.distance[2]), c(prague[1], prague[2]))
    d.lat <- distHaversine(c(prague[1] - ppp.deg.distance[1], prague[2]), c(prague[1], prague[2]))
    # 0.001 a 0.01 jsou 100% korelované

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
        bg[[as.character(adj)]] <- generateRPall(raster_stack.bias, nback, random = FALSE)
    }

    # return(list(
    #     # zdrojové rastery neukládám, zbytečně velké...
    #     "bg" = bg, "bw" = ppp.deg.distance, "d.lon" = d.lon, "d.lat" = d.lat
    # ))

    return(bg)
}



#
# run per species models and varying BG
#
e.mx.all <- list()
tune.args <- ndop.fs$tuneArgs
ll <- c("longitude", "latitude")
druhy <- unique(as.vector(unlist(ndop.sp.selected$DRUH)))

set.seed(123)
druhy <- sample(druhy)
speciesParts <- split(druhy, ceiling(seq_along(druhy) / ndop.fs$speciesPerGroup))

# "tradiční" TGOB - použít jako BG všechny unikátní per_pixel presence všech druhů (NDOP pokrývá cca 85% území s cca 7300 pixely (z 9700)), takže je to zde upočitatelné
tgob.trad <- ndop.ff.au.more.sf.POLE.pp.c %>%
    group_by(POLE) %>%
    slice_head(n = 1)

# "tradiční" TGOB per pole a druh - pro thin
tgob.trad.sp <- ndop.ff.au.more.sf.POLE.pp.c %>%
    group_by(POLE, DRUH) %>%
    slice_head(n = 1)

sp.group <- ndop.sp.selected.regrouped %>% filter(nthGroup == ndop.fs$speciesPart)



null.default <- generateRPall(pcamap6.lsdClip[[1]], nback = ndop.fs$bg)
first <- TRUE
for (druh in as.vector(sp.group$DRUH)) { # speciesParts[[ndop.fs$speciesPart]]
    collector <- list()
    collector.thin <- list()
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(druh)
    print(sp.group)
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    e.mx.all <- list()
    bg.col <- list()

    pres <- ndop.ff.au.more.sf.POLE.pp.c %>% filter(DRUH == druh)
    pres.au <- unique(as.vector(unlist(pres$AUTOR)))
    # "skill TGOB" - vyberu pro TGOB pouze nálezy těch autorů, kteří už daný druh pozorovali
    # Ostatní jsou nejistí, respektive jejich sampling effort je irelevantní, protože není jisté, zda-li:
    # - druh vůbec poznají,
    # - jsou ochotní ho zaznamenávat nebo
    # - vůbec nemapují v oblasti výskytu druhu
    tgob.puv <- ndop.ff.au.more.sf.POLE.pp.c %>% filter(AUTOR %in% pres.au)
    tgob <- ndop.ff.au.more.sf %>% filter(AUTOR %in% pres.au)

    tgob.extract <- ndop.ff.au.more.sf.extract %>% filter(AUTOR %in% pres.au)

    pres.unique <- pres %>%
        group_by(POLE) %>%
        slice_head(n = 1)
    # skill verze tradičního TGOB - všechny unikátní pixely do BG // TODO zbytečné, stačilo by použít níže tgob.count!!!
    tgob.unique <- tgob %>%
        group_by(POLE) %>%
        slice_head(n = 1)

    # # počet (unikátních, zrušit znásobení spoluautory) nálezů na pixel
    # tgob.count <- tgob %>%
    #     group_by(POLE) %>%
    #     mutate(pp = length(unique(ID_NALEZ))) %>%
    #     slice_head(n = 1)
    #
    # pl <- tgob.count %>% dplyr::select(pp)

    df.temp <- as.data.frame(st_coordinates(pres.unique))
    names(df.temp) <- ll




    # bf.v0 <- rasterize(pl, pcamap6.lsdClip[[1]], field = "pp", fun = "last", background = NA, mask = FALSE) # nakonec přímo nepoužívám, protože tento (normalizovaný (0-1) per_pixel počet nálezů) odpovídá smoothingu adjust=0.01 (asi ještě ověřit individuálně pro všechny druhy?)
    #
    # # normalizace
    # r.min <- minValue(bf.v0)
    # r.max <- maxValue(bf.v0)
    # bf.v0 <- ((bf.v0 - r.min) / (r.max - r.min))
    # bf.v0 <- setMinMax(bf.v0)

    # pcamap6.skill <- stack(raster::mask(pcamap6, bf.v0)) # raster jen v místech dostupnych skill autorum daneho druhu

    pcamap6.skill <- mask(pcamap6.lsdClip[[1]], tgob.unique)

    sq2rad.dist <- ndop.fs$sq2rad[1] + 0.00001 # delší strana kvadrátu, přičtení drobné vzdálenosti
    sq2rad.dist.range <- sq2rad.dist * ndop.fs$sq2radDist

    thinned.versions <- as.character(17:22) # ????????????? budou všechny
    tgob.trad.tm <- tgob.trad.sp %>% filter(DRUH == druh)
    tgob.trad.tm.df <- as.data.frame(st_coordinates(tgob.trad.tm))
    colnames(tgob.trad.tm.df) <- c("x", "y")

    prmtr <- list(
        "base" = list(
            "r.cz" = pcamap6.lsdClip[[1]], "r.ssos" = pcamap6.skill,
            "v.all" = ndop.ff.au.more.sf.extract, "v.ssos" = tgob.extract,
            "fix.cz" = tgob.trad$geometry, "fix.ssos" = tgob.unique$geometry,
            "un" = null.default
        ),
        "nback" = list("5000" = ndop.fs$bg, "ssos" = nrow(tgob.unique)), # , "p.thin"=NULL "p" = nrow(pres.unique), "p.thin" = NULL
        "tg" = list("tgob" = ndop.ff.au.more.sf, "ssos" = tgob) # pro KSS
    )

    # prmtr <- list(
    #   "base" = list("r.cz" = NULL, "r.ssos" = NULL, "v.all"=NULL, "v.ssos"=NULL, "fix.cz" = NULL, "fix.ssos" = NULL),
    #   "nback" = list("5000" = NULL, "ssos" = NULL, "p"=NULL, "p.thin"=NULL ),
    #   "tg" = list("tgob" = NULL, "ssos" = NULL)
    # )

    # prmtr <- list(
    #   "base" = list("r.cz" = 11, "r.ssos" = 11111111, "v.all" = 222222, "v.ssos" = 22222222222, "fix.cz" = 2, "fix.ssos" = 22),
    #   "nback" = list("5000" = 3333, "ssos" = 3, "p" = 33, "p.thin" = 3333333),
    #   "tg" = list("tgob" = 444444, "ssos" = 4)
    # )


    for (p.base in names(prmtr[["base"]])) {
        print("p.base ++++++++++++++++++++++++++++++++++")

        for (p.nback in names(prmtr[["nback"]])) {
            print("p.nback ++++++++++++++++")
            for (p.tg in names(prmtr[["tg"]])) {
                print("p.tg +++++++")
                print("KSS*****************************")
                collector.temp <- list()
                print(c(p.base, p.nback, p.tg))
                ###
                # KSS
                ###
                id <- paste(c(p.base, p.nback, p.tg), collapse = "_")
                #  kss
                if (sum(p.base %in% c("r.cz", "r.ssos")) > 0) {
                    collector.temp <- smoothingRP(prmtr[["base"]][[p.base]], ndop.fs$adjusts, prmtr[["tg"]][[p.tg]], nback = prmtr[["nback"]][[p.nback]])
                    collector[[id]] <- collector.temp
                }
            }
            ###
            # RSS
            ###
            if (sum(p.base %in% c("v.all", "v.ssos", "un")) > 0) {
                print("RSS*****************************")
                collector.temp <- list()
                id <- paste(c(p.base, p.nback), collapse = "_")
                if (sum(p.base %in% c("v.all", "v.ssos")) > 0) {
                    collector.temp[["0"]] <- generateRPall(prmtr[["base"]][[p.base]], nback = prmtr[["nback"]][[p.nback]])
                } else {
                    collector.temp[["0"]] <- prmtr[["base"]][[p.base]]
                }
                if (p.base == "un") {
                    # thinning  presencí pro UN
                    for (thinDist in sq2rad.dist.range) {
                        occs.thinned <- ecospat.occ.desaggregation(xy = tgob.trad.tm.df, min.dist = thinDist, by = NULL)
                        thinDist <- as.character(thinDist)
                        collector.thin[[id]][[thinDist]] <- (occs.thinned %>% st_as_sf(coords = c("x", "y"), crs = 4326))$geometry
                    }
                }
                collector[[id]] <- collector.temp
            }
        }
        ###
        # area
        ###
        if (sum(p.base %in% c("fix.cz", "fix.ssos")) > 0) {
            print("area*****************************")
            collector.temp <- list()
            id <- p.base
            collector.temp[["0"]] <- prmtr[["base"]][[p.base]]
            collector[[id]] <- collector.temp
        }
    }


    # přidání clip verze všeho (bez area SSOS a RSS SSOS)

    for (id in names(collector)) {
        id.names <- unlist(strsplit(id, "_"))

        if (sum(c("fix.ssos", "v.ssos", "clip") %in% id.names) > 0) {
            next
        }
        for (adjust in names(collector[[id]])) {
            for (duplORnot in names(collector[[id]][[adjust]])) {
                print("clip*****************************")
                ex.temp <- extract(pcamap6.skill, st_coordinates(collector[[id]][[adjust]][[duplORnot]]))
                collector[[paste0(id, "_clip")]][[adjust]][[duplORnot]] <- collector[[id]][[adjust]][[duplORnot]][!is.na(ex.temp)]
            }
        }
    }
    gc()

    # run vsech variant BG s ENMeval
    for (id in names(collector)) {
        id.names <- unlist(strsplit(id, "_"))
        print(id)
        for (adjust in names(collector[[id]])) {
            print(adjust)
            for (duplORnot in names(collector[[id]][[adjust]])) {
                print(duplORnot)
                bg.temp <- as.data.frame(st_coordinates(collector[[id]][[adjust]][[duplORnot]]))
                names(bg.temp) <- ll
                print("základní:")
                e.mx.all[[druh]][[id]][[adjust]][[duplORnot]] <- ENMevaluate(
                    occs = df.temp,
                    envs = pcamap6,
                    bg = bg.temp,
                    algorithm = "maxnet", partitions = "randomkfold",
                    partition.settings = list("kfolds" = 3),
                    tune.args = tune.args
                )
                if (sum(c("un") %in% id.names) > 0) {
                    print("thin:")
                    for (thinDist in names(collector.thin[[id]])) {
                        df.temp.thin <- as.data.frame(st_coordinates(collector.thin[[id]][[thinDist]]))
                        names(df.temp.thin) <- ll
                        # přidat thinning presence
                        print("měním presenční dataset pro thinnovací verze")

                        e.mx.all[[druh]][[paste0(id, "_", as.character(thinDist))]][[adjust]][[duplORnot]] <- ENMevaluate(
                            occs = df.temp.thin,
                            envs = pcamap6,
                            bg = bg.temp,
                            algorithm = "maxnet", partitions = "randomkfold",
                            partition.settings = list("kfolds" = 3),
                            tune.args = tune.args
                        )
                    }
                }
            }
        }
    }


    gc()


    saveRDS(e.mx.all, paste0(path.wd.prep, "skillXXX_", druh, "_", ndop.fs$speciesPart, ".rds"))
    gc()
}


end_time <- Sys.time()
print(end_time - start_time)

print(sp.group)