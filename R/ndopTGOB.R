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


# *** Running ENMeval v2.0.4 with maxnet from maxnet package v0.1.4 *** -- na WS1
# "C:/Users/svc.kagup/AppData/Local/Programs/R/R-4.2.1/bin/x64/Rscript.exe" "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/R/ndopTGOB.R" 1
start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "ENMeval", "ecospat") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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

# NDOP originální nedofiltrované pro TGOB
ndop <- readRDS(paste0(path.ndop, "ndop.rds"))
ndop.POLE <- readRDS(paste0(path.ndop, "ndop.POLE.rds"))
# NDOP dofiltrované pro SSOS
# NDOP druhy presence
ndopP <- readRDS(paste0(path.ndop, "ndopP.rds"))
ndopP.POLE <- readRDS(paste0(path.ndop, "ndopP.POLE.rds"))

ndop.topAUTOR100 <- readRDS(paste0(path.ndop, "ndop.topAUTOR100.rds"))
ndop.topAUTOR10 <- readRDS(paste0(path.ndop, "ndop.topAUTOR10.rds"))

ndop.topDRUH50 <- readRDS(paste0(path.ndop, "ndop.topDRUH50.rds"))
ndop.topDRUH10 <- readRDS(paste0(path.ndop, "ndop.topDRUH10.rds"))


# LSD PA
# lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds"))  %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci
lsd.pa.min <- readRDS(paste0(path.lsd, "lsd.pa.min.rds"))

# statistiky NDOP druhy
ndop.stat.res.sp <- readRDS(paste0(path.ndop, "ndop.stat.res.sp.rds"))



# # souřadnice na SF - původní počty pro PPP
# ndop.ff.au.more.sf <- readRDS(paste0(path.ndop, "ndop.ff.au.more.sf.rds"))
# ndop.ff.au.more.sf.POLE <- readRDS(paste0(path.ndop, "ndop.ff.au.more.sf.POLE.rds"))
# ndop.ff.au.more.sf %<>% arrange(ID_NALEZ)
# ndop.ff.au.more.sf.POLE %<>% arrange(ID_NALEZ)
# # sum(ndop.ff.au.more.sf$ID_NALEZ == ndop.ff.au.more.sf.POLE$ID_NALEZ)
# ndop.ff.au.more.sf$POLE <- ndop.ff.au.more.sf.POLE$POLE


# # per pixel - unikátní hodnoty druhů na pixel a autora
# ndop.ff.au.more.sf.POLE.pp <- readRDS(paste0(path.ndop, "ndop.ff.au.more.sf.POLE.pp.rds"))
# # jen body (centroidy pixelů)
# ndop.ff.au.more.sf.POLE.pp.c <- readRDS(paste0(path.ndop, "ndop.ff.au.more.sf.POLE.pp.c.rds"))

# # vybrané druhy (musím znovu přepočítat po odebrání LSD)
# ndop.ff.au.more.sf.POLE.pp.sp <- readRDS(paste0(path.ndop, "ndop.ff.au.more.sf.POLE.pp.sp.rds"))

# ### LSD pro filtraci druhů které lze ověřit
# lsd.pa.min <- readRDS(paste0(path.prep, "lsd/lsd.pa.min.rds"))
# lsd.POLE <- st_read(paste0(path.prep, "lsd/overview-selected-2rad.shp"))


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
    "tuneArgs" = list(fc = c("L", "LQ"), rm = c(1, 2, 3, 5, 10)),
    "bgRatio" = 1 / 4,
    "speciesPerGroup" = 7, "speciesOccMin" = 30,
    "sq2rad" = c((1 / 6) / 4, 0.1 / 4), # kvadráty KFME 2rad, xy velikost ve stupních
    "sq2radDist" = c(1:5),
    "speciesPart" = cmd_arg, "version" = "v1",
    "versionNames" = vn, "versionSmooting" = vf
)
# dopočty
ndop.fs$bg <- round(czechia.sq * ndop.fs$bgRatio)

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


generateRPall <- function(RasterLayer_OR_sf_POINT, random = TRUE, nBackRatio = 1 / 3) {
    # počítám s tím, že dostanu už unikátní body per pixel!!!
    # pak můžu vždy poměrově vybrat část k "ploše"

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

    total <- nrow(background.points)
    proportion <- total * nBackRatio

    probs <- NULL
    if (random == FALSE) {
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

smoothingRP <- function(raster.template, adjusts, occs.sf, nBackRatio = 1 / 3) {
    rcrs <- crs(raster.template)
    ext <- extent(raster.template)
    ow <- owin(xrange = c(ext@xmin, ext@xmax), yrange = c(ext@ymin, ext@ymax))

    res.ndop.coords <- st_coordinates(occs.sf)
    res.ndop.coords.ppp <- ppp(res.ndop.coords[, 1], res.ndop.coords[, 2], window = ow)

    # raster_stack.bias.col <- list()
    bg <- list()
    # ppp.deg.distance <- bw.scott(res.ndop.coords.ppp)
    # prague <- c(14.4, 50.0)
    # # vzdálenosti v jednotkovém ppp na jedné lat a lon
    # d.lon <- distHaversine(c(prague[1], prague[2] - ppp.deg.distance[2]), c(prague[1], prague[2]))
    # d.lat <- distHaversine(c(prague[1] - ppp.deg.distance[1], prague[2]), c(prague[1], prague[2]))
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
        bg[[as.character(adj)]] <- generateRPall(raster_stack.bias, nBackRatio, random = FALSE)
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
druhy <- unique(as.vector(unlist(ndop.stat.res.sp.selected$DRUH)))



# "tradiční" TGOB - použít jako BG všechny unikátní per_pixel presence všech druhů (NDOP pokrývá cca 85% území s cca 7300 pixely (z 9700)), takže je to zde upočitatelné
tgob.trad <- ndop.POLE %>%
    group_by(POLE) %>%
    slice_head(n = 1) %>%
    st_centroid() %>%
    dplyr::select(-everything())

# "tradiční" TGOB per pole a druh - pro thin
tgob.trad.sp <- ndopP.POLE %>%
    group_by(POLE, DRUH) %>%
    st_centroid() %>%
    slice_head(n = 1)

sp.group <- ndop.stat.res.sp.selected.regrouped %>% filter(nthGroup == ndop.fs$speciesPart)



null.default <- generateRPall(predictors[[1]], nBackRatio = ndop.fs$bgRatio)
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

    pres <- ndopP %>% filter(DRUH == druh)
    pres.au <- unique(as.vector(unlist(pres$AUTOR)))
    # "skill TGOB" - vyberu pro TGOB pouze nálezy těch autorů, kteří už daný druh pozorovali
    # Ostatní jsou nejistí, respektive jejich sampling effort je irelevantní, protože není jisté, zda-li:
    # - druh vůbec poznají,
    # - jsou ochotní ho zaznamenávat nebo
    # - vůbec nemapují v oblasti výskytu druhu

    tgob <- ndopP %>% filter(AUTOR %in% pres.au)

    pres.unique <- pres %>%
        group_by(POLE) %>%
        slice_head(n = 1)

    # skill verze tradičního TGOB - všechny unikátní pixely do BG
    tgob.unique <- ndopP.POLE %>%
        filter(AUTOR %in% pres.au) %>%
        group_by(POLE) %>%
        st_centroid() %>%
        slice_head(n = 1) %>%
        dplyr::select(-everything())

    df.temp <- as.data.frame(st_coordinates(pres.unique))
    names(df.temp) <- ll


    predictors.ssos <- mask(predictors[[1]], tgob.unique)

    sq2rad.dist <- ndop.fs$sq2rad[1] + 0.00001 # delší strana kvadrátu, přičtení drobné vzdálenosti
    sq2rad.dist.range <- sq2rad.dist * ndop.fs$sq2radDist

    # buffer kolem presencí
    pres.unique.buffer <- st_as_sf(st_union(st_buffer(pres.unique, 10000)))
    predictors.ssos.buffer <- mask(predictors[[1]], pres.unique.buffer)


    tgob.trad.tm <- tgob.trad.sp %>% filter(DRUH == druh)
    tgob.trad.tm.df <- as.data.frame(st_coordinates(tgob.trad.tm))
    colnames(tgob.trad.tm.df) <- c("x", "y")



    #
    # ssos - omezení autorů s "nepřirozeným"" poměrem sběru daného druhu
    #

    ### celkové
    ssos.DRUH <- as_tibble(tgob) %>%
        dplyr::select(-geometry) %>%
        group_by(DRUH) %>%
        summarise(POLE.n = n_distinct(POLE))
    ssos.DRUH.sum.all <- sum(ssos.DRUH$POLE.n)
    ssos.DRUH.sum <- as.numeric(ssos.DRUH %>% filter(DRUH == druh) %>% dplyr::select(POLE.n))
    # "přirozený" poměr
    ssos.DRUH.ratio <- (ssos.DRUH.sum / ssos.DRUH.sum.all) * 0.5 # snížím hranici na 30%, odstraní extrémy (lidi co druh značí jen občas)

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
        mutate(ratio = POLE.n.sum / POLE.n.sum.all) %>%
        mutate(ratioLess = ifelse(ratio < ssos.DRUH.ratio, 1, 0))

    # unname(quantile(ssos.DRUH.AUTOR.ratio$ratio, probs=0.25)
    ssos.natural.au <- ssos.DRUH.AUTOR.ratio %>% filter(ratioLess == 0)
    pres.au <- unique(ssos.natural.au$AUTOR)
    ssos <- ndop %>% filter(AUTOR %in% pres.au)
    predictors.ssos2 <- mask(predictors[[1]], ssos)


    tgob.unique2 <- ndopP.POLE %>%
        filter(AUTOR %in% pres.au) %>%
        group_by(POLE) %>%
        st_centroid() %>%
        slice_head(n = 1) %>%
        dplyr::select(-everything())

    prmtr <- list(
        "base" = list(
            "r.cz" = predictors[[1]], "r.ssos" = predictors.ssos, "r.ssos.2" = predictors.ssos2,
            "r.buffer" = predictors.ssos.buffer,
            "v.all" = tgob.trad, "v.ssos" = tgob.unique, "v.ssos.2" = tgob.unique2, # vybírám 1/3 (random nebo kss)
            "fix.cz" = tgob.trad, "fix.ssos" = tgob.unique, "fix.ssos.2" = tgob.unique2, # beru celé !!! nemám ve výsledcích, kde to neprojde? v generateRPall???
            "un" = null.default
        ),
        # "nback" = list("5000" = ndop.fs$bg, "ssos" = nrow(tgob.unique)), # , "p.thin"=NULL "p" = nrow(pres.unique), "p.thin" = NULL
        "tg" = list(
            "tgob" = ndop, "ssos" = tgob, "ssos.2" = ssos,
            "top.10a" = ndop.topAUTOR10, "top.100a" = ndop.topAUTOR100, "top.10s" = ndop.topDRUH10, "top.50s" = ndop.topDRUH50
        ) # pro KSS
    )


    for (p.base in names(prmtr[["base"]])) {
        print("p.base ++++++++++++++++++++++++++++++++++")

        # for (p.nback in names(prmtr[["nback"]])) {
        # print("p.nback ++++++++++++++++")
        for (p.tg in names(prmtr[["tg"]])) {
            print("p.tg +++++++")
            print("KSS*****************************")
            collector.temp <- list()
            print(c(p.base, p.tg))
            ###
            # KSS
            ###

            #  kss
            if (sum(p.base %in% c("r.cz", "r.ssos", "r.ssos.2")) > 0) {
                if (sum(p.base %in% c("r.cz")) > 0) {
                    id <- paste(c(p.base, p.tg), collapse = "_")
                    collector.temp <- smoothingRP(prmtr[["base"]][[p.base]], ndop.fs$adjusts, prmtr[["tg"]][[p.tg]], nBackRatio = ndop.fs$bgRatio)
                    collector[[id]] <- collector.temp
                } else {
                    if (sum(p.tg %in% c("tgob", "ssos", "ssos.2")) > 0) {
                        id <- paste(c(p.base, p.tg), collapse = "_")
                        collector.temp <- smoothingRP(prmtr[["base"]][[p.base]], ndop.fs$adjusts, prmtr[["tg"]][[p.tg]], nBackRatio = ndop.fs$bgRatio)
                        collector[[id]] <- collector.temp
                    }
                }
            }
        }

        ###
        # buffer
        ###


        if (sum(p.base %in% c("r.buffer")) > 0) {
            id <- paste(c(p.base, "kss"), collapse = "_")
            collector.temp <- smoothingRP(prmtr[["base"]][[p.base]], 1, pres.unique, nBackRatio = ndop.fs$bgRatio)
            collector[[id]] <- collector.temp
        }


        collector.temp <- list()
        if (sum(p.base %in% c("r.buffer")) > 0) {
            id <- paste(c(p.base, "rss"), collapse = "_")
            collector.temp[["0"]] <- generateRPall(prmtr[["base"]][[p.base]], nBackRatio = ndop.fs$bgRatio)
            collector[[id]] <- collector.temp
        }

        collector.temp <- list()


        ###
        # RSS
        ###
        if (sum(p.base %in% c("v.all", "v.ssos", "v.ssos.2", "un")) > 0) {
            print("RSS*****************************")
            collector.temp <- list()
            id <- paste(c(p.base, "rss"), collapse = "_")
            if (sum(p.base %in% c("v.all", "v.ssos", "v.ssos.2")) > 0) {
                collector.temp[["0"]] <- generateRPall(prmtr[["base"]][[p.base]], nBackRatio = ndop.fs$bgRatio)
            } else {
                collector.temp[["0"]] <- prmtr[["base"]][[p.base]] # !!!!! tady je problém, pro un_ssos dávám rovněž 5000 !!!!! měl bych to tu pro každý druh zvlášť nasamplovat podle počtu areA
            }
            collector[[id]] <- collector.temp
            if (p.base == "un") {
                # thinning  presencí pro UN
                for (thinDist in sq2rad.dist.range) {
                    occs.thinned <- ecospat.occ.desaggregation(xy = tgob.trad.tm.df, min.dist = thinDist, by = NULL)
                    thinDist <- as.character(thinDist)
                    collector.thin[[id]][[thinDist]] <- (occs.thinned %>% st_as_sf(coords = c("x", "y"), crs = 4326))$geometry
                }
            }
        }
        # }
        ###
        # area
        ###
        if (sum(p.base %in% c("fix.cz", "fix.ssos", "fix.ssos.2")) > 0) {
            print("area*****************************")
            collector.temp <- list()
            id <- p.base
            collector.temp[["0"]] <- prmtr[["base"]][[p.base]]
            collector[[id]] <- collector.temp
        }
    }


    # přidání clip verze všeho (bez area SSOS a RSS SSOS)

    for (id in names(collector)) {
        id.names <- unlist(strsplit(id, "_"))[1]

        if (sum(c("r.buffer") %in% id.names) > 0) {
            next
        }
        if (sum(c("r.ssos.2", "v.ssos.2", "fix.ssos.2") %notin% id.names) == 3) {
            for (adjust in names(collector[[id]])) {
                print("clip1*****************************")
                ex.temp <- extract(predictors.ssos, st_coordinates(collector[[id]][[adjust]]))
                collector[[paste0(id, "_clip.ssos")]][[adjust]] <- collector[[id]][[adjust]][!is.na(ex.temp), ]
            }
        }
        if (sum(c("r.ssos", "v.ssos", "fix.ssos") %notin% id.names) == 3) {
            for (adjust in names(collector[[id]])) {
                print("clip2*****************************")
                ex.temp <- extract(predictors.ssos2, st_coordinates(collector[[id]][[adjust]]))
                collector[[paste0(id, "_clip.ssos.2")]][[adjust]] <- collector[[id]][[adjust]][!is.na(ex.temp), ]
            }
        }
    }

  
    gc()


    # run vsech variant BG s ENMeval
    for (id in names(collector)) {
        id.names <- unlist(strsplit(id, "_"))[1]
        print(id)
        for (adjust in names(collector[[id]])) {
            print(adjust)

            bg.temp <- as.data.frame(st_coordinates(collector[[id]][[adjust]]))
            names(bg.temp) <- ll
            print("základní:")
            e.mx.all[[druh]][[id]][[adjust]] <- ENMevaluate(
                occs = df.temp,
                envs = predictors,
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

                    e.mx.all[[druh]][[paste0(id, "_", as.character(thinDist))]][[adjust]] <- ENMevaluate(
                        occs = df.temp.thin,
                        envs = predictors,
                        bg = bg.temp,
                        algorithm = "maxnet", partitions = "randomkfold",
                        partition.settings = list("kfolds" = 3),
                        tune.args = tune.args
                    )
                }
            }
        }
    }

    saveRDS(e.mx.all, paste0(path.tgob, "ssos_", druh, "_", ndop.fs$speciesPart, ".rds"))
    gc()
}


end_time <- Sys.time()
print(end_time - start_time)

print(sp.group)