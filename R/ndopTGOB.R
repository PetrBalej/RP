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
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "geosphere", "ENMeval") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

# nastavit working directory
path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/projects-data/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.rgee <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/rgee/" # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
# source(paste0(path.rgee, "R/export_raster/functions.R"))
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

# per pixel - unikátní hodnoty druhů na pixel a autora
ndop.ff.au.more.sf.POLE.pp <- readRDS(paste0(path.wd.prep.ndop, "ndop.ff.au.more.sf.POLE.pp.rds"))
# jen body (centroidy pixelů)
ndop.ff.au.more.sf.POLE.pp.c <- readRDS(paste0(path.wd.prep.ndop, "ndop.ff.au.more.sf.POLE.pp.c.rds"))
# vybrané druhy
ndop.sp.selected <- readRDS(paste0(path.wd.prep.ndop, "ndop.sp.selected.rds"))


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
    "adjusts" = c(0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5),
    "tuneArgs" = list(fc = c("L", "LQ"), rm = c(1, 2, 3, 5, 10)),
    "bg" = 5000, "speciesPerGroup" = 9, "speciesOccMin" = 30,
    "speciesPart" = cmd_arg, "version" = "v1",
    "versionNames" = vn, "versionSmooting" = vf
)

############
# execution
############

ndop.sp.selected %<>% filter(n >= ndop.fs$speciesOccMin)

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

generateRPall <- function(raster_stack.bias, nback = 5000, random = FALSE) {
    # https://github.com/danlwarren/ENMTools/blob/004a4a1e182127900a5f62bc015770479bcd0415/R/check.bg.R#L138-L144
    # raster::sampleRandom ani sample(.., replace=FALSE) nejsou vhodné/nefungují pro menší N (clipnuté rastery skillem)
    background.points <- as.data.frame(rasterToPoints(raster_stack.bias))
    out <- list()
    if (random == FALSE) {
        # a1)
        inds <- sample(1:nrow(background.points),
            size = nback,
            prob = background.points[, 3],
            replace = TRUE
        )
        bg.temp <- background.points[inds, 1:2]
        out["dupl"] <- bg.temp %>% st_as_sf(coords = c("x", "y"), crs = 4326)

        # a2)
        out["uniq"] <- bg.temp %>%
            distinct() %>%
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

        # b2)
        out["uniq"] <- bg.temp %>%
            distinct() %>%
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
    ppp.deg.distance <- bw.scott.iso(res.ndop.coords.ppp)
    prague <- c(14.4, 50.0)
    # vzdálenosti v jednotkovém ppp na jedné lat a lon
    d.lon <- distHaversine(c(prague[1], prague[2] - ppp.deg.distance), c(prague[1], prague[2]))
    d.lat <- distHaversine(c(prague[1] - ppp.deg.distance, prague[2]), c(prague[1], prague[2]))
    # 0.001 a 0.01 jsou 100% korelované
    for (adj in adjusts) {
        raster_stack.bias <- resample(raster(density.ppp(res.ndop.coords.ppp, sigma = bw.scott.iso(res.ndop.coords.ppp), adjust = adj)), raster.template, method = "bilinear")
        crs(raster_stack.bias) <- rcrs
        raster_stack.bias <- mask(crop(raster_stack.bias, extent(raster.template)), raster.template)
        raster_stack.bias <- setMinMax(raster_stack.bias)

        # normalizace
        r.min <- minValue(raster_stack.bias)
        r.max <- maxValue(raster_stack.bias)
        raster_stack.bias <- ((raster_stack.bias - r.min) / (r.max - r.min))
        raster_stack.bias <- setMinMax(raster_stack.bias)
        bg[[as.character(adj * 100)]] <- generateRPall(raster_stack.bias, nback)
    }

    return(list(
        # zdrojové rastery neukládám, zbytečně velké...
        "bg" = bg, "bw" = ppp.deg.distance, "d.lon" = d.lon, "d.lat" = d.lat
    ))
}


# klasický celkový TGOB pro všechny druhy naráz
bf.all <- smoothingRP(pcamap6[[1]], ndop.fs$adjusts, ndop.ff.au.more.sf, nback = ndop.fs$bg)
bf.all$bg[["0"]] <- generateRPall(pcamap6[[1]], nback = ndop.fs$bg, random = TRUE) # přidám random body pro "null" model

saveRDS(bf.all, paste0(path.wd.prep, "bf.all-", cmd_arg, ".rds")) # asi zbytečné ukládat duplicitně pro každou skupinu...

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

sp.group <- ndop.sp.selected.regrouped %>% filter(nthGroup == ndop.fs$speciesPart)

first <- TRUE
for (druh in as.vector(sp.group$DRUH)) { # speciesParts[[ndop.fs$speciesPart]]
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
    tgob <- ndop.ff.au.more.sf.POLE.pp.c %>% filter(AUTOR %in% pres.au)
    pres.unique <- pres %>%
        group_by(POLE) %>%
        slice_head(n = 1)
    # skill verze tradičního TGOB - všechny unikátní pixely do BG // TODO zbytečné, stačilo by použít níže tgob.count!!!
    tgob.unique <- tgob %>%
        group_by(POLE) %>%
        slice_head(n = 1)

    # počet (unikátních, zrušit znásobení spoluautory) nálezů na pixel
    tgob.count <- tgob %>%
        group_by(POLE) %>%
        mutate(pp = length(unique(ID_NALEZ))) %>%
        slice_head(n = 1)

    pl <- tgob.count %>% dplyr::select(pp)

    df.temp <- as.data.frame(st_coordinates(pres.unique))
    names(df.temp) <- ll

    bf.v0 <- rasterize(pl, pcamap6[[1]], field = "pp", fun = "last", background = NA, mask = FALSE) # nakonec přímo nepoužívám, protože tento (normalizovaný (0-1) per_pixel počet nálezů) odpovídá smoothingu adjust=0.01 (asi ještě ověřit individuálně pro všechny druhy?)

    # normalizace
    r.min <- minValue(bf.v0)
    r.max <- maxValue(bf.v0)
    bf.v0 <- ((bf.v0 - r.min) / (r.max - r.min))
    bf.v0 <- setMinMax(bf.v0)

    pcamap6.skill <- stack(raster::mask(pcamap6, bf.v0)) # raster jen v místech dostupnych skill autorum daneho druhu


    ############ 1) tgob.all
    bg.col[["1"]] <- bf.all

    ############ 2) tgob.all orez skillem
    # ořezat celkovy BG skill arealem
    bf.all.skill <- list()
    for (adjust in names(bf.all$bg)) {
        for (duplORnot in names(bf.all$bg[[adjust]])) {
            ex.temp <- extract(pcamap6.skill[[1]], st_coordinates(bf.all$bg[[adjust]][[duplORnot]]))
            bf.all.skill$bg[[adjust]][[duplORnot]] <- bf.all$bg[[adjust]][[duplORnot]][!is.na(ex.temp)]
        }
    }
    bg.col[["2"]] <- bf.all.skill


    # TODO následující (3-8) zjednodušit co cyklu!

    ############ 3+4+5) tgob.skill cela CR s variantami počtu TGOB
    bf.skill.v1 <- smoothingRP(pcamap6[[1]], ndop.fs$adjusts, tgob, nback = ndop.fs$bg)
    bf.skill.v1$bg[["0"]] <- generateRPall(pcamap6[[1]], nback = ndop.fs$bg, random = TRUE)
    bg.col[["3"]] <- bf.skill.v1
    bf.skill.v2 <- smoothingRP(pcamap6[[1]], ndop.fs$adjusts, tgob, nback = nrow(tgob))
    bf.skill.v2$bg[["0"]] <- generateRPall(pcamap6[[1]], nback = nrow(tgob), random = TRUE)
    bg.col[["4"]] <- bf.skill.v2
    bf.skill.v3 <- smoothingRP(pcamap6[[1]], ndop.fs$adjusts, tgob, nback = nrow(pres.unique))
    bf.skill.v3$bg[["0"]] <- generateRPall(pcamap6[[1]], nback = nrow(pres.unique), random = TRUE)
    bg.col[["5"]] <- bf.skill.v3

    ############ 6+7+8) tgob.skill orez jen na skill s variantami počtu TGOB
    bf.skill.v1s <- smoothingRP(pcamap6.skill[[1]], ndop.fs$adjusts, tgob, nback = ndop.fs$bg)
    bf.skill.v1s$bg[["0"]] <- generateRPall(pcamap6.skill[[1]], nback = ndop.fs$bg, random = TRUE)
    bg.col[["6"]] <- bf.skill.v1s
    bf.skill.v2s <- smoothingRP(pcamap6.skill[[1]], ndop.fs$adjusts, tgob, nback = nrow(tgob))
    bf.skill.v2s$bg[["0"]] <- generateRPall(pcamap6.skill[[1]], nback = nrow(tgob), random = TRUE)
    bg.col[["7"]] <- bf.skill.v2s
    bf.skill.v3s <- smoothingRP(pcamap6.skill[[1]], ndop.fs$adjusts, tgob, nback = nrow(pres.unique))
    bf.skill.v3s$bg[["0"]] <- generateRPall(pcamap6.skill[[1]], nback = nrow(pres.unique), random = TRUE)
    bg.col[["8"]] <- bf.skill.v3s

    ############ 9) tradiční TGOB ze všech bodů nefitované + omezené skillem
    bf.trad <- list()
    bf.trad$bg[["tgobAll"]][["uniq"]] <- tgob.trad$geometry
    bf.trad$bg[["tgobSkill"]][["uniq"]] <- tgob.unique$geometry
    bg.col[["9"]] <- bf.trad

    gc()
    # run vsech variant BG s ENMeval
    for (v1 in names(bg.col)) {
        for (v2 in names(bg.col[[v1]]$bg)) {
            for (v3 in names(bg.col[[v1]]$bg[[v2]])) {
                bg.temp <- as.data.frame(st_coordinates(bg.col[[v1]]$bg[[v2]][[v3]]))
                names(bg.temp) <- ll

                e.mx.all[[druh]][[v1]][[v2]][[v3]] <- ENMevaluate(
                    occs = df.temp, envs = pcamap6,
                    bg.coords = bg.temp,
                    algorithm = "maxnet", partitions = "randomkfold",
                    partition.settings = list("kfolds" = 3),
                    tune.args = tune.args
                )
            }
        }
    }

    saveRDS(bg.col, paste0(path.wd.prep, "skillBG_", druh, "_", ndop.fs$speciesPart, ".rds"))
    saveRDS(e.mx.all, paste0(path.wd.prep, "skill_", druh, "_", ndop.fs$speciesPart, ".rds"))
    gc()
}


end_time <- Sys.time()
print(end_time - start_time)

print(sp.group)