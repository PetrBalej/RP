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
path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/"
setwd(path.wd)
path.data <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/projects-data/"
path.rgee <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/rgee/" # samsung500ntfs # paste0(path.expand("~"), "/Downloads/rgee2/rgee")
# source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "dataPrep/ndopTGOB/")
path.wd.prep.ndop <- paste0(path.wd, "dataPrep/ndop/")


############
# inputs
############

sitmap_2rad.czechia <- readRDS(paste0(path.wd, "dataPrep/sitmap_2rad/sitmap_2rad-czechia.rds")) # 4326

pcamap6 <- readRDS(paste0(path.wd, "dataPrep/observerSkill/pcamap6.rds"))
vifstep2 <- saveRDS(paste0(path.wd, "dataPrep/observerSkill/vifstep2.rds"))
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
# názvy variant TGB?
ndop.fs <- list("adjusts" = c(0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5), "bg" = 5000, "speciesPerGroup" = 9, "speciesPart" = cmd_arg, "version" = "v1")


generateRP <- function(raster_stack.bias, nback = 5000) {
    # https://github.com/danlwarren/ENMTools/blob/004a4a1e182127900a5f62bc015770479bcd0415/R/check.bg.R#L138-L144
    background.points <- as.data.frame(rasterToPoints(raster_stack.bias))
    inds <- sample(1:nrow(background.points),
        size = nback,
        prob = background.points[, 3],
        replace = TRUE
    )
    bg.temp <- background.points[inds, 1:2]
    return(bg.temp %>% st_as_sf(coords = c("x", "y"), crs = 4326))
}

generateBiasFile <- function(raster.template, adjusts, occs.sf, nback = 5000) {
    #
    # celkový TGOB
    #

    # raster.template[!is.na(raster.template)] <- 1
    # raster.template[is.na(raster.template)] <- 0
    # raster.template.matrix <- as.matrix(raster.template)
    # raster.template.matrix[raster.template.matrix == 1] <- "TRUE"
    # raster.template.matrix[raster.template.matrix == 0] <- "FALSE"

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

        bg[[as.character(adj * 100)]] <- generateRP(raster_stack.bias, nback)

        #  raster_stack.bias.col[[as.character(adj * 100)]] <- raster_stack.bias
    }
    # celkový bias file
    # saveRDS(raster_stack.bias.col, paste0(path.wd.prep, "raster_stack.bias.col.rds"))
    return(list(
        # "r" = raster_stack.bias.col,
        "bg" = bg, "bw" = ppp.deg.distance, "d.lon" = d.lon, "d.lat" = d.lat
    ))
}


bf.all <- generateBiasFile(pcamap6[[1]], ndop.fs$adjusts, ndop.ff.au.more.sf, nback = ndop.fs$bg)
# saveRDS(bf.all, paste0(path.wd.prep, "bf.all.rds"))
# plot(bf.all$bg[["1"]])




#
# run per species models and varying BG
#
e.mx.all <- list()
# tune.args <- list(fc = c("L", "LQ"), rm = c(0.5, 1:5, 10))
tune.args <- list(fc = c("L", "LQ"), rm = c(1:3))
ll <- c("longitude", "latitude")
druhy <- unique(as.vector(unlist(ndop.sp.selected$DRUH)))

set.seed(123) 
druhy <- sample(druhy)
speciesParts <- split(druhy, ceiling(seq_along(druhy) / ndop.fs$speciesPerGroup))


first <- TRUE
for (druh in speciesParts[[ndop.fs$speciesPart]]) {
    print(druh)

    pres <- ndop.ff.au.more.sf.POLE.pp.c %>% filter(DRUH == druh)
    pres.au <- unique(as.vector(unlist(pres$AUTOR)))
    tgob <- ndop.ff.au.more.sf.POLE.pp.c %>% filter(AUTOR %in% pres.au)
    pres.unique <- pres %>%
        group_by(POLE) %>%
        slice_head(n = 1)
    tgob.unique <- tgob %>%
        group_by(POLE) %>%
        slice_head(n = 1)

    # počet (unikátních, zrušit znásobení spoluautory) nálezů na pixel
    tgob.count <- tgob %>%
        group_by(POLE) %>%
        mutate(pp = length(unique(ID_NALEZ))) %>%
        slice_head(n = 1)

    # pl <- tgob.count %>% dplyr::select(pp)

    df.temp <- as.data.frame(st_coordinates(pres.unique))
    names(df.temp) <- ll


    bf.v0 <- rasterize(pl, pcamap6[[1]], field = "pp", fun = "last", background = NA, mask = FALSE)

    # normalizace
    r.min <- minValue(bf.v0)
    r.max <- maxValue(bf.v0)
    bf.v0 <- ((bf.v0 - r.min) / (r.max - r.min))
    bf.v0 <- setMinMax(bf.v0)
    bf.v0.bg <- generateRP(bf.v0, ndop.fs$bg)


    #
    # 0) per pixel weighted  BG (podle BG skills)
    #
    print("0) ---------------------------------------")
    bg.temp <- as.data.frame(st_coordinates(bf.v0.bg))
    names(bg.temp) <- ll

    e.mx.all[[druh]][[0]][["0"]] <- ENMevaluate(
        occs = df.temp, envs = pcamap6, # pcamap6
        bg.coords = bg.temp,
        algorithm = "maxnet", partitions = "randomkfold",
        partition.settings = list("kfolds" = 3),
        tune.args = tune.args
    )


    #
    # 1) per pixel unique  BG (podle BG skills) - získat unikátní pixely kde někdo byl
    #
    print("1) ---------------------------------------")
    bg.temp <- as.data.frame(st_coordinates(tgob.count))
    names(bg.temp) <- ll
    e.mx.all[[druh]][[1]][["0"]] <- ENMevaluate(
        occs = df.temp, envs = pcamap6,
        bg.coords = bg.temp,
        algorithm = "maxnet", partitions = "randomkfold",
        partition.settings = list("kfolds" = 3),
        tune.args = tune.args
    )

    #
    # 2) celkový TGB
    #
    print("2) ---------------------------------------")
    for (adjust in names(bf.all$bg)) {
        bg.temp <- as.data.frame(bf.all$bg[[adjust]])
        names(bg.temp) <- ll

        e.mx.all[[druh]][[2]][[adjust]] <- ENMevaluate(
            occs = df.temp, envs = pcamap6,
            bg.coords = bg.temp,
            algorithm = "maxnet", partitions = "randomkfold",
            partition.settings = list("kfolds" = 3),
            tune.args = tune.args
        )
    }

    #
    # 3)  TGB (podle BG skills)
    #
    print("3) ---------------------------------------")

    bf.ss <- generateBiasFile(pcamap6[[1]], ndop.fs$adjusts, tgob, nback = ndop.fs$bg)


    for (adjust in names(bf.ss$bg)) {
        bg.temp <- as.data.frame(bf.ss$bg[[adjust]])
        names(bg.temp) <- ll

        e.mx.all[[druh]][[3]][[adjust]] <- ENMevaluate(
            occs = df.temp, envs = pcamap6,
            bg.coords = bg.temp,
            algorithm = "maxnet", partitions = "randomkfold",
            partition.settings = list("kfolds" = 3),
            tune.args = tune.args
        )
    }

    gc()
    saveRDS(e.mx.all, paste0(path.wd.prep, "ee_", druh, ".rds"))


    # df.temp <- data.frame(
    #     "species" = druh,
    #     "p" = nrow(pres), "pu" = nrow(pres.unique),
    #     "tgob" = nrow(tgob), "tgobu" = nrow(tgob.unique)
    # )
    # if (first) {
    #     first <- FALSE
    #     df.out <- df.temp
    # } else {
    #     df.out %<>% add_row(df.temp)
    # }



    # ggpl <- ggplot() +
    #     geom_sf(data = SPH_STAT.source.s, fill = "LightBlue") +
    #     geom_sf(data = pl, aes(colour = pp), pch = 15, size = 2) +
    #     labs(color = "TGOB (log10)", caption = paste0("occs source: AOPK NDOP (years: ", min(ndop.fs$years), "-", max(ndop.fs$years), "; months: ", min(ndop.fs$months), "-", max(ndop.fs$months), ") | author: Petr Balej | WGS84 | grid: KFME (SitMap_2Rad) | generated: ", Sys.Date())) +
    #     scale_colour_gradient(low = "white", high = "red", trans = "log10") +
    #     theme_light() +
    #     geom_sf(data = pres.unique) +
    #     theme(
    #         plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    #         plot.subtitle = element_text(hjust = 0.5)
    #     ) +
    #     ggtitle(
    #         label = (druh),
    #         subtitle = paste0("TGOB (same species observers occs per pixel) | ", nrow(pres.unique), " / ", nrow(tgob.unique))
    #     )

    # png(paste0(path.wd.prep, "tgob/tgob-", str_replace_all(druh, "[\\/\\:]", "_"), ".png"), width = 1200, height = 700)

    # print(ggpl)
    # dev.off()
}


# saveRDS(df.out, paste0(path.wd.prep, "df.out.rds"))
# write.csv(df.out, paste0(path.wd.prep, "df.out.csv"))




# uložení rozdělených dat per species pro model

# vygenerovat i mapy presencí druhu (počet) a backgroundu podle odpovídajících autorů (počty i per pixel)



############
# execution
############



# minimum x presenci a y absencí - možné až po získání prediktorů a "děr"

# + remove spatialy autocorrelated squares - geographically/environmentally
# alespoň reprezentativnost vůči altitude a landcoveru?

# do protokolu časy i adresu skriptu, který to vygeneroval + název vygenerovaného souboru

# https://github.com/PetrBalej/rgee/blob/master/R/clean/export_raster.R
# https://github.com/PetrBalej/rgee/blob/master/R/igaD.R

end_time <- Sys.time()
print(end_time - start_time)