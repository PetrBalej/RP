start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr", "rgee")

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

path.wd <- paste0(gcfl(), "/../")

# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/"
setwd(path.wd)
path.data <- "../projects-data/"
path.rgee <- "../../rgee/" # samsung500ntfs # paste0(path.expand("~"), "/Downloads/rgee2/rgee")
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "../dataPrep/predictors/20192022/")

# použít Google Drive?
use_google_drive <- TRUE

ee_Initialize(email = "balej.petr@gmail.com", drive = use_google_drive)

# ee_user_info()

# při odpojení nebo zneplatnění původního přihlašovacího tokenu
# ee_clean_credentials()

# tempdir() # dočasný adresář pro aktuální R session


# # # # # # # # # # # # # # # # # # # # # #
# nastavení základních parametrů [start]  #
# # # # # # # # # # # # # # # # # # # # # #

## projekce výsledných rastrů, (většinou) není možné použít WGS84 (4326) - nesouhlasí pak rozlišení (odlišný počet rows/cols) pokud se používají různé zdroje (L8, WorldClim, SRTM, ...)
# 3035: ETRS89-LAEA	- Lambertovo azimutální stejnoploché zobrazení
# 32633: UTM zone 33N - použito Mercatorovo válcové konformní zobrazení (UTM zobrazení), základní poledník 15°
res_proj_epsg <- 4326

## výsledná velikost pixelu v m
scale <- 30

## výběr regionu

# definice obálek (bounding box) různě velkých území pro testování
sz_cechy <- list(
    xmin = 13.0,
    xmax = 13.5,
    ymin = 50.0,
    ymax = 50.5
)
cesko <- list(
    xmin = 12.0,
    xmax = 19.0,
    ymin = 48.5,
    ymax = 51.2
)
str_evropa <-
    list(
        xmin = 8.5,
        xmax = 22.0,
        ymin = 46.0,
        ymax = 53.5
    )
str_evropa2 <-
    list(
        xmin = 8.6,
        xmax = 21.9,
        ymin = 46.4,
        ymax = 53.1
    )

# výběr konkrétního území
bb <- cesko
# bb <- paste0(git_project_path, "/shp/ne_50m_admin_0_countries/czechia/cz_4326.shp")

# bb <- list(
#   xmin = 12.7,
#   xmax = 13.2,
#   ymin = 47.1,
#   ymax = 47.5
# )
## časové rozsahy

# pro paper1-observert skills

# rozsah snímků od/do
years_range <- list(from = "2019-01-01", to = "2022-12-31")

# rozsah jedné sezóny v měsících (podvýběr z vybraného období years_range výše), možno i více sezón
season_months_range <- list(c(4, 4), c(5, 5), c(6, 6))

# # # # # # # # # # # # # # # # # # # # # #
# nastavení základních parametrů [konec]  #
# # # # # # # # # # # # # # # # # # # # # #

# parametry použitých datasetů z GEE - export z gee_datasets/gee-pouzite-datasety.xlsx
gee_datasets_path_csv <-
    paste0(path.rgee, "gee_datasets/gee-pouzite-datasety.csv")

# načtení csv s datasety z GEE
gdl <- gee_datasets_list(gee_datasets_path_csv)

if (!is.character(bb)) {
    xmin <- bb$xmin
    xmax <- bb$xmax
    ymin <- bb$ymin
    ymax <- bb$ymax

    bb_geometry <- NULL
    bb_geometry_rectangle <- ee$Geometry$Rectangle(
        coords = c(xmin, ymin, xmax, ymax),
        proj = paste0("EPSG:", res_proj_epsg),
        geodesic = FALSE
    )
} else {
    if (file.exists(bb)) {
        bb_geometry_ee <- st_read(bb) %>% sf_as_ee()
        bb_geometry <- bb_geometry_ee$geometry() # [[1]]$getInfo()
        bb_geometry_rectangle <- bb_geometry_ee$geometry()$bounds()
    } else {
        stop(paste0("Shapefile ", bb, " not exist!"))
    }
}


sitmap_2rad <- list()
# Define a geometry
sitmap_2rad$N <- st_read(paste0(path.data, "aopk/sitmap_2rad/rgee/sitmap_2rad_4326_czechia_N.shp"))
sitmap_2rad$S <- st_read(paste0(path.data, "aopk/sitmap_2rad/rgee/sitmap_2rad_4326_czechia_S.shp"))
output_df <- list()

for (half in names(sitmap_2rad)) {
    print(half)
    for (season in season_months_range) {
        print(season)
        ################################################################
        # L8 _SR 'LANDSAT/LC08/C01/T1_SR' - raw bandy
        ################################################################

        # aplikace základních geografických, časových a odmračňovacích/odstiňovacích/aerosol/radio filtrů
        l8_sr_collection <-
            ee$ImageCollection(gdl$landsat2$geeSnippet)$
                filterBounds(bb_geometry_rectangle)$
                filterDate(years_range$from, years_range$to)$
                filter(
                ee$Filter$calendarRange(season[1], season[2], "month")
            )$map(mask_L8_sr2_qa)$map(mask_L8_sr2_radsat)$map(mask_L8_sr2_aerosol)$map(applyScaleFactors)





        bands_all <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B10")
        bn <- l8_sr_collection$first()$bandNames()$getInfo()
        bn.renamed <- gsub("SR_|ST_", "", bn)

        l8_sr_collection <- l8_sr_collection$select(bn, bn.renamed)

        #####
        # raw
        #####
        output_df[[paste0("l8_", scale, "_", season[1], "_raw_stdev")]] <- ee_extract(
            x = l8_sr_collection$select(bands_all)$median(),
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$stdDev(),
            sf = FALSE
        )

        output_df[[paste0("l8_", scale, "_", season[1], "_raw_mean")]] <- ee_extract(
            x = l8_sr_collection$select(bands_all)$median(),
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$median(),
            sf = FALSE
        )

        #####
        # ndvi
        #####
        output_df[[paste0("l8_", scale, "_", season[1], "_ndvi_stdev")]] <- ee_extract(
            x = l8_sr_collection$select(c("B5", "B4"))$median()$normalizedDifference(c("B5", "B4")),
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$stdDev(),
            sf = FALSE
        )

        output_df[[paste0("l8_", scale, "_", season[1], "_ndvi_mean")]] <- ee_extract(
            x = l8_sr_collection$select(c("B5", "B4"))$median()$normalizedDifference(c("B5", "B4")),
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$median(),
            sf = FALSE
        )


        #####
        # mndwi
        #####
        # https://www.linkedin.com/pulse/ndvi-ndbi-ndwi-calculation-using-landsat-7-8-tek-bahadur-kshetri
        output_df[[paste0("l8_", scale, "_", season[1], "_mndwi_stdev")]] <- ee_extract(
            x = l8_sr_collection$select(c("B3", "B6"))$median()$normalizedDifference(c("B3", "B6")),
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$stdDev(),
            sf = FALSE
        )

        output_df[[paste0("l8_", scale, "_", season[1], "_mndwi_mean")]] <- ee_extract(
            x = l8_sr_collection$select(c("B3", "B6"))$median()$normalizedDifference(c("B3", "B6")),
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$median(),
            sf = FALSE
        )


        #####
        # evi
        #####
        #   # https://www.usgs.gov/core-science-systems/nli/landsat/landsat-enhanced-vegetation-index?qt-science_support_page_related_con=0
        #   # In Landsat 8, EVI = 2.5 * ((Band 5 – Band 4) / (Band 5 + 6 * Band 4 – 7.5 * Band 2 + 1)).

        evi_prep <- l8_sr_collection$select(c("B5", "B4", "B2"))$median()

        evi_res <- evi_prep$expression(
            "2.5 * ((B5 - B4) / (B5 + 6 * B4 - 7.5 * B2 + 1))",
            list(
                B5 = evi_prep$select("B5"),
                B4 = evi_prep$select("B4"),
                B2 = evi_prep$select("B2")
            )
        )

        output_df[[paste0("l8_", scale, "_", season[1], "_evi_stdev")]] <- ee_extract(
            x = evi_res,
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$stdDev(),
            sf = FALSE
        )

        output_df[[paste0("l8_", scale, "_", season[1], "_evi_mean")]] <- ee_extract(
            x = evi_res,
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$median(),
            sf = FALSE
        )

        #####
        # msavi
        #####
        #   # https://www.usgs.gov/core-science-systems/nli/landsat/landsat-modified-soil-adjusted-vegetation-index
        #   # In Landsat 8, MSAVI = (2 * Band 5 + 1 – sqrt ((2 * Band 5 + 1)2 – 8 * (Band 5 – Band 4))) / 2.

        msavi_prep <- l8_sr_collection$select(c("B5", "B4"))$median()
        msavi_res <- msavi_prep$expression(
            "(2 * B5 + 1 - sqrt(pow((2 * B5 + 1), 2) - 8 * (B5 - B4))) / 2",
            list(
                B5 = msavi_prep$select("B5"),
                B4 = msavi_prep$select("B4")
            )
        )

        output_df[[paste0("l8_", scale, "_", season[1], "_msavi_stdev")]] <- ee_extract(
            x = msavi_res,
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$stdDev(),
            sf = FALSE
        )

        output_df[[paste0("l8_", scale, "_", season[1], "_msavi_mean")]] <- ee_extract(
            x = msavi_res,
            y = sitmap_2rad[[half]],
            scale = scale,
            fun = ee$Reducer$median(),
            sf = FALSE
        )
    }
    saveRDS(output_df, file = paste0(path.wd.prep, "kfme16-l8-", half, "_czechia_wc_l8_4-6.rds"))
    output_df <- list()
    gc()

    # scale_bio <- 927.67
    # ################################################################
    # # Worldclim/Bioclim 'WORLDCLIM/V1/BIO'
    # ################################################################

    wc <- applyScaleFactorsBio(ee$Image(gdl$worldclim$geeSnippet))

    bands_all <- wc$bandNames()$getInfo()
    # bands_all <- c("bio01", "bio02")


    # # Extract values
    output_df[[paste0("wc_", scale, "_", season[1], "_stdev")]] <- ee_extract(
        x = wc$select(bands_all[1:10]),
        y = sitmap_2rad[[half]],
        scale = scale,
        fun = ee$Reducer$stdDev(),
        sf = FALSE
    )
    print("BIOa")
    output_df[[paste0("wc_", scale, "_", season[1], "_mean")]] <- ee_extract(
        x = wc$select(bands_all[1:10]),
        y = sitmap_2rad[[half]],
        scale = scale,
        fun = ee$Reducer$median(),
        sf = FALSE
    )
    saveRDS(output_df, file = paste0(path.wd.prep, "kfme16-bioA-", half, "_czechia_wc_l8_4-6.rds"))
    output_df <- list()
    gc()


    output_df[[paste0("wc_", scale, "_", season[1], "_stdev")]] <- ee_extract(
        x = wc$select(bands_all[11:19]),
        y = sitmap_2rad[[half]],
        scale = scale,
        fun = ee$Reducer$stdDev(),
        sf = FALSE
    )
    print("BIOb")
    output_df[[paste0("wc_", scale, "_", season[1], "_mean")]] <- ee_extract(
        x = wc$select(bands_all[11:19]),
        y = sitmap_2rad[[half]],
        scale = scale,
        fun = ee$Reducer$median(),
        sf = FALSE
    )
    saveRDS(output_df, file = paste0(path.wd.prep, "kfme16-bioB-", half, "_czechia_wc_l8_4-6.rds"))
    output_df <- list()
    gc()


    # saveRDS(output_df, file = paste0(path.wd.prep, "test-kfme16-N_czechia_wc_l8_2018-2021_4-6.rds"))
    # saveRDS(output_df, file = paste0(path.wd.prep, "clean/predictors/kfme16-", half, "_czechia_wc_l8_4-6.rds"))
}

end_time <- Sys.time()

# časové rozmezí a celkový čas generování
print(paste("start:", start_time))
print(paste("konec:", end_time))
print(end_time - start_time)


#########
# spojení obou "polovin" ČR a vepsání vypočtených variant hodnot pixelů do rasterů
#########



stop()
# po dokopírování do N a S adresářů

rds_list_N <-
    list.files(
        path.wd.prep,
        pattern = paste0("[N][_]"),
        ignore.case = TRUE,
        full.names = TRUE
    )
rds_list_S <-
    list.files(
        path.wd.prep,
        pattern = paste0("[S][_]"),
        ignore.case = TRUE,
        full.names = TRUE
    )


first <- TRUE
parts.out <- list()
# sloučit N
for (nN in rds_list_N) {
    nNf <- readRDS(nN)
    nS <- str_replace(nN, "-N_", "-S_")
    nSf <- readRDS(nS)
    print(nN)

    for (l1 in names(nNf)) {
        print(l1)
        nNf[[l1]] %<>% add_row(nSf[[l1]])
    }

    parts.out[[nN]] <- nNf
}


# dopřidání bioclim dvou částí prediktorů
bioName <- "wc_30_6_mean"
bioName.cv <- "wc_30_6_stdev"
parts.a <- append(
    list(
        "wc_30_6_mean" = append(parts.out[[1]][[bioName]], parts.out[[2]][[bioName]][, -1]),
        "wc_30_6_stdev" = append(parts.out[[1]][[bioName.cv]], parts.out[[2]][[bioName.cv]][, -1])
    ),
    parts.out[[3]]
)


kfme16 <- parts.a

# ####### nové načítání L8 a WC
# kfme16.N <- readRDS(paste0(path.wd.prep, "kfme16-N_czechia_wc_l8_4-6.rds"))
# kfme16.S <- readRDS(paste0(path.wd.prep, "kfme16-S_czechia_wc_l8_4-6.rds"))

# # sloučit N a S
# for (l1 in names(kfme16.N)) {
#     kfme16.N[[l1]] %<>% add_row(kfme16.S[[l1]])
# }

# kfme16 <- kfme16.N

sf.grid <- st_read(paste0(path.wd, "dataPrep/sitmap_2rad/sitmap_2rad-czechia.shp"))


# vzorový KFME raster pro 16 subkvadrátů
r <- raster(
    xmn = 12.0834157383896326, # set minimum x coordinate
    xmx = 18.8750299183168408, # set maximum x coordinate
    ymn = 48.5500817521775048, # set minimum y coordinate
    ymx = 51.0750755551042701, # set maximum y coordinate
    res = c((1 / 6) / 4, 0.1 / 4), # resolution in c(x,y) direction
    vals = runif(16463, 0, 1)
)


for (l1 in names(kfme16)) {
    kfme16.t <- as_tibble(kfme16[[l1]])

    print(l1)

    names(kfme16[[l1]])[-1]

    for (band in names(kfme16[[l1]])[-1]) {
        print(band)
        kfme16.t.select <- kfme16.t %>% dplyr::select(all_of(c("POLE", band)))

        m <- merge(sf.grid, kfme16.t.select, by = "POLE")

        rr <- rasterize(m, r, field = band)
        crs(rr) <- res_proj_epsg

        writeRaster(rr, paste0(path.wd.prep, "kfme16-", l1, "-", band, ".tif"), format = "GTiff", overwrite = TRUE)

        # dodělávka CV
        if (grepl("stdev", l1, fixed = TRUE)) {
            l1.mean <- gsub("stdev", "mean", l1)
            print(l1.mean)
            kfme16.t.mean <- as_tibble(kfme16[[l1.mean]])
            kfme16.t.mean.select <- kfme16.t.mean %>% dplyr::select(all_of(c("POLE", band)))

            m.mean <- merge(sf.grid, kfme16.t.mean.select, by = "POLE")

            rr.mean <- rasterize(m.mean, r, field = band)
            crs(rr.mean) <- res_proj_epsg

            writeRaster(rr / rr.mean, paste0(path.wd.prep, "kfme16-", gsub("stdev", "cv", l1), "-", band, ".tif"), format = "GTiff", overwrite = TRUE)
        }
    }
}

#
# vizuální kontrola - odstranění artefaktů a nedokonalých spojení bandů
#
raster_stack <-
    rasters_dir_stack(
        path.wd.prep,
        "tif"
    )

tif_names <-
    list.files(
        path.wd.prep,
        pattern = paste0("tif$"),
        ignore.case = TRUE,
        full.names = FALSE
    )

# kontrola po tom co budu rozlišovat mean a cv
names(raster_stack) <- str_replace_all(tif_names, c("\\.tif" = "", "kfme16-" = "", "raw-" = "", "_30" = "", "-nd" = "", "-constant" = "", "6_mean-bio" = "0_mean-bio", "6_stdev-bio" = "0_stdev-bio", "6_cv-bio" = "0_cv-bio"))


raster_stack <- dropLayer(raster_stack, grep("B10|stdev", names(raster_stack), ignore.case = TRUE))
length(names(raster_stack)) # 52

rr <- writeRaster(raster_stack, paste0(path.wd.prep, "preds_base.grd"), format = "raster", overwrite = TRUE)
hdr(rr, format = "ENVI")
saveRDS(raster_stack, paste0(path.wd.prep, "preds_base.rds"))





pdf(paste0(path.wd.prep, "check_visually.pdf"), width = 12, height = 9)
for (l in sort(names(raster_stack))) {
    # plot(normalize(raster_stack[[l]]), main=l)
    plot(raster_stack[[l]], main = l)
}
dev.off()

# dropnout květen - artefakty s nepovedeným švem v JČ a některé bio také artefakty nebo nepřirozené gradienty
# dropnout červen (2014-2017) - artefakty s nepovedeným švem v Kralický Sněžník a některé bio také artefakty nebo nepřirozené gradienty
raster_stack <- dropLayer(raster_stack, grep("_5_|bio03|bio08|bio09|bio15|cv.bio", names(raster_stack), ignore.case = TRUE))
length(names(raster_stack)) # 60


pdf(paste0(path.wd.prep, "check_visually_final.pdf"), width = 12, height = 9)
for (l in sort(names(raster_stack))) {
    # plot(normalize(raster_stack[[l]]), main=l)
    plot(raster_stack[[l]], main = l)
}
dev.off()

rr <- writeRaster(raster_stack, paste0(path.wd.prep, "preds_final.grd"), format = "raster", overwrite = TRUE)
hdr(rr, format = "ENVI")
saveRDS(raster_stack, paste0(path.wd.prep, "preds_final.rds"))

stop()

#
# VIF
#

# raster_stack <- readRDS(paste0(path.wd.prep, "preds_final.rds"))

library(usdm)

# vifstep
vs <- vifstep(raster_stack, th = 3)
# vifcor - jen pro dofiltr, asi není nutné
vs.vc <- vifcor(raster_stack[[vs@results$Variables]], th = 0.7)
print(vs.vc@results$Variables)

raster_stack_vifs <- raster_stack[[vs.vc@results$Variables]]
# raster_stack_vifs <- readRDS(paste0(path.wd.prep, "preds_final_vifs.rds"))

raster_stack_vifs <- raster::mask(raster_stack_vifs, sum(raster_stack_vifs))
raster_stack_vifs <- raster_stack_vifs2
replaceExtreme <- function(x) {
    # především pro úlety z CV
    quantiles <- quantile(x, c(.0001, .9999))
    x[x < quantiles[1]] <- quantiles[1]
    x[x > quantiles[2]] <- quantiles[2]
    return(x)
}

# odstranění extrémů a normalizace hodnot 0-1
for (rn in names(raster_stack_vifs)) {
    raster_stack_vifs[[rn]] <- replaceExtreme(raster_stack_vifs[[rn]])

    raster_stack_vifs[[rn]] <- raster_stack_vifs[[rn]] / cellStats(raster_stack_vifs[[rn]], stat = sum)

    r.min <- minValue(raster_stack_vifs[[rn]])
    r.max <- maxValue(raster_stack_vifs[[rn]])
    raster_stack_vifs[[rn]] <- ((raster_stack_vifs[[rn]] - r.min) / (r.max - r.min))
    raster_stack_vifs[[rn]] <- setMinMax(raster_stack_vifs[[rn]])
}
# #  hist(raster_stack_vifs[[1]])
plot(raster_stack_vifs)
vif(raster_stack_vifs)



rr <- writeRaster(raster_stack_vifs, paste0(path.wd.prep, "preds_final_vifs.grd"), format = "raster", overwrite = TRUE)
hdr(rr, format = "ENVI")
saveRDS(raster_stack_vifs, paste0(path.wd.prep, "preds_final_vifs.rds"))



# library(RStoolbox)
# # PCA vs.vc
# pcamap <- rasterPCA(raster_stack, spca = TRUE)
# summary(pcamap$model)