cmd_arg <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "geosphere", "ENMeval", "sdm") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

# "C:\Program Files\R\R-4.2.1\bin\x64\Rscript.exe" "C:\Users\petr\Documents\2023-03-20\RP\RP\R\ndopTGOBeval.R"

# nastavit working directory
path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/projects-data/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.rgee <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/rgee/" # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "dataPrep/ndopTGOBeval/")
path.wd.prep.ndop <- paste0(path.wd, "dataPrep/ndop/")
path.wd.prep.tgobEval <- paste0(path.wd, "dataPrep/ndopTGOBeval/")
path.temp.res <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/modelsExport/merged/"
source(paste0(path.wd, "shared.R"))
############
# inputs
############
# uložit si results, rastery a udělat evaluate; ukládat i počty presencí a bg (z modelů) a WKT presencí (NDOP i LSD)
lsd.pa.centroids <- readRDS(paste0(path.wd, "dataPrep/lsd/lsd.pa.centroids.rds"))

# synonyma (mám funkci v rgee)
# evaluace všech predikcí z ENMeval přes LSD PA pomocí sdm::evaluates(observed, predicted)


############
# settings
############

ndop.fs <- list("months" = c(4:6), "years" = c(2019:2022), "precision" = 1000, "version" = "v1")

############
# execution
############

# #
# # spojení částí
# #

# # načtení RDS Z ENMeval
# rds_list <-
#     list.files(
#         path.wd.prep.tgobEval,
#         pattern = "^out.*\\.rds$",
#         ignore.case = TRUE,
#         full.names = TRUE
#     )

# first <- TRUE
# # merge částí do jednoho rds
# for (fpath in rds_list) {
#     rds.l <- list()
#     # dočasně rozdělené do dvou částí, nutno spojit:
#     rds.out.temp <- readRDS(fpath)
#     rds.r.temp <- readRDS(str_replace(fpath, "out.t", "rds.r"))

#     if (first) {
#         first <- FALSE
#         rds.out <- rds.out.temp
#         rds.r <- rds.r.temp
#     } else {
#         rds.out %<>% add_row(rds.out.temp)
#         rds.r <- append(rds.r, rds.r.temp)
#     }
# }

# rds.out <- as_tibble(rds.out)
# saveRDS(rds.out, paste0(path.wd.prep, "all_rds.out.rds"))
# saveRDS(rds.r, paste0(path.wd.prep, "all_rds.r.rds"))




rds.out <- readRDS(paste0(path.wd.prep, "all_rds.out.rds"))
rds.r <- readRDS(paste0(path.wd.prep, "all_rds.r.rds"))
rds.out %<>% filter(!is.na(AUC)) # odstranit modely kde nešla provést LSD evaluace

topXpc <- 0.01
topXrows <- 3
# výběr nejlepších X procent
rds.out.best <- rds.out %>%
    group_by(species) %>%
    mutate(AUC.max = max(AUC)) %>%
    ungroup() %>%
    filter(AUC >= (AUC.max - topXpc))


# kolik bylo úspěšných modelů pro zvolený treshold
mCounts <- rds.out.best %>%
    group_by(species) %>%
    summarise(n = length(species), AUC.max = max(AUC.max))


for (sp in unique(rds.out.best$species)) {
    rds.out.ndop <- rds.out %>%
        filter(species == sp & version == 1) %>%
        slice_head(n = 1) # potřebuju presence bez thinningu

    sp.selected <- rds.out.best %>%
        filter(species == sp) %>%
        slice_max(AUC, n = topXrows) %>%
        dplyr::select(species, AUC, tune.args, version, adjust, duplORnot, AUC.max, occs.wkt, p.wkt, a.wkt)
    print(sp.selected)
    r.stack.list <- list()
    for (i in 1:nrow(sp.selected)) {
        sp.selected.row <- sp.selected[i, ]
        r.stack.list[[i]] <- rds.r[[sp]][[sp.selected.row$version]][[sp.selected.row$adjust]][[sp.selected.row$duplORnot]][[sp.selected.row$tune.args]]
    }
    r.stack <- stack(r.stack.list)
    r.median <- calc(r.stack, fun = median)

    r.min <- minValue(r.median)
    r.max <- maxValue(r.median)
    r.median <- ((r.median - r.min) / (r.max - r.min))
    r.median <- setMinMax(r.median)

    pu <- st_as_sf(rds.out.ndop, wkt = "occs.wkt", crs = st_crs(4326))
    pp <- st_as_sf(sp.selected, wkt = "p.wkt", crs = st_crs(4326))
    pa <- st_as_sf(sp.selected, wkt = "a.wkt", crs = st_crs(4326))

    ggpl <- ggplot() +
        geom_raster(data = as.data.frame(r.median, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
        scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
        labs(fill = "hab. suitab. \n(cloglog)", caption = paste0("occs source: AOPK NDOP (years: ", min(ndop.fs$years), "-", max(ndop.fs$years), "; months: ", min(ndop.fs$months), "-", max(ndop.fs$months), ") | author: Petr Balej | WGS84 | grid: KFME (SitMap_2Rad) | generated: ", Sys.Date())) +
        theme_light() +
        geom_sf(data = pu, pch = 1, size = 2, color = "#34282C") +
        geom_sf(data = pp, pch = 20, size = 1.5, color = "green") +
        geom_sf(data = pa, pch = 20, size = 1.5, color = "#DC381F") +
        theme(
            plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5)
        ) +
        ggtitle(
            label = sp,
            subtitle = paste0("AUC=", sp.selected.row$AUC.max, " | topXpc=", topXpc, "; topXrows=", topXrows)
        )

    png(paste0(path.wd.prep, "preds/", sp, "_topXpc-", topXpc, "_topXrows-", topXrows, ".png"), width = 1400, height = 800)
    print(ggpl)
    dev.off()
}


end_time <- Sys.time()
print(end_time - start_time)