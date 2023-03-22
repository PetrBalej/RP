start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "geosphere", "ENMeval", "sdm", "ggplot2", "grid", "gridExtra") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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

vn <- versionNames()

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



# ## napárování 20-22 rasterů
# rr.orig <- readRDS(paste0(path.wd.prep, "all_rds.r-copy.rds"))
# rr.v2022 <- readRDS(paste0(path.wd.prep, "all_rds.rc.rds"))
# rr.orig.names <- names(rr.orig)
# for (sp in rr.orig.names) {
#     print(sp)
#     for (version in names(rr.v2022[[sp]])) {
#         print(version)
#         rr.orig[[sp]][[version]] <- rr.v2022[[sp]][[version]]
#     }
# }
# # saveRDS(rr.orig, paste0(path.wd.prep, "all_rds.r.rds"))




rds.out <- readRDS(paste0(path.wd.prep, "all_rds.out.rds"))
rds.out %<>% filter(!is.na(AUC)) # odstranit modely kde nešla provést LSD evaluace


# porovnat 2 samostatné random null:
rds.out.null.compare <- rds.out %>%
    filter((version == 1 | version == 3) & adjust == 0) %>%
    group_by(species, version, duplORnot) %>%
    mutate(compareGroup = paste0(species, "-", version, "-", duplORnot)) %>%
    ungroup() %>%
    group_by(compareGroup) %>%
    slice_max(AUC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(species, duplORnot) %>%
    mutate(Diff = AUC - lag(AUC)) %>%
    na.omit()

visDiff1 <- rds.out.null.compare %>%
    filter(duplORnot == "dupl") %>%
    dplyr::select(Diff)
visDiff2 <- rds.out.null.compare %>%
    filter(duplORnot == "uniq") %>%
    dplyr::select(Diff)

boxplot(list("dupl" = visDiff1$Diff, "uniq" = visDiff2$Diff))
hist(abs(visDiff1$Diff))
hist(abs(visDiff2$Diff))

#
# dif vůči random
#

# nulová bg random verze
vs.null <- rds.out %>%
    filter(version == 3 & adjust == 0) %>% # vybírám 3, ale mohl bych i 1, reálně musím udělat více replikací (změn backgroundu) a zprůměrovat je
    group_by(species, duplORnot) %>%
    slice_max(AUC, n = 1, with_ties = FALSE) # zde neřeším více modelů se stejně vysokým AUC


# porovnání nej z jednotlivých základních verzí - neřeším adjusty ani šířku thinu!!!  (4176 / 2 / 116 = 18 verzí - mám ale 19 verzí - všechny nemají dupl/unique varianty - nemůžu dělit 2...)

vs.all <- rds.out %>%
    # filter(adjust != 0) %>% # ponechávám všechny interní random verze ("0") - neměl bych dělat porovnání s interními random verzemi?
    group_by(species, version, duplORnot) %>%
    slice_max(AUC, n = 1, with_ties = FALSE) # %>% # zde neřeším více modelů se stejně vysokým AUC
# mutate(compareGroup = paste0(species, "-", version, "-", duplORnot))

vs.all %<>% left_join(vs.null %>% dplyr::select(AUC), by = c("species", "duplORnot"), suffix = c("", ".null"))

vs.all.versions <- unique(vs.all$version)

first <- TRUE
for (v in vs.all.versions) {
    vs.all.temp <- vs.all %>%
        ungroup() %>%
        filter(version == v) %>%
        mutate(AUC.diff = AUC - AUC.null)

    if (first) {
        first <- FALSE
        vs.diff <- vs.all.temp
    } else {
        vs.diff %<>% add_row(vs.all.temp)
    }
}

# pivot_wider
vs.all.compare <- vs.diff %>%
    ungroup() %>%
    group_by(version, duplORnot)
# %>% summarise(AUC.diff.median = median(AUC.diff))

vs.all.compare.order <- vs.all.compare %>%
    summarise(AUC.diff.median = median(AUC.diff)) %>%
    ungroup() %>%
    ungroup() %>%
    arrange(AUC.diff.median)


group_ordered <- with(
    vs.all.compare, # Order boxes by median
    reorder(
        version,
        AUC.diff,
        median
    )
)
data_ordered <- vs.all.compare # Create data with reordered group levels
data_ordered$group <- factor(data_ordered$version,
    levels = levels(group_ordered)
)

# # # # # # udělat místo toho per species rozdíly dupl/uniq a ty dát do boxplotů
ggplot(vs.all.compare, aes(x = version, y = AUC.diff, fill = duplORnot)) +
    geom_boxplot()



stop()
rds.r <- readRDS(paste0(path.wd.prep, "all_rds.r.rds"))

topXpc <- 0.01
topXrows <- 1 # vybírám už jen nejlepší predikci pro každou verzi a duplikaci, neberu všechny...
# výběr nejlepších X procent

#
### # validace nad independent LSD (AUC)
#
# nejlepší z nejlepších
rds.out.best <- rds.out %>%
    group_by(species, version, duplORnot) %>%
    slice_max(AUC, n = topXrows, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(species) %>%
    filter(AUC >= (max(AUC) - topXpc)) %>%
    ungroup()

# "null" - random 5000 a všechny presence (verze 1 adjust 0 nebo verze 3 adjust 0 )
rds.out.null <- rds.out %>%
    filter(version == 1 & adjust == 0) %>%
    group_by(species, version, duplORnot) %>%
    slice_max(AUC, n = topXrows, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(species) %>%
    filter(AUC >= (max(AUC) - topXpc)) %>%
    ungroup()
#
### # validace nad interním blokem z NDOP (auc.val.avg)
#
# nejlepší z nejlepších
rds.out.best.ndop <- rds.out %>%
    group_by(species, version, duplORnot) %>%
    slice_max(auc.val.avg, n = topXrows, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(species) %>%
    filter(auc.val.avg >= (max(auc.val.avg) - topXpc)) %>%
    ungroup()

# "null" - random 5000 a všechny presence (verze 1 adjust 0 nebo verze 3 adjust 0 )
rds.out.null.ndop <- rds.out %>%
    filter(version == 1 & adjust == 0) %>%
    group_by(species, version, duplORnot) %>%
    slice_max(auc.val.avg, n = topXrows, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(species) %>%
    filter(auc.val.avg >= (max(auc.val.avg) - topXpc)) %>%
    ungroup()





# kolik bylo úspěšných modelů pro zvolený treshold
mCounts <- rds.out.best %>%
    group_by(species) %>%
    summarise(n = length(species), AUC.max = max(AUC.max))


for (sp in unique(rds.out.best$species)) {
    rds.out.ndop <- rds.out %>%
        filter(species == sp & version == 1) %>%
        slice_head(n = 1) # potřebuju presence bez thinningu

    ### # LSD validace
    sp.selected <- rds.out.best %>%
        filter(species == sp) %>%
        dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
    print(sp.selected)

    sp.selected.null <- rds.out.null %>%
        filter(species == sp) %>%
        dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
    print(sp.selected.null)

    ### # NDOP validace
    sp.selected.ndop <- rds.out.best.ndop %>%
        filter(species == sp) %>%
        dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
    print(sp.selected)

    sp.selected.null.ndop <- rds.out.null.ndop %>%
        filter(species == sp) %>%
        dplyr::select(species, AUC, auc.val.avg, tune.args, version, adjust, duplORnot, occs.wkt, p.wkt, a.wkt)
    print(sp.selected.null.ndop)
    #
    ### # validace nad independent LSD (AUC)
    #
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

    r.stack.list.null <- list()
    for (i in 1:nrow(sp.selected.null)) {
        sp.selected.row.null <- sp.selected.null[i, ]
        r.stack.list.null[[i]] <- rds.r[[sp]][[sp.selected.row.null$version]][[sp.selected.row.null$adjust]][[sp.selected.row.null$duplORnot]][[sp.selected.row.null$tune.args]]
    }
    r.stack.null <- stack(r.stack.list.null)
    r.median.null <- calc(r.stack.null, fun = median)

    r.min <- minValue(r.median.null)
    r.max <- maxValue(r.median.null)
    r.median.null <- ((r.median.null - r.min) / (r.max - r.min))
    r.median.null <- setMinMax(r.median.null)

    #
    ### # validace nad interním blokem z NDOP (auc.val.avg)
    #
    r.stack.list.ndop <- list()
    for (i in 1:nrow(sp.selected.ndop)) {
        sp.selected.row.ndop <- sp.selected.ndop[i, ]
        r.stack.list.ndop[[i]] <- rds.r[[sp]][[sp.selected.row.ndop$version]][[sp.selected.row.ndop$adjust]][[sp.selected.row.ndop$duplORnot]][[sp.selected.row.ndop$tune.args]]
    }
    r.stack.ndop <- stack(r.stack.list.ndop)
    r.median.ndop <- calc(r.stack.ndop, fun = median)

    r.min <- minValue(r.median.ndop)
    r.max <- maxValue(r.median.ndop)
    r.median.ndop <- ((r.median.ndop - r.min) / (r.max - r.min))
    r.median.ndop <- setMinMax(r.median.ndop)

    r.stack.list.null.ndop <- list()
    for (i in 1:nrow(sp.selected.null.ndop)) {
        sp.selected.row.null.ndop <- sp.selected.null.ndop[i, ]
        r.stack.list.null.ndop[[i]] <- rds.r[[sp]][[sp.selected.row.null.ndop$version]][[sp.selected.row.null.ndop$adjust]][[sp.selected.row.null.ndop$duplORnot]][[sp.selected.row.null.ndop$tune.args]]
    }
    r.stack.null.ndop <- stack(r.stack.list.null.ndop)
    r.median.null.ndop <- calc(r.stack.null.ndop, fun = median)

    r.min <- minValue(r.median.null.ndop)
    r.max <- maxValue(r.median.null.ndop)
    r.median.null.ndop <- ((r.median.null.ndop - r.min) / (r.max - r.min))
    r.median.null.ndop <- setMinMax(r.median.null.ndop)




    pu <- st_as_sf(rds.out.ndop, wkt = "occs.wkt", crs = st_crs(4326))
    pp <- st_as_sf(sp.selected, wkt = "p.wkt", crs = st_crs(4326))
    pa <- st_as_sf(sp.selected, wkt = "a.wkt", crs = st_crs(4326))

    #
    ### # validace nad independent LSD (AUC)
    #
    ggpl1 <- ggplot() +
        geom_raster(data = as.data.frame(r.median, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
        scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
        labs(fill = "hab. suitab. \n(cloglog)") +
        theme_light() +
        geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
        geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
        geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
        theme(
            legend.title = element_text(hjust = 0.5, size = 9),
            # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            # plot.subtitle = element_text(hjust = 0.5)
            plot.title = element_text(hjust = 0.5, size = 10)
        ) +
        ggtitle(
            label = paste0("AUC(NDOP)=", round(max(sp.selected$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected$AUC), digits = 2))
        )


    ggpl2 <- ggplot() +
        geom_raster(data = as.data.frame(r.median.null, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
        scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
        labs(fill = "hab. suitab. \n(cloglog)") +
        theme_light() +
        geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
        geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
        geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
        theme(
            legend.title = element_text(hjust = 0.5, size = 9),
            # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            # plot.subtitle = element_text(hjust = 0.5)
            plot.title = element_text(hjust = 0.5, size = 10)
        ) +
        ggtitle(
            label = paste0("AUC(NDOP)=", round(max(sp.selected.null$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.null$AUC), digits = 2))
        )

    #
    ### # validace nad interním blokem z NDOP (auc.val.avg)
    #
    ggpl3 <- ggplot() +
        geom_raster(data = as.data.frame(r.median.ndop, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
        scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
        labs(fill = "hab. suitab. \n(cloglog)") +
        theme_light() +
        geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
        geom_sf(data = pp, pch = 20, size = 0.5, color = "green", show.legend = "point") +
        geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F", show.legend = "point") +
        theme(
            legend.title = element_text(hjust = 0.5, size = 9),
            # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            # plot.subtitle = element_text(hjust = 0.5)
            plot.title = element_text(hjust = 0.5, size = 10)
        ) +
        ggtitle(
            label = paste0("AUC(NDOP)=", round(max(sp.selected.ndop$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.ndop$AUC), digits = 2))
        )

    ggpl4 <- ggplot() +
        geom_raster(data = as.data.frame(r.median.null.ndop, xy = TRUE) %>% na.omit(), aes(x = x, y = y, fill = layer)) +
        scale_fill_gradient2(low = "white", mid = "#DBE9FA", high = "#4682B4", midpoint = 0.5) +
        labs(fill = "hab. suitab. \n(cloglog)") +
        theme_light() +
        geom_sf(data = pu, pch = 20, size = 0.8, color = "#34282C") +
        geom_sf(data = pp, pch = 20, size = 0.5, color = "green") +
        geom_sf(data = pa, pch = 20, size = 0.5, color = "#DC381F") +
        theme(
            legend.title = element_text(hjust = 0.5, size = 9),
            # plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            # plot.subtitle = element_text(hjust = 0.5)
            plot.title = element_text(hjust = 0.5, size = 10)
        ) +
        ggtitle(
            label = paste0("AUC(NDOP)=", round(max(sp.selected.null.ndop$auc.val.avg), digits = 2), " | AUC(LSD)=", round(max(sp.selected.null.ndop$AUC), digits = 2))
        )


    png(paste0(path.wd.prep, "preds/", sp, "_topXpc-", topXpc, "_topXrows-", topXrows, ".png"), width = 1500, height = 900)
    print(
        grid.arrange(
            ggpl4,
            ggpl3,
            ggpl2,
            ggpl1,
            nrow = 2, ncol = 2,
            top = textGrob(paste0(sp, "\n no bias correction                                                                                             best bias correction"), gp = gpar(fontsize = 20, fontface = "bold")),
            left = textGrob("test indep. (AUC LSD selection)                  test dep. (AUC NDOP selection)", gp = gpar(fontsize = 20, fontface = "bold"), rot = 90),
            bottom = textGrob(
                paste0("NDOP: černá | LSD: zelená=presence, červená=absence) | topXpc=", topXpc, "; topXrows=", topXrows, "\n", paste(unlist(vn[which(names(vn) %in% unique(sp.selected$version))]), collapse = " | "), "\n occs source: AOPK NDOP (years: ", min(ndop.fs$years), "-", max(ndop.fs$years), "; months: ", min(ndop.fs$months), "-", max(ndop.fs$months), ") | author: Petr Balej | WGS84 | grid: KFME (SitMap_2Rad) | generated: ", Sys.Date())
            )
        )
    )
    dev.off()
}


end_time <- Sys.time()
print(end_time - start_time)