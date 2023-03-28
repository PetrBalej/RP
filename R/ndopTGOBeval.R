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
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "geosphere", "ENMeval", "sdm") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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

# "C:\Program Files\R\R-4.2.1\bin\x64\Rscript.exe" "C:\Users\petr\Documents\2023-03-20\RP\RP\R\ndopTGOBeval.R"
# Rscript /mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/R/ndopTGOBeval.R 

# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- paste0(path.wd, "../projects-data/") # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.rgee <- paste0(path.wd, "../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "../dataPrep/ndopTGOBeval/")
path.wd.prep.ndop <- paste0(path.wd, "../dataPrep/ndop/")
path.wd.prep.tgobEval <- paste0(path.wd, "../dataPrep/ndopTGOBeval/")
path.temp.res <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/modelsExport29/merged/"
source(paste0(path.wd, "R/shared.R"))
############
# inputs
############
# uložit si results, rastery a udělat evaluate; ukládat i počty presencí a bg (z modelů) a WKT presencí (NDOP i LSD)
lsd.pa.centroids <- readRDS(paste0(path.wd, "../dataPrep/lsd/lsd.pa.centroids.rds"))

# synonyma (mám funkci v rgee)
# evaluace všech predikcí z ENMeval přes LSD PA pomocí sdm::evaluates(observed, predicted)


############
# settings
############

ndop.fs <- list("version" = "v1")

############
# execution
############

# sjednocení synonym LSD
lsd.pa.centroids %<>% rename(species = TaxonNameLAT)
lsd.pa.centroids.syn <- synonyms_unite(lsd.pa.centroids)
lsd.pa.centroids.syn.species <- unlist(unique(lsd.pa.centroids.syn$species))
# načtení RDS Z ENMeval
rds_list <-
    list.files(
        path.temp.res,
        pattern = paste0("^ssos_"), # !!! skillc je samostatná dodělávka všech 116 druhů dohromady 20-22 variant
        ignore.case = TRUE,
        full.names = TRUE
    )

# do nastavení!!!
tmp.group <- 4
tmp.perGroup <- ceiling(length(rds_list) / tmp.group)
tmp.group.parts <- split(rds_list, ceiling(seq_along(rds_list) / tmp.perGroup))
# omezím původní seznam
rds_list <- tmp.group.parts[[cmd_arg]]
print(rds_list)



# cyklíme...
for (fpath in rds_list) {
first <- TRUE
out.t <- list()
rds.r <- list()
    rds.l <- list()
    # dočasně rozdělené do dvou částí, nutno spojit:
    rds.l1 <- readRDS(fpath)
    # rds.l2 <- readRDS(str_replace(fpath, "skill", "skillb")) # !!! zakomentovat

    sp <- names(rds.l1)
    if (length(sp) > 1) {
        print("Více druhů v jednom souboru!")
        stop()
    }
    #  rds.l[[sp]] <- append(rds.l1[[sp]], rds.l2[[sp]])
    rds.l[[sp]] <- rds.l1[[sp]] # !!!

    # nezávislá PA data z LSD pro daný druh
    lsd.temp <- lsd.pa.centroids.syn %>% filter(species == sp)

    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(sp)
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    for (version in names(rds.l[[sp]])) {
        print("version: -------------------------")
        print(version)
        print("version: -------------------------")
        for (adjust in names(rds.l[[sp]][[version]])) {
            print("adjust: ____________")
            print(adjust)
            for (duplORnot in names(rds.l[[sp]][[version]][[adjust]])) {
                print(duplORnot)

                second <- TRUE
                for (layer in names(rds.l[[sp]][[version]][[adjust]][[duplORnot]]@predictions)) {
                    print(layer)
                    r.temp <- rds.l[[sp]][[version]][[adjust]][[duplORnot]]@predictions[[layer]]
                    rds.r[[sp]][[version]][[adjust]][[duplORnot]][[layer]] <- r.temp
                    ex.predicted <- extract(r.temp, st_coordinates(lsd.temp))
                    ev <- NA
                    #
                    # dopočíst a přidat další metriky z performance() (rgee)!!!
                    # tam ale není nic co zohledňuje tn tp (bez poměru k fp)
                    #

                    if (inherits(try({
                        ev <- sdm::evaluates(lsd.temp$presence, ex.predicted)
                    }), "try-error")) {
                        ev <- NA
                    }
                    if (is.na(ev)) {
                        next
                    }
                    ev.temp <- as.data.frame(ev@statistics[-3])
                    ev.temp <- merge(ev.temp, t(as.data.frame(ev@statistics$COR)))
                    ev.temp <- merge(ev.temp, ev@threshold_based[2, ]) # max(se+sp)
                    ev.temp <- as.data.frame(ev@statistics[-3])
                    ev.temp <- merge(ev.temp, t(as.data.frame(ev@statistics$COR)))
                    # max(se+sp)
                    ev.temp <- merge(ev.temp, ev@threshold_based[2, ])
                    ev.temp[["tune.args"]] <- layer
                    if (second) {
                        second <- FALSE
                        out.t.second <- ev.temp
                    } else {
                        out.t.second %<>% add_row(ev.temp)
                    }
                }
                rds.l[[sp]][[version]][[adjust]][[duplORnot]]@results[["species"]] <- sp
                rds.l[[sp]][[version]][[adjust]][[duplORnot]]@results[["version"]] <- version
                rds.l[[sp]][[version]][[adjust]][[duplORnot]]@results[["adjust"]] <- adjust
                rds.l[[sp]][[version]][[adjust]][[duplORnot]]@results[["duplORnot"]] <- duplORnot

                temp.t <- rds.l[[sp]][[version]][[adjust]][[duplORnot]]@results %>% left_join(out.t.second, by = "tune.args")

                occs <- rds.l[[sp]][[version]][[adjust]][[duplORnot]]@occs
                temp.t$occs.n <- nrow(occs)
                temp.t$occs.wkt <- st_as_text(st_combine(occs %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)))
                temp.t$bg.n <- nrow(rds.l[[sp]][[version]][[adjust]][[duplORnot]]@bg)
                # LSD
                temp.t$p.n <- nrow(lsd.temp %>% filter(presence == 1))
                temp.t$a.n <- nrow(lsd.temp %>% filter(presence == 0))
                temp.t$p.wkt <- st_as_text(st_combine(lsd.temp %>% filter(presence == 1)))
                temp.t$a.wkt <- st_as_text(st_combine(lsd.temp %>% filter(presence == 0)))

                if (first) {
                    first <- FALSE
                    out.t <- temp.t
                    col.names <- TRUE
                    if (file.exists(paste0(path.wd.prep, "ssos.out.t.csv"))) {
                        # pokud už paralelně v jiném procesu soubor vznikl, nevkládám znovu názvy sloupců
                        col.names <- FALSE
                    }
                    write.table(temp.t, paste0(path.wd.prep, "ssos.out.t.csv"), sep = ",", row.names = FALSE, col.names = col.names, append = TRUE)
                } else {
                    out.t %<>% add_row(temp.t)
                    write.table(temp.t, paste0(path.wd.prep, "ssos.out.t.csv"), sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
                }


                gc()
            }
        }
    }


    # ukládat per species - nutno - jinak je výsledné rds s rastery predikcí extrémně velké a nenačitatelné/neuložitelné
    saveRDS(out.t, paste0(path.wd.prep, "ssos.out.t.c_", sp, "_", as.character(cmd_arg), ".rds")) # !!! nové názvy, nepřepsat původní
    saveRDS(rds.r, paste0(path.wd.prep, "ssos.rds.r.c_", sp, "_", as.character(cmd_arg), ".rds")) # !!!
}



# stop()
# vn <- versionNames()
# vn9 <- unlist(strsplit(vn[["9"]], "\\+"))
# vn[["9"]] <- vn9[1]
# vn[["10"]] <- vn9[2]


end_time <- Sys.time()
print(end_time - start_time)
