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

path.wd <- paste0(gcfl(), "/")
# "C:\Program Files\R\R-4.2.1\bin\x64\Rscript.exe" "D:\PersonalWork\Balej\v2\RP\RP\R\ndopTGOBeval.R" 1

# Rscript /mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/R/ndopTGOBeval.R
# "C:\Program Files\R\R-4.1.2\bin\x64\Rscript.exe" "C:\Users\balej\Documents\2023-04-06\RP\RP\R\ndopTGOBeval.R" 1

# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- paste0(path.wd, "../../projects-data/") # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.prep <- paste0(path.wd, "../../dataPrep/")
source(paste0(path.wd, "shared.R"))
path.rgee <- paste0(path.wd, "../../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.lsd <- paste0(path.prep, "lsd/")
path.models <- paste0(path.prep, "ndopTGOB/")
path.eval <- paste0(path.prep, "ndopTGOBeval/")
############
# inputs
############

# uložit si results, rastery a udělat evaluate; ukládat i počty presencí a bg (z modelů) a WKT presencí (NDOP i LSD)
lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds")) %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci
# lsd.pa.centroids <- readRDS(paste0(path.wd, "../dataPrep/lsd/lsd.pa.centroids.rds"))

# synonyma (mám funkci v rgee)
# evaluace všech predikcí z ENMeval přes LSD PA pomocí sdm::evaluates(observed, predicted)


############
# settings
############

ndop.fs <- list("groups" = 37, "version" = "v1")

############
# execution
############

# sjednocení synonym LSD

lsd.pa.centroids.species <- unlist(unique(lsd.pa.centroids$TaxonNameLAT))
# načtení RDS Z ENMeval
rds_list <-
    list.files(
        path.models,
        pattern = paste0("^4ssos_"), # !!! skillc je samostatná dodělávka všech 116 druhů dohromady 20-22 variant
        ignore.case = TRUE,
        full.names = TRUE
    )
# rds_list <- readRDS("D:/PersonalWork/Balej/v2/RP/dataPrep/missing.rds")

# rozdělení RDS do více skupin pro paralelní zpracování
perGroup <- ceiling(length(rds_list) / ndop.fs$groups)
group.parts <- split(rds_list, ceiling(seq_along(rds_list) / perGroup))
# omezím původní seznam
rds_list <- group.parts[[cmd_arg]]
print(rds_list)

# cyklíme...
for (fpath in rds_list) {
    first <- TRUE
    out.t <- list()
    rds.r <- list()
    rds.l <- list()
    # dočasně rozdělené do dvou částí, nutno spojit:
    rds.l <- readRDS(fpath)
    bn <- basename(fpath)
    # rds.l2 <- readRDS(str_replace(fpath, "skill", "skillb")) # !!! zakomentovat

    sp <- names(rds.l)
    if (length(sp) > 1) {
        print("Více druhů v jednom souboru!")
        stop()
    }
    # nezávislá PA data z LSD pro daný druh
    lsd.temp <- lsd.pa.centroids %>% filter(TaxonNameLAT == sp)

    if (nrow(lsd.temp) == 0) {
        print("************************************************")
        print("není v LSD, nemůžu validovat:")
        print(sp)
        next
    } else {
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print(sp)
    }

    for (version in names(rds.l[[sp]])) {
        print("version: -------------------------")
        print(version)
        print("version: -------------------------")
        for (adjust in names(rds.l[[sp]][[version]])) {
            print("adjust: ____________")
            print(adjust)
            second <- TRUE



            # FC
            for (FC in names(rds.l[[sp]][[version]][[adjust]])) {
                # RM
                for (RM in names(rds.l[[sp]][[version]][[adjust]][[FC]])) {
                    if (is.logical(rds.l[[sp]][[version]][[adjust]][[FC]][[RM]])) {
                        # nékdy daný adjust nebyl upočitatelný a nic tam není...
                        next
                    }

                    for (layer in names(rds.l[[sp]][[version]][[adjust]][[FC]][[RM]]@predictions)) {
                        ev <- NA

                        print(layer)
                        r.temp <- rds.l[[sp]][[version]][[adjust]][[FC]][[RM]]@predictions[[layer]]
                        # rds.r[[sp]][[version]][[adjust]][[layer]] <- r.temp
                        ex.predicted <- extract(r.temp, st_coordinates(lsd.temp))

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
                    rds.l[[sp]][[version]][[adjust]][[FC]][[RM]]@results[["species"]] <- sp
                    rds.l[[sp]][[version]][[adjust]][[FC]][[RM]]@results[["version"]] <- version
                    rds.l[[sp]][[version]][[adjust]][[FC]][[RM]]@results[["adjust"]] <- adjust

                    temp.t <- rds.l[[sp]][[version]][[adjust]][[FC]][[RM]]@results %>% left_join(out.t.second, by = "tune.args")

                    occs <- rds.l[[sp]][[version]][[adjust]][[FC]][[RM]]@occs
                    temp.t$occs.n <- nrow(occs)

                    temp.t$bg.n <- nrow(rds.l[[sp]][[version]][[adjust]][[FC]][[RM]]@bg)
                    # LSD
                    temp.t$p.n <- nrow(lsd.temp %>% filter(presence == 1))
                    temp.t$a.n <- nrow(lsd.temp %>% filter(presence == 0))
                    # temp.t$p.wkt <- st_as_text(st_combine(lsd.temp %>% filter(presence == 1)))
                    # temp.t$a.wkt <- st_as_text(st_combine(lsd.temp %>% filter(presence == 0)))
                    # temp.t$occs.wkt <- st_as_text(st_combine(occs %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)))

                    if (first) {
                        first <- FALSE
                        out.t <- temp.t
                    } else {
                        out.t %<>% add_row(temp.t)
                    }
                    gc()
                }
            }
        }
    }
    # ukládat per species - nutno - jinak je výsledné rds s rastery predikcí extrémně velké a nenačitatelné/neuložitelné
    # saveRDS(out.t, paste0(path.eval, "ssos.t_", sp, "_", as.character(cmd_arg), ".rds")) # !!! nové názvy, nepřepsat původní
    # saveRDS(rds.r, paste0(path.eval, "ssos.r_", sp, "_", as.character(cmd_arg), ".rds")) # !!!

    saveRDS(out.t, paste0(path.eval, "t_", bn)) # !!! nové názvy, nepřepsat původní
    # neukládám, mám v původním rds s modely, bylo by duplicitní
    # saveRDS(rds.r, paste0(path.eval, "r_",  bn)) # !!!
}


# stop()
# vn <- versionNames()
# vn9 <- unlist(strsplit(vn[["9"]], "\\+"))
# vn[["9"]] <- vn9[1]
# vn[["10"]] <- vn9[2]

end_time <- Sys.time()
print(end_time - start_time)