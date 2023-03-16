start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "geosphere", "ENMeval", "sdm") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "dataPrep/ndopTGOBeval/")
path.wd.prep.ndop <- paste0(path.wd, "dataPrep/ndop/")
source(paste0(path.wd, "shared.R"))
path.wd.prep.tgobEval <- paste0(path.wd, "dataPrep/ndopTGOBeval/")
path.temp.res <- "/home/petr/Downloads/delete/paper1/ndopTGOB - kopie/"
source(paste0(path.wd, "shared.R"))
############
# inputs
############

lsd.pa.centroids <- readRDS(paste0(path.wd, "dataPrep/lsd/lsd.pa.centroids.rds"))

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
        pattern = paste0("^skill_"),
        ignore.case = TRUE,
        full.names = TRUE
    )

# cyklíme...
for (fpath in rds_list) {
    rds.r <- list()
    # zabránit zbytečnému načtení RDS, pokud druh není v LSD
    fname <- unlist(strsplit(basename(fpath), "_")) # na druhé pozici je druh
    species.syn <- synonyms_unite(as.data.frame(list("species" = fname[2])))
    if (sum(lsd.pa.centroids.syn.species %in% species.syn[1, 1]) < 1) {
        print("Není v LSD:")
        print(fname[2])
        next
    }

    rds.l <- readRDS(fpath)
    sp <- names(rds.l)

    if (length(sp) > 1) {
        print("Více druhů v jednom souboru!")
        stop()
    }

    lsd.temp <- lsd.pa.centroids.syn %>% filter(species == sp)

    print("++++++++++++++++++++++++++++++++++++++++++++++++")
    print(sp)

    for (version in names(rds.l[[sp]])) {
        print("version: -------------------------")
        print(version)
        for (adjust in names(rds.l[[sp]][[version]])) {
            print(adjust)
            for (duplORnot in names(rds.l[[sp]][[version]][[adjust]])) {
                print(duplORnot)
                for (layer in names(rds.l[[sp]][[version]][[adjust]][[duplORnot]]@predictions)) {
                    print(layer)
                    r.temp <- rds.l[[sp]][[version]][[adjust]][[duplORnot]]@predictions[[layer]]
                    ex.predicted <- extract(r.temp, st_coordinates(lsd.temp))

                    if (inherits(try(
                        ev <- sdm::evaluates(lsd.temp$presence, ex.predicted)
                    ), "try-error")) {
                        print("E3error")
                        ev <- NA
                    }
                    rds.r[[sp]][[version]][[adjust]][[duplORnot]][[layer]] <- ev
                }
                gc()
            }
        }
    }
    saveRDS(rds.r, paste0(path.wd.prep, sp, ".rds"))
}

stop()
vn <- versionNames()
vn9 <- unlist(strsplit(vn[["9"]], "\\+"))
vn[["9"]] <- vn9[1]
vn[["10"]] <- vn9[2]


end_time <- Sys.time()
print(end_time - start_time)

print(sp.group)
