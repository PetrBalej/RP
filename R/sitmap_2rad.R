start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

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

# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/"
setwd(path.wd)
path.data <- "../projects-data/"
path.rgee <- "../../rgee/" # samsung500ntfs # paste0(path.expand("~"), "/Downloads/rgee2/rgee")
# source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "dataPrep/sitmap_2rad/")



############
# inputs
############
sitmap_2rad.source <- st_read(paste0(path.data, "aopk/sitmap_2rad/sitmap_2rad.shp")) # 5514
SPH_STAT.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_STAT.shp"))
SPH_KRAJ.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_KRAJ.shp"))


############
# execution
############

#
# sitmap_2rad
#

# potřebuju jen POLE
sitmap_2rad.filter <- sitmap_2rad.source %>% dplyr::select(POLE)
sitmap_2rad.filter %<>% st_transform(st_crs(4326))


# potřebuju jen geometrii
SPH_STAT.filter <- SPH_STAT.source %<>% dplyr::select(-everything())
st_crs(SPH_STAT.filter) <- 4326

# jen kvadráty v ČR
sitmap_2rad.filter %<>% filter(st_intersects(geometry, SPH_STAT.filter, sparse = FALSE))
nrow(sitmap_2rad.filter)

#
# overview
#

png(paste0(path.wd.prep, "overview.png"), width = 1000, height = 800)
plot(sitmap_2rad.filter %>% dplyr::select(-everything()), main = nrow(sitmap_2rad.filter))
dev.off()

st_write(sitmap_2rad.filter, paste0(path.wd.prep, "sitmap_2rad-czechia.shp"))
saveRDS(sitmap_2rad.filter, paste0(path.wd.prep, "sitmap_2rad-czechia.rds"))


end_time <- Sys.time()
print(end_time - start_time)