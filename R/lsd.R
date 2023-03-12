start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

# nastavit working directory
path.wd <- "C:/Users/petr/Documents/pc02/projects/"
setwd(path.wd)
path.data <- "C:/Users/petr/Documents/pc02/projects-data/"
path.rgee <- "C:/Users/petr/Documents/pc02/rgee20230223/" # samsung500ntfs # paste0(path.expand("~"), "/Downloads/rgee2/rgee")
# source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "dataPrep/lsd/")

############
# inputs
############
lsd.source <- read_delim(paste0(path.data, "lsd/export-20230223.csv"), delim = ";")
sitmap_2rad.source <- st_read(paste0(path.data, "aopk/sitmap_2rad/sitmap_2rad.shp")) # 5514
SPH_STAT.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_STAT.shp"))
SPH_KRAJ.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_KRAJ.shp"))

############
# settings
############
lsd.fs <- list("minMin" = 50, "months" = c(4:6), "years" = c(2018:2022), "minSurveys" = 4, "version" = "v1")

############
# execution
############

#
# LSD
#

nrow(lsd.source)

# základní filtr
lsd.filter <- lsd.source %>%
  filter(grepl("^[0-9]{4}[a-d]{2}", SiteName) &
    !is.na(Lon) &
    !is.na(Lat) &
    UncertainIdent != 1) %>%
  mutate(POLE = substr(SiteName, start = 1, stop = 6))

nrow(lsd.filter)

# alespoň více než x minut (vyřešeny časy přes půlnoc do nového dne)
lsd.filter %<>% mutate(dt = ifelse((TimeEnd - TimeStart) / 60 > 0, (TimeEnd - TimeStart) / 60, (TimeEnd - TimeStart + 60 * 60 * 24) / 60)) %>% filter(dt >= lsd.fs[["minMin"]])

nrow(lsd.filter)

# jaro+roky
lsd.filter %<>% filter(Year %in% lsd.fs[["years"]] & Month %in% lsd.fs[["months"]])
nrow(lsd.filter)

# aspoň X návštěv na pole...
POLE.ObsListsID.freq <- lsd.filter %>%
  group_by(POLE) %>%
  summarise(visits = n_distinct(ObsListsID)) %>%
  # arrange(desc(visits)) %>%
  filter(visits >= lsd.fs[["minSurveys"]])

nrow(POLE.ObsListsID.freq)

lsd.filter %<>% filter(POLE %in% POLE.ObsListsID.freq$POLE)
nrow(lsd.filter)

#
# sitmap_2rad
#

# potřebuju jen POLE
sitmap_2rad.filter <- sitmap_2rad.source %>% dplyr::select(POLE)
sitmap_2rad.filter %<>% st_transform(st_crs(4326))


# potřebuju jen geometrii
SPH_STAT.filter <- SPH_STAT.source %<>% dplyr::select(-everything())
st_crs(SPH_STAT.filter) <- 4326
SPH_KRAJ.filter <- SPH_KRAJ.source %<>% dplyr::select(-everything())
st_crs(SPH_KRAJ.filter) <- 4326

# jen kvadráty v ČR
sitmap_2rad.filter %<>% filter(st_intersects(geometry, SPH_STAT.filter, sparse = FALSE))
nrow(sitmap_2rad.filter)



#
# průběžný overview
#
sitmap_2rad.filter.ow <- sitmap_2rad.filter %>% filter(POLE %in% POLE.ObsListsID.freq$POLE)

png(paste0(path.wd.prep, "overview.png"), width = 1000, height = 800)
plot(SPH_KRAJ.filter, main = paste0(nrow(sitmap_2rad.filter), " / ", nrow(sitmap_2rad.filter.ow)))
par(new = TRUE)
plot(sitmap_2rad.filter.ow %>% dplyr::select(-everything()), add = TRUE, col = "red")
dev.off()

st_write(SPH_STAT.filter, paste0(path.wd.prep, "overview-stat.shp"))
st_write(SPH_KRAJ.filter, paste0(path.wd.prep, "overview-kraj.shp"))
st_write(sitmap_2rad.filter.ow, paste0(path.wd.prep, "overview-selected-2rad.shp"))

saveRDS(lsd.filter, paste0(path.wd.prep, "overview-lsd.filter.rds"))


# Zuur, A.F., Ieno, E.N. and Elphick, C.S., 2010. A protocol for data exploration to avoid common statistical problems. Methods in ecology and evolution, 1(1), pp.3-14.
# - but a more stringent approach is to use values as low as 3 as we did here. High, or even moderate, collinearity is especially problematic when ecological signals are weak. In that case, even a VIF of 2 may cause nonsignificant parameter estimates, compared to the situation without collinearity.


lsd.filter.oneTaxonOccPerPOLE <- lsd.filter %>%
  group_by(POLE, TaxonNameLAT) %>%
  slice_head(n = 1)
# aspoň X návštěv na pole...
spCounts <- lsd.filter.oneTaxonOccPerPOLE %>%
  group_by(TaxonNameLAT) %>%
  summarise(spCount = n_distinct(ObsListsID)) %>%
  arrange(desc(spCount)) %>%
  filter(spCount >= 10)





stop()
# TODO: cross check - odpovídají názvy kvadrátů ze SiteName reálně poloze (Lat, Lon)? - udělat i prostorový průnik se sítí
lsd.filter %>% left_join(sitmap_2rad.source, by = "POLE")

st.sitmap_2rad.centroids <- st_centroid(st.sitmap_2rad)

# minimum x presenci a y absencí - možné až po získání prediktorů a "děr"

# + remove spatialy autocorrelated squares - geographically/environmentally
# alespoň reprezentativnost vůči altitude a landcoveru?

# do protokolu časy i adresu skriptu, který to vygeneroval + název vygenerovaného souboru

# https://github.com/PetrBalej/rgee/blob/master/R/clean/export_raster.R
# https://github.com/PetrBalej/rgee/blob/master/R/igaD.R

end_time <- Sys.time()
print(end_time - start_time)