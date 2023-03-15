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
path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/projects-data/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.rgee <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/rgee/" # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
# source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "dataPrep/lsd/")

############
# inputs
############
lsd.source <- read_delim(paste0(path.data, "lsd/export-20230223.csv"), delim = ";", locale = locale(encoding = "Windows-1250")) #  ISO-8859-2

sitmap_2rad.czechia <- readRDS(paste0(path.wd, "dataPrep/sitmap_2rad/sitmap_2rad-czechia.rds"))
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

# filtrace problematických druhů (název)
lsd.filter %<>%
  filter(!str_detect(TaxonNameLAT, "×")) %>%
  filter(!str_detect(TaxonNameLAT, "/")) %>%
  filter(!str_detect(TaxonNameLAT, "sp\\.")) %>%
  filter(TaxonNameLAT != "Ondatra zibethica")
nrow(lsd.filter)


# potřebuju jen geometrii
SPH_STAT.filter <- SPH_STAT.source %<>% dplyr::select(-everything())
st_crs(SPH_STAT.filter) <- 4326
SPH_KRAJ.filter <- SPH_KRAJ.source %<>% dplyr::select(-everything())
st_crs(SPH_KRAJ.filter) <- 4326


#
# průběžný overview
#
sitmap_2rad.czechia.ow <- sitmap_2rad.czechia %>% filter(POLE %in% POLE.ObsListsID.freq$POLE)

png(paste0(path.wd.prep, "overview.png"), width = 1000, height = 800)
plot(SPH_KRAJ.filter, main = paste0(nrow(sitmap_2rad.czechia), " / ", nrow(sitmap_2rad.czechia.ow)))
par(new = TRUE)
plot(sitmap_2rad.czechia.ow %>% dplyr::select(-everything()), add = TRUE, col = "red")
dev.off()

st_write(SPH_STAT.filter, paste0(path.wd.prep, "overview-stat.shp"))
st_write(SPH_KRAJ.filter, paste0(path.wd.prep, "overview-kraj.shp"))
st_write(sitmap_2rad.czechia.ow, paste0(path.wd.prep, "overview-selected-2rad.shp"))

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

# cross check - odpovídají názvy kvadrátů ze SiteName reálně poloze (Lat, Lon)?
lsd.filter.anti <- lsd.filter %>% anti_join(sitmap_2rad.czechia, by = "POLE")
if (nrow(lsd.filter.anti > 0)) {
  print("Všechny kódy kvadrátů neodpovídají síti!")
  print(lsd.filter.anti)
}

lsd.filter.oneTaxonOccPerPOLE.simple <- lsd.filter.oneTaxonOccPerPOLE %>% dplyr::select(POLE, TaxonNameLAT)

#
# odvození absencí
#
"%notin%" <- Negate("%in%")
lsd.filter.oneTaxonOccPerPOLE.simple %<>% ungroup()
lsd.filter.oneTaxonOccPerPOLE.simple$presence <- 1

lsd.filter.species <- unique(lsd.filter.oneTaxonOccPerPOLE.simple$TaxonNameLAT)
lsd.pa <- lsd.filter.oneTaxonOccPerPOLE.simple
for (sp in lsd.filter.species) {
  sp.presences <- lsd.filter.oneTaxonOccPerPOLE.simple %>% filter(TaxonNameLAT == sp)
  sp.absences <- lsd.filter.oneTaxonOccPerPOLE.simple %>%
    filter(POLE %notin% sp.presences$POLE) %>%
    distinct(POLE)

  sp.absences$TaxonNameLAT <- sp
  sp.absences$presence <- 0

  lsd.pa %<>% add_row(sp.absences)
}

lsd.pa.polygons <- st_as_sf(lsd.pa %>% left_join(sitmap_2rad.czechia, by = "POLE"))
saveRDS(lsd.pa.polygons, paste0(path.wd.prep, "lsd.pa.polygons.rds"))
lsd.pa.centroids <- lsd.pa.polygons
lsd.pa.centroids <- st_as_sf(lsd.pa.centroids %<>% mutate(geometry = st_centroid(geometry)))
saveRDS(lsd.pa.centroids, paste0(path.wd.prep, "lsd.pa.centroids.rds"))

# + remove spatialy autocorrelated squares - geographically/environmentally
# alespoň reprezentativnost vůči altitude a landcoveru?

# do protokolu časy i adresu skriptu, který to vygeneroval + název vygenerovaného souboru

# https://github.com/PetrBalej/rgee/blob/master/R/clean/export_raster.R
# https://github.com/PetrBalej/rgee/blob/master/R/igaD.R

end_time <- Sys.time()
print(end_time - start_time)