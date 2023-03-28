start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "ecospat") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- "../projects-data/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.rgee <- "../../rgee/" # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "../dataPrep/lsd/") # hodinovka / lsd

############
# inputs
############
# hodinovka/export-20230326.csv lsd/export-20230223.csv
lsd.source <- read_delim(paste0(path.data, "lsd/export-20230223.csv"), delim = ";", locale = locale(encoding = "Windows-1250")) #  ISO-8859-2

sitmap_2rad.czechia <- readRDS(paste0(path.wd, "dataPrep/sitmap_2rad/sitmap_2rad-czechia.rds"))
SPH_STAT.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_STAT.shp"))
SPH_KRAJ.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_KRAJ.shp"))

############
# settings
############
lsd.fs <- list("minMin" = 50, "months" = c(4:6), "years" = c(2019:2022), "minSurveys" = 4, "minPA" = 10, "version" = "v1")

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
    # !is.na(Lon) &
    # !is.na(Lat) &
    UncertainIdent != 1) %>%
  mutate(POLE = substr(SiteName, start = 1, stop = 6))

nrow(lsd.filter)

# alespoň více než x minut (vyřešeny časy přes půlnoc do nového dne)
lsd.filter %<>% mutate(dt = ifelse((TimeEnd - TimeStart) / 60 > 0, (TimeEnd - TimeStart) / 60, (TimeEnd - TimeStart + 60 * 60 * 24) / 60)) %>% filter(dt >= lsd.fs[["minMin"]])

nrow(lsd.filter)

# jaro+roky
lsd.filter %<>% filter(Year %in% lsd.fs[["years"]] & Month %in% lsd.fs[["months"]])
nrow(lsd.filter)

# filtrace problematických druhů (název)
lsd.filter %<>%
  filter(!str_detect(TaxonNameLAT, "×")) %>%
  filter(!str_detect(TaxonNameLAT, "/")) %>%
  filter(!str_detect(TaxonNameLAT, "sp\\.")) %>%
  filter(!str_detect(TaxonNameLAT, "f\\. domestica")) %>%
  filter(TaxonNameLAT != "Ondatra zibethica")
nrow(lsd.filter)

# sjednocení poddruhů na druhy (nerozlišuju, nemělo by mít význam z hlediska nároku různých poddruhů?)
lsd.filter %<>% rowwise() %>% mutate(TaxonNameLAT = paste(strsplit(TaxonNameLAT, " ")[[1]][1:2], collapse = " "))
lsd.filter <- droplevels(lsd.filter)


# prostorový join (nelze věřit názvu ceo je v Sitename?)
lsd.filter.coords <- lsd.filter %>%
  filter(!is.na(ObsLon) & !is.na(ObsLat)) %>%
  st_as_sf(coords = c("ObsLon", "ObsLat"), crs = 4326)

lsd.filter.coords %<>% st_join(sitmap_2rad.czechia, suffix = c("", "_2rad")) # získám prostorově odvozený název kvadrátu, ne ten na pevno zadaný

lsd.filter.coords <- as_tibble(lsd.filter.coords) # odstraním sf

lsd.filter.coords.ne <- lsd.filter.coords %>% filter(POLE != POLE_2rad)

# st_write(lsd.filter.coords.ne, "delete-hodinovky-nesoulad.shp")

if (nrow(lsd.filter.coords.ne) > 0) {
  # kde je pravda, v souřadnicích nebo SiteName???
  # - větší chyby jsou v překlepech v názvech kvadrátů (často mimo i 10 km pokud je překlep v 1řádu), jinak jde často jen o rozdíl začátku a konce trasy, kdy lidi přecházejí v rámci hodinovky přes dva kvadráty 2 řádu...
  # Jak to pro Atlas řešil T. Telenský??? Dát echo správci AVIFu... - chyby jsou ale propsané i do NDOP???
  print("manuálně zadané kódy kvadrátů neodpovídají prostorovým! Použiju odvozené prostorem... ")
  print(lsd.filter.coords.ne)
}

# join čistě podle názvu
lsd.filter.codes <- lsd.filter %>%
  filter(is.na(ObsLon) & is.na(ObsLat)) %>%
  dplyr::select(-c(ObsLat, ObsLon))
lsd.filter.codes.join <- lsd.filter.codes %>% left_join(sitmap_2rad.czechia, by = "POLE")
# cross check - odpovídají názvy kvadrátů ze SiteName reálně poloze (Lat, Lon)?
lsd.filter.codes.anti <- lsd.filter.codes %>% anti_join(sitmap_2rad.czechia, by = "POLE")
if (nrow(lsd.filter.codes.anti) > 0) {
  print("Všechny kódy kvadrátů neodpovídají síti! Použiju pouze ty, které se protly (nejsou alespoň mimo ČR)")
  print(lsd.filter.codes.anti)
}
lsd.filter.codes.join$POLE_2rad <- NA
lsd.filter.codes.join$POLE_2rad %<>% as.character

# spojím obě verze
lsd.filter.codes.join %<>% add_row(lsd.filter.coords)
lsd.filter <- lsd.filter.codes.join
lsd.filter %<>% rename(POLEpuv = POLE) %>% mutate(POLE = ifelse(is.na(POLE_2rad), POLEpuv, POLE_2rad))
nrow(lsd.filter)


# aspoň X návštěv na pole...
POLE.ObsListsID.freq <- lsd.filter %>%
  group_by(POLE) %>%
  summarise(visits = n_distinct(ObsListsID)) %>%
  # arrange(desc(visits)) %>%
  filter(visits >= lsd.fs$minSurveys)
nrow(POLE.ObsListsID.freq)


lsd.filter %<>% filter(POLE %in% POLE.ObsListsID.freq$POLE)
nrow(lsd.filter)


# potřebuju jen geometrii
SPH_STAT.filter <- SPH_STAT.source %<>% dplyr::select(-everything())
st_crs(SPH_STAT.filter) <- 4326
SPH_KRAJ.filter <- SPH_KRAJ.source %<>% dplyr::select(-everything())
st_crs(SPH_KRAJ.filter) <- 4326
st_write(SPH_STAT.filter, paste0(path.wd.prep, "overview-stat.shp"))
st_write(SPH_KRAJ.filter, paste0(path.wd.prep, "overview-kraj.shp"))

#
# průběžný overview
#
sitmap_2rad.czechia.ow <- sitmap_2rad.czechia %>%
  filter(POLE %in% POLE.ObsListsID.freq$POLE) %>%
  distinct()

png(paste0(path.wd.prep, "overview.png"), width = 1000, height = 800)
plot(SPH_KRAJ.filter, main = paste0(nrow(sitmap_2rad.czechia), " / ", nrow(sitmap_2rad.czechia.ow)))
par(new = TRUE)
plot(sitmap_2rad.czechia.ow %>% dplyr::select(-everything()), add = TRUE, col = "red")
dev.off()
# st_write(sitmap_2rad.czechia.ow, paste0(path.wd.prep, "overview-selected-2rad.shp"))
# saveRDS(lsd.filter, paste0(path.wd.prep, "overview-lsd.filter.rds"))


sitmap_2rad.czechia.ow.df <- as.data.frame(st_coordinates(st_centroid(sitmap_2rad.czechia.ow)))[, 1:2]
names(sitmap_2rad.czechia.ow.df) <- c("x", "y")
lsd.thinned <- ecospat.occ.desaggregation(xy = sitmap_2rad.czechia.ow.df, min.dist = ((1 / 6) / 4) * 1.5, by = NULL)

lsd.thinned %<>% st_as_sf(coords = c("x", "y"), crs = 4326)
lsd.thinned$thin <- 1
sitmap_2rad.czechia.ow %<>% st_join(lsd.thinned, suffix = c("", "_thin")) %>% filter(thin == 1)

png(paste0(path.wd.prep, "overview-thin.png"), width = 1000, height = 800)
plot(SPH_KRAJ.filter, main = paste0(nrow(sitmap_2rad.czechia), " / ", nrow(sitmap_2rad.czechia.ow)))
par(new = TRUE)
plot(sitmap_2rad.czechia.ow %>% dplyr::select(-everything()), add = TRUE, col = "red")
dev.off()
st_write(sitmap_2rad.czechia.ow, paste0(path.wd.prep, "sitmap_2rad.czechia.ow.shp"))
saveRDS(sitmap_2rad.czechia.ow, paste0(path.wd.prep, "sitmap_2rad.czechia.ow.rds"))

lsd.filter %<>% filter(POLE %in% sitmap_2rad.czechia.ow$POLE)
lsd.filter <- st_as_sf(as_tibble(lsd.filter) %>% dplyr::select(-geometry) %>% left_join(sitmap_2rad.czechia, by = "POLE"))
nrow(lsd.filter)

saveRDS(lsd.filter, paste0(path.wd.prep, "lsd.filter.rds"))
st_write(lsd.filter, paste0(path.wd.prep, "lsd.filter.shp"))


# Zuur, A.F., Ieno, E.N. and Elphick, C.S., 2010. A protocol for data exploration to avoid common statistical problems. Methods in ecology and evolution, 1(1), pp.3-14.
# - but a more stringent approach is to use values as low as 3 as we did here. High, or even moderate, collinearity is especially problematic when ecological signals are weak. In that case, even a VIF of 2 may cause nonsignificant parameter estimates, compared to the situation without collinearity.


# per POLE a druh
lsd.filter.oneTaxonOccPerPOLE.simple <- lsd.filter %>%
  group_by(POLE, TaxonNameLAT) %>%
  slice_head(n = 1) %>%
  dplyr::select(POLE, TaxonNameLAT)

#
# odvození absencí
#
"%notin%" <- Negate("%in%")
lsd.filter.oneTaxonOccPerPOLE.simple %<>% ungroup()
lsd.filter.oneTaxonOccPerPOLE.simple$presence <- 1

lsd.filter.species <- unique(lsd.filter.oneTaxonOccPerPOLE.simple$TaxonNameLAT)
lsd.pa <- lsd.filter.oneTaxonOccPerPOLE.simple
for (sp in lsd.filter.species) {
  print(sp)
  sp.presences <- lsd.filter.oneTaxonOccPerPOLE.simple %>% filter(TaxonNameLAT == sp)
  sp.absences <- lsd.filter.oneTaxonOccPerPOLE.simple %>%
    filter(POLE %notin% sp.presences$POLE) %>%
    group_by(POLE) %>%
    slice_head(n = 1)

  sp.absences$TaxonNameLAT <- sp
  sp.absences$presence <- 0

  lsd.pa %<>% add_row(sp.absences)
}


lsd.pa.polygons <- st_as_sf(lsd.pa)
saveRDS(lsd.pa.polygons, paste0(path.wd.prep, "lsd.pa.polygons.rds"))
lsd.pa.centroids <- lsd.pa.polygons
lsd.pa.centroids <- st_as_sf(st_centroid(lsd.pa.polygons))
saveRDS(lsd.pa.centroids, paste0(path.wd.prep, "lsd.pa.centroids.rds"))

# alespoň 10 presencí nebo 10 absencí
lsd.pa.min <- as_tibble(lsd.pa.centroids) %>%
  group_by(TaxonNameLAT) %>%
  summarize(np = sum(presence), na = sum(presence == 0)) %>%
  filter(np >= lsd.fs$minPA & na >= lsd.fs$minPA)
saveRDS(lsd.pa.min, paste0(path.wd.prep, "lsd.pa.min.rds"))

# + remove spatialy autocorrelated squares - geographically/environmentally
# alespoň reprezentativnost vůči altitude a landcoveru?

# do protokolu časy i adresu skriptu, který to vygeneroval + název vygenerovaného souboru

# https://github.com/PetrBalej/rgee/blob/master/R/clean/export_raster.R
# https://github.com/PetrBalej/rgee/blob/master/R/igaD.R

end_time <- Sys.time()
print(end_time - start_time)