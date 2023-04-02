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

path.wd <- paste0(gcfl(), "/")
# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- paste0(path.wd, "../../projects-data/")
path.prep <- paste0(path.wd, "../../dataPrep/")
path.rgee <- paste0(path.wd, "../../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.lsd <- paste0(path.prep, "../dataPrep/lsd/") # hodinovka / lsd


############
# inputs
############
# hodinovka/export-20230326.csv lsd/export-20230223.csv
lsd.source <- read_delim(paste0(path.data, "lsd/export-20230223.csv"), delim = ";", locale = locale(encoding = "Windows-1250")) #  ISO-8859-2

sitmap_2rad.czechia <- readRDS(paste0(path.prep, "sitmap_2rad/sitmap_2rad-czechia.rds"))
SPH_STAT.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_STAT.shp"))
SPH_KRAJ.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_KRAJ.shp"))

############
# settings
############
lsd.fs <- list(
  "minMin" = 50, "months" = c(4:6), "years" = c(2019:2022),
  "minSurveys" = 4, "minPA" = 10, "version" = "v1",
  "filterProblematicSpecies" = c("×", "/", "sp\\.", "f\\. domestica")
)

############
# execution
############

"%notin%" <- Negate("%in%")

#
# LSD
#

lsd.source <- droplevels(lsd.source)
nrow(lsd.source)

# základní filtr
lsd <- lsd.source %>%
  filter(grepl("^[0-9]{4}[a-d]{2}", SiteName) &
    # !is.na(Lon) &
    # !is.na(Lat) &
    UncertainIdent != 1) %>%
  mutate(POLE = substr(SiteName, start = 1, stop = 6))

nrow(lsd)

# alespoň více než x minut (vyřešeny časy přes půlnoc do nového dne)
lsd %<>% mutate(dt = ifelse((TimeEnd - TimeStart) / 60 > 0, (TimeEnd - TimeStart) / 60, (TimeEnd - TimeStart + 60 * 60 * 24) / 60)) %>% filter(dt >= lsd.fs[["minMin"]])

nrow(lsd)

# jaro+roky
lsd %<>% filter(Year %in% lsd.fs[["years"]] & Month %in% lsd.fs[["months"]])
nrow(lsd)
length(unique(lsd$TaxonNameLAT))

#
# odstraňuji problematické druhy, pro vstup jako presence jednotlivých druhů do SDM
#
lsd %<>% filter(TaxonNameLAT != "Ondatra zibethica")
if (length(lsd.fs$filterProblematicSpecies) > 0) {
  # filtrace problematických druhů (název)
  for (bad in lsd.fs$filterProblematicSpecies) {
    lsd %<>%
      filter(!str_detect(TaxonNameLAT, bad))
  }
  print("LSD po odstranění problematických druhů:")
  nrow(lsd)
  length(unique(lsd$TaxonNameLAT))
}

#
# sjednocení poddruhů na druhy (nerozlišuju, nemělo by mít význam z hlediska nároku různých poddruhů?)
#
lsd %<>% rowwise() %>% mutate(TaxonNameLAT = paste(strsplit(TaxonNameLAT, " ")[[1]][1:2], collapse = " "))
print("NDOP po downgrade poddruhů:")
nrow(lsd)
length(unique(lsd$TaxonNameLAT))

#
# synonymizace (sjednocení různé taxonomie)
#
lsd <- synonyms_unite(lsd, spCol = "TaxonNameLAT")
print("NDOP po synonymizaci poddruhů:")
nrow(lsd)
length(unique(lsd$TaxonNameLAT))

bad <- unname(unlist(nepuvodni_problematicke()))
lsd %<>% filter(TaxonNameLAT %notin% bad)
nrow(lsd)
length(unique(lsd$TaxonNameLAT))


# prostorový join (nelze věřit názvu co je v Sitename?)
lsd.coords <- lsd %>%
  filter(!is.na(ObsLon) & !is.na(ObsLat)) %>%
  st_as_sf(coords = c("ObsLon", "ObsLat"), crs = 4326)

lsd.coords %<>% st_join(sitmap_2rad.czechia, suffix = c("", "_2rad")) # získám prostorově odvozený název kvadrátu, ne ten na pevno zadaný

lsd.coords <- as_tibble(lsd.coords) # odstraním sf

lsd.coords.ne <- lsd.coords %>% filter(POLE != POLE_2rad)

# st_write(lsd.coords.ne, "delete-hodinovky-nesoulad.shp")

if (nrow(lsd.coords.ne) > 0) {
  # kde je pravda, v souřadnicích nebo SiteName???
  # - větší chyby jsou v překlepech v názvech kvadrátů (často mimo i 10 km pokud je překlep v 1řádu), jinak jde často jen o rozdíl začátku a konce trasy, kdy lidi přecházejí v rámci hodinovky přes dva kvadráty 2 řádu...
  # Jak to pro Atlas řešil T. Telenský??? Dát echo správci AVIFu... - chyby jsou ale propsané i do NDOP???
  print("manuálně zadané kódy kvadrátů neodpovídají prostorovým! Použiju odvozené prostorem... ")
  print(lsd.coords.ne)
}

# join čistě podle názvu
lsd.codes <- lsd %>%
  filter(is.na(ObsLon) & is.na(ObsLat)) %>%
  dplyr::select(-c(ObsLat, ObsLon))
lsd.codes.join <- lsd.codes %>% left_join(sitmap_2rad.czechia, by = "POLE")
# cross check - odpovídají názvy kvadrátů ze SiteName reálně poloze (Lat, Lon)?
lsd.codes.anti <- lsd.codes %>% anti_join(sitmap_2rad.czechia, by = "POLE")
if (nrow(lsd.codes.anti) > 0) {
  print("Všechny kódy kvadrátů neodpovídají síti! Použiju pouze ty, které se protly (nejsou alespoň mimo ČR)")
  print(lsd.codes.anti)
}
lsd.codes.join$POLE_2rad <- NA
lsd.codes.join$POLE_2rad %<>% as.character

# spojím obě verze
lsd.codes.join %<>% add_row(lsd.coords)
lsd <- lsd.codes.join
lsd %<>% rename(POLEpuv = POLE) %>% mutate(POLE = ifelse(is.na(POLE_2rad), POLEpuv, POLE_2rad))
nrow(lsd)


# aspoň X návštěv na pole...
POLE.ObsListsID.freq <- lsd %>%
  group_by(POLE) %>%
  summarise(visits = n_distinct(ObsListsID)) %>%
  # arrange(desc(visits)) %>%
  filter(visits >= lsd.fs$minSurveys)
nrow(POLE.ObsListsID.freq)


lsd %<>% filter(POLE %in% POLE.ObsListsID.freq$POLE)
nrow(lsd)


# potřebuju jen geometrii
SPH_STAT.filter <- SPH_STAT.source %<>% dplyr::select(-everything())
st_crs(SPH_STAT.filter) <- 4326
SPH_KRAJ.filter <- SPH_KRAJ.source %<>% dplyr::select(-everything())
st_crs(SPH_KRAJ.filter) <- 4326
st_write(SPH_STAT.filter, paste0(path.lsd, "overview-stat.shp"))
st_write(SPH_KRAJ.filter, paste0(path.lsd, "overview-kraj.shp"))

#
# průběžný overview
#
sitmap_2rad.czechia.ow <- sitmap_2rad.czechia %>%
  filter(POLE %in% POLE.ObsListsID.freq$POLE) %>%
  distinct()

png(paste0(path.lsd, "overview.png"), width = 1000, height = 800)
plot(SPH_KRAJ.filter, main = paste0(nrow(sitmap_2rad.czechia), " / ", nrow(sitmap_2rad.czechia.ow)))
par(new = TRUE)
plot(sitmap_2rad.czechia.ow %>% dplyr::select(-everything()), add = TRUE, col = "red")
dev.off()
# st_write(sitmap_2rad.czechia.ow, paste0(path.lsd, "overview-selected-2rad.shp"))
# saveRDS(lsd, paste0(path.lsd, "overview-lsd.rds"))

#
# thinning blízkých čtverců (odstranění potenciální spatial autocorrelation)
#
sitmap_2rad.czechia.ow.df <- as.data.frame(st_coordinates(st_centroid(sitmap_2rad.czechia.ow)))[, 1:2]
names(sitmap_2rad.czechia.ow.df) <- c("x", "y")
lsd.thinned <- ecospat.occ.desaggregation(xy = sitmap_2rad.czechia.ow.df, min.dist = ((1 / 6) / 4) * 1.5, by = NULL)

lsd.thinned %<>% st_as_sf(coords = c("x", "y"), crs = 4326)
lsd.thinned$thin <- 1
sitmap_2rad.czechia.ow %<>% st_join(lsd.thinned, suffix = c("", "_thin")) %>% filter(thin == 1)

png(paste0(path.lsd, "overview-thin.png"), width = 1000, height = 800)
plot(SPH_KRAJ.filter, main = paste0(nrow(sitmap_2rad.czechia), " / ", nrow(sitmap_2rad.czechia.ow)))
par(new = TRUE)
plot(sitmap_2rad.czechia.ow %>% dplyr::select(-everything()), add = TRUE, col = "red")
dev.off()
st_write(sitmap_2rad.czechia.ow, paste0(path.lsd, "sitmap_2rad.czechia.ow.shp"))
saveRDS(sitmap_2rad.czechia.ow, paste0(path.lsd, "sitmap_2rad.czechia.ow.rds"))

lsd %<>% filter(POLE %in% sitmap_2rad.czechia.ow$POLE)
lsd <- st_as_sf(as_tibble(lsd) %>% dplyr::select(-geometry) %>% left_join(sitmap_2rad.czechia, by = "POLE"))
nrow(lsd)

saveRDS(lsd, paste0(path.lsd, "lsd.rds"))
st_write(lsd, paste0(path.lsd, "lsd.shp"))


# Zuur, A.F., Ieno, E.N. and Elphick, C.S., 2010. A protocol for data exploration to avoid common statistical problems. Methods in ecology and evolution, 1(1), pp.3-14.
# - but a more stringent approach is to use values as low as 3 as we did here. High, or even moderate, collinearity is especially problematic when ecological signals are weak. In that case, even a VIF of 2 may cause nonsignificant parameter estimates, compared to the situation without collinearity.


# per POLE a druh
lsd.oneTaxonOccPerPOLE.simple <- lsd %>%
  group_by(POLE, TaxonNameLAT) %>%
  slice_head(n = 1) %>%
  dplyr::select(POLE, TaxonNameLAT)

#
# odvození absencí
#
"%notin%" <- Negate("%in%")
lsd.oneTaxonOccPerPOLE.simple %<>% ungroup()
lsd.oneTaxonOccPerPOLE.simple$presence <- 1

lsd.species <- unique(lsd.oneTaxonOccPerPOLE.simple$TaxonNameLAT)
lsd.pa <- lsd.oneTaxonOccPerPOLE.simple
for (sp in lsd.species) {
  print(sp)
  sp.presences <- lsd.oneTaxonOccPerPOLE.simple %>% filter(TaxonNameLAT == sp)
  sp.absences <- lsd.oneTaxonOccPerPOLE.simple %>%
    filter(POLE %notin% sp.presences$POLE) %>%
    group_by(POLE) %>%
    slice_head(n = 1)

  sp.absences$TaxonNameLAT <- sp
  sp.absences$presence <- 0

  lsd.pa %<>% add_row(sp.absences)
}


lsd.pa.polygons <- st_as_sf(lsd.pa)
saveRDS(lsd.pa.polygons, paste0(path.lsd, "lsd.pa.polygons.rds"))
lsd.pa.centroids <- lsd.pa.polygons
lsd.pa.centroids <- st_as_sf(st_centroid(lsd.pa.polygons))
saveRDS(lsd.pa.centroids, paste0(path.lsd, "lsd.pa.centroids.rds"))

# alespoň 10 presencí nebo 10 absencí
lsd.pa.min <- as_tibble(lsd.pa.centroids) %>%
  group_by(TaxonNameLAT) %>%
  summarize(np = sum(presence), na = sum(presence == 0)) %>%
  filter(np >= lsd.fs$minPA & na >= lsd.fs$minPA)
saveRDS(lsd.pa.min, paste0(path.lsd, "lsd.pa.min.rds"))

# + remove spatialy autocorrelated squares - geographically/environmentally
# alespoň reprezentativnost vůči altitude a landcoveru?

# do protokolu časy i adresu skriptu, který to vygeneroval + název vygenerovaného souboru

# https://github.com/PetrBalej/rgee/blob/master/R/clean/export_raster.R
# https://github.com/PetrBalej/rgee/blob/master/R/igaD.R

end_time <- Sys.time()
print(end_time - start_time)