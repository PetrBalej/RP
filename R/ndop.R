start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

# nastavit working directory
path.wd <- "C:/Users/petr/Downloads/RP/RP/"
setwd(path.wd)
path.data <- "C:/Users/petr/Downloads/RP/projects-data/"
path.rgee <- "C:/Users/petr/Documents/pc02/rgee20230227/" # samsung500ntfs # paste0(path.expand("~"), "/Downloads/rgee2/rgee")
# source(paste0(path.rgee, "R/export_raster/functions.R"))
path.wd.prep <- paste0(path.wd, "dataPrep/ndop/")



############
# inputs
############
ndop.csv.path <- paste0(path.data, "aopk/ndop/download20230307/")
sitmap_2rad.czechia <- readRDS(paste0(path.wd, "dataPrep/sitmap_2rad/sitmap_2rad-czechia.rds")) # 4326
SPH_STAT.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_STAT.shp"))

############
# settings
############
ndop.fs <- list("months" = c(4:6), "years" = c(2018:2022), "version" = "v1")


st_crs(SPH_STAT.source) <- 4326



ndop_ugc <-
  function(years_range = list(from = "2017-01-01", to = "2019-12-31"),
           season_months_range = list(from = 4, to = 7),
           import_path_ndop = "/../ndop/csv",
           res_crs = 5514,
           presicion = 100) {
    col_types <- cols_only(
      DATUM_DO = col_date("%Y%m%d"),
      DATUM_OD = col_date("%Y%m%d"),
      ZDROJ = "c", SITMAP = "i", X = "d", Y = "d", DAT_SADA = "c",
      CXPRESNOST = "d", PROJEKT = "c", KAT_TAX = "f", DRUH = "f",
      GARANCE = "c", VALIDACE = "c", ID_NALEZ = "n", NEGATIV = "i",
      VEROH = "i", AUTOR = "f"
    )

    # načte všechny *.csv z import_path_ndop
    csv_ndop <-
      list.files(
        path = import_path_ndop,
        pattern = "*.csv",
        full.names = T
      ) %>%
      map_df(~ read_csv(., col_types = cols(.default = "c"))) %>%
      filter(X != "<i>Skrytá lokalizace</i>") %>%
      filter(CXPRESNOST != "") %>%
      distinct(ID_NALEZ, .keep_all = TRUE) %>%
      # filter(DAT_SADA != "iNaturalist - data ČR") %>%
      # filter(AUTOR != "iNaturalist uživatel") %>%
      type_convert(
        col_types = col_types,
        locale = locale("cs", decimal_mark = ",")
      )


    # základní dofiltrovaní nálezů z NDOPu
    csv_ndop_filter <- csv_ndop %>%
      filter(
        DATUM_OD >= years_range$from & DATUM_OD <= years_range$to &
          DATUM_DO >= years_range$from & DATUM_DO <= years_range$to &
          between(
            month(DATUM_OD),
            season_months_range$from,
            season_months_range$to
          ) &
          between(
            month(DATUM_DO),
            season_months_range$from,
            season_months_range$to
          ) &
          (VEROH == 1 | VALIDACE == "věrohodný záznam") &
          NEGATIV == 0 &
          CXPRESNOST <= presicion
      )


    if (is.null(res_crs)) {
      res_crs <- 5514
    }

    if (res_crs != 5514) {
      csv_ndop_filter %<>%
        mutate(
          csv_ndop_filter %>% st_as_sf(coords = c("X", "Y"), crs = 5514) %>%
            st_transform(res_crs) %>%
            st_coordinates() %>% as_tibble()
        )
    }

    if (res_crs == 3035 || res_crs == 5514) {
      csv_ndop_filter$X %<>% as.integer
      csv_ndop_filter$Y %<>% as.integer
    } else {
      options(pillar.sigfig = 8) # jen pro případnou vizualizaci
    }

    # přidání sloupců s WGS84 souřadnicemi, výběr záznamů z polygonu a potřebných sloupců
    csv_ndop_filter %<>%
      dplyr::select(ID_NALEZ, DRUH, X, Y, KAT_TAX, CXPRESNOST, DATUM_OD, DATUM_DO, AUTOR)

    return(csv_ndop_filter)
  }


# zvážit jiný filtr garance?
ndop.ff <- ndop_ugc(
  years_range = list(from = "2018-01-01", to = "2022-12-31"),
  season_months_range = list(from = 4, to = 6),
  import_path_ndop = ndop.csv.path,
  res_crs = 5514,
  presicion = 1000
)


# ndop.ff$AUTOR %<>% as.factor
ndop.ff %<>% filter(KAT_TAX == "Ptáci")

# odstraním diakritiku a převedu na malá písmena
ndop.ff %<>% mutate(AUTOR = tolower(stri_trans_general(str = AUTOR, id = "Latin-ASCII")))
ndop.ff <- droplevels(ndop.ff)

au.base <- unique(as.vector(unlist(ndop.ff$AUTOR)))
length(au.base)

AUTOR.unique <- unique(trimws(unlist(strsplit(as.vector(unlist(au.base)), "[,]"))))
length(AUTOR.unique)

AUTOR.unique.one <- AUTOR.unique[stri_count_fixed(AUTOR.unique, " ") == 0]
print(AUTOR.unique.one)
# ponechat jen autory dvou- a víceslovné (u jednoslovných hrozí zgroupování nepříslušných autorů podle křestního jména)
AUTOR.unique.more <- AUTOR.unique[stri_count_fixed(AUTOR.unique, " ") > 0]



# počet nálezů na autora
ndop.ff.occs <- ndop.ff %>%
  group_by(AUTOR) %>%
  summarise(occs = n_distinct(ID_NALEZ)) %>%
  arrange(desc(occs))


## "lucie" ale nabere všechny lucie... - check pro jednoslovné - problematické
# nálezy per unique autor (mohou být duplicitně) - jde mi jen o druhovou znalost na autora
ndop.ff.au.first <- TRUE
for (au in AUTOR.unique.more) {
  print(au)
  ndop.ff.temp <- ndop.ff %>%
    filter(str_detect(AUTOR, au)) %>%
    mutate(autor.unique = au)

  if (ndop.ff.au.first) {
    ndop.ff.au.first <- FALSE
    ndop.ff.au <- ndop.ff.temp
  } else {
    ndop.ff.au %<>% add_row(ndop.ff.temp)
  }
}

# počet nálezů na autora (po rozdělení čárkou a dopárování v případě stejného jména)
ndop.ff.au.occs <- ndop.ff.au %>%
  group_by(autor.unique) %>%
  summarise(occs = n_distinct(ID_NALEZ)) %>%
  arrange(desc(occs))

# počet druhů na autora
ndop.ff.au.sps <- ndop.ff.au %>%
  group_by(autor.unique) %>%
  summarise(sps = n_distinct(DRUH)) %>%
  arrange(desc(sps))



# saveRDS(ndop.ff, paste0(path.wd.prep, "ndop.ff.rds"))
# saveRDS(ndop.ff.occs, paste0(path.wd.prep, "ndop.ff.occs.rds"))
# saveRDS(ndop.ff.au, paste0(path.wd.prep, "ndop.ff.au.rds"))
# saveRDS(ndop.ff.au.occs, paste0(path.wd.prep, "ndop.ff.au.occs.rds"))

boxplot(ndop.ff.au.occs$occs)
boxplot(ndop.ff.au.sps$sps)


# trvá dlouho
ndop.occsXsps <- ndop.ff.au.occs %>% left_join(ndop.ff.au.sps, by = "autor.unique")
# koreluje počet druhů s počtem nálezů (per autor?)
plot(sps ~ log(occs), data = ndop.occsXsps)



# per pixel

ndop.ff.au.sf <- ndop.ff.au %>% st_as_sf(coords = c("X", "Y"), crs = 5514)

ndop.ff.au.sf %<>% st_transform(st_crs(4326))

ndop.ff.au.sf.POLE <- sitmap_2rad.czechia %>% st_join(ndop.ff.au.sf)

ndop.ff.au.sf.POLE %<>% filter(!is.na(ID_NALEZ))

# saveRDS(ndop.ff.au.sf.POLE, paste0(path.wd.prep, "ndop.ff.au.sf.POLE.rds"))

# per pixel
ndop.ff.au.sf.POLE.pp <- ndop.ff.au.sf.POLE %>%
  group_by(DRUH, POLE, autor.unique) %>%
  slice_head(n = 1)


# per pixel
ndop.ff.au.sf.POLE.pp.c <- ndop.ff.au.sf.POLE.pp %>% st_centroid()

# saveRDS(ndop.ff.au.sf.POLE.pp, paste0(path.wd.prep, "ndop.ff.au.sf.POLE.pp.rds"))
# saveRDS(ndop.ff.au.sf.POLE.pp.c, paste0(path.wd.prep, "ndop.ff.au.sf.POLE.pp.c.rds"))

# ndop.ff.au.sf.POLE.pp.c  <- readRDS(paste0(path.wd.prep, "ndop.ff.au.sf.POLE.pp.c.rds"))
# ndop.ff.au.sf.POLE <- readRDS(paste0(path.wd.prep, "ndop.ff.au.sf.POLE.rds"))

SPH_STAT.source.s <- st_simplify(SPH_STAT.source, preserveTopology = FALSE, dTolerance = 1000)
druhy <- unique(as.vector(unlist(ndop.ff.au.sf.POLE$DRUH)))
first <- TRUE
for (druh in druhy) {
  print(druh)
  pres <- ndop.ff.au.sf.POLE.pp.c %>% filter(DRUH == druh)
  pres.au <- unique(as.vector(unlist(pres$autor.unique)))
  tgob <- ndop.ff.au.sf.POLE.pp.c %>% filter(autor.unique %in% pres.au)
  pres.unique <- pres %>%
    group_by(POLE) %>%
    slice_head(n = 1)
  tgob.unique <- tgob %>%
    group_by(POLE) %>%
    slice_head(n = 1)

  # počet (unikátních, zrušit znásobení spoluautory) nálezů na pixel
  tgob.count <- tgob %>%
    group_by(POLE) %>%
    mutate(pp = length(unique(ID_NALEZ))) %>%
    slice_head(n = 1)
  length(unique(tgob.count$POLE))


  df.temp <- data.frame(
    "species" = druh,
    "p" = nrow(pres), "pu" = nrow(pres.unique),
    "tgob" = nrow(tgob), "tgobu" = nrow(tgob.unique)
  )
  if (first) {
    first <- FALSE
    df.out <- df.temp
  } else {
    df.out %<>% add_row(df.temp)
  }

  pl <- tgob.count %>% dplyr::select(pp)

  ggpl <- ggplot() +
    geom_sf(data = SPH_STAT.source.s, fill = "#F0F0F0") +
    geom_sf(data = pl, aes(colour = pp), pch = 15, size = 2, ) +
    labs(color = "TGB (log10)", caption = "NDOP (2018-2022)") +
    scale_colour_gradient(low = "white", high = "red", trans = "log10") +
    theme_light() +
    geom_sf(data = pres.unique) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(druh, " \n", nrow(pres.unique), " / ", nrow(tgob.unique)))

  png(paste0(path.wd.prep, "tgob-", str_replace_all(druh, "[\\/\\:]", "_"), ".png"), width = 1500, height = 700)

  print(ggpl)
  dev.off()
}

# saveRDS(df.out, paste0(path.wd.prep, "df.out.rds"))
# write.csv(df.out, paste0(path.wd.prep, "df.out.csv"))

# uložení rozdělených dat per species pro model


plot(sitmap_2rad.czechia)

# vygenerovat i mapy presencí druhu (počet) a backgroundu podle odpovídajících autorů (počty i per pixel)



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