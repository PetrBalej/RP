start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "lubridate") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

# nastavit working directory
path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/"
setwd(path.wd)
path.data <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/projects-data/"
path.rgee <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/rgee/" # samsung500ntfs # paste0(path.expand("~"), "/Downloads/rgee2/rgee")
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
ndop.fs <- list("months" = c(4:6), "years" = c(2019:2022), "precision" = 1000, "version" = "v1")


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
  years_range = list(from = paste0(min(ndop.fs$years), "-01-01"), to = paste0(max(ndop.fs$years), "-12-31")),
  season_months_range = list(from = min(ndop.fs$months), to = max(ndop.fs$months)),
  import_path_ndop = ndop.csv.path,
  res_crs = 5514,
  presicion = ndop.fs$precision
)


# ndop.ff$AUTOR %<>% as.factor
ndop.ff %<>% filter(KAT_TAX == "Ptáci")

# odstraním diakritiku a převedu na malá písmena
ndop.ff %<>% mutate(AUTOR = tolower(stri_trans_general(str = AUTOR, id = "Latin-ASCII")))
ndop.ff <- droplevels(ndop.ff)


au.base <- unique(as.vector(unlist(ndop.ff$AUTOR)))
length(au.base)
## nebudu dělat explode (čárkou) více spoluautorů - dohromady mohou mít jako tým jiné určovací schopnosti, které může jednotlivý autor postrádat a tím bych mohl individuálnímu autoru přiřadit druhy, které sám nepozná
# AUTOR.unique <- unique(trimws(unlist(strsplit(as.vector(unlist(au.base)), "[,]"))))
# length(AUTOR.unique)
au.base.one <- au.base[stri_count_fixed(au.base, " ") == 0]
print(au.base.one)
write.csv(au.base.one, file = paste0(path.wd.prep, "autors-check-single-remove.csv"), row.names = FALSE)
# ponechat jen autory dvou- a víceslovné (u jednoslovných hrozí zgroupování nepříslušných autorů podle křestního jména a nemusí být důvěryhodné)
au.base.more <- au.base[stri_count_fixed(au.base, " ") > 0]

# dofiltr podezřelých a problematických autorů
length(au.base.more)
write.csv(au.base.more, file = paste0(path.wd.prep, "autors-check.csv"), row.names = FALSE)

nrow(ndop.ff) # 425933
saveRDS(ndop.ff, paste0(path.wd.prep, "ndop.ff.rds"))


# pouze víceslovní autoři
ndop.ff.au.more <- ndop.ff %>% filter(AUTOR %in% au.base.more)
nrow(ndop.ff.au.more) # 425008



# počet nálezů na autora
ndop.ff.au.more.occs <- ndop.ff.au.more %>%
  group_by(AUTOR) %>%
  summarise(occs = n_distinct(ID_NALEZ)) %>%
  arrange(desc(occs))
saveRDS(ndop.ff.au.more.occs, paste0(path.wd.prep, "ndop.ff.au.more.occs.rds"))


# počet druhů na autora
ndop.ff.au.more.sps <- ndop.ff.au.more %>%
  group_by(AUTOR) %>%
  summarise(sps = n_distinct(DRUH)) %>%
  arrange(desc(sps))
saveRDS(ndop.ff.au.more.sps, paste0(path.wd.prep, "ndop.ff.au.more.sps.rds"))


boxplot(ndop.ff.au.more.occs$occs)
boxplot(ndop.ff.au.more.sps$sps)

ggplot(ndop.ff.au.more.occs, aes(y = occs)) +
  geom_boxplot() +
  scale_y_log10()
ggplot(ndop.ff.au.more.occs, aes(x = occs)) +
  geom_histogram() +
  scale_x_log10()


# koreluje počet druhů s počtem nálezů (per autor?)
ndop.occsXsps <- ndop.ff.au.more.occs %>% left_join(ndop.ff.au.more.sps, by = "AUTOR")

# Default scatter plot
sp <- ggplot(ndop.occsXsps, aes(x = occs, y = sps)) +
  geom_point()
sp + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10")



# souřadnice na SF
ndop.ff.au.more.sf <- ndop.ff.au.more %>% st_as_sf(coords = c("X", "Y"), crs = 5514)
ndop.ff.au.more.sf %<>% st_transform(st_crs(4326))

# namapování nálezů na POLE
ndop.ff.au.more.sf.POLE <- sitmap_2rad.czechia %>% st_join(ndop.ff.au.more.sf)
ndop.ff.au.more.sf.POLE %<>% filter(!is.na(ID_NALEZ))


saveRDS(ndop.ff.au.more.sf, paste0(path.wd.prep, "ndop.ff.au.more.sf.rds"))
saveRDS(ndop.ff.au.more.sf.POLE, paste0(path.wd.prep, "ndop.ff.au.more.sf.POLE.rds"))



# per pixel - unikátní hodnoty druhů na pixel a autora
ndop.ff.au.more.sf.POLE.pp <- ndop.ff.au.more.sf.POLE %>%
  group_by(DRUH, POLE, AUTOR) %>%
  slice_head(n = 1)
nrow(ndop.ff.au.more.sf.POLE.pp) # 210539
saveRDS(ndop.ff.au.more.sf.POLE.pp, paste0(path.wd.prep, "ndop.ff.au.more.sf.POLE.pp.rds"))

# jen body (centroidy pixelů)
ndop.ff.au.more.sf.POLE.pp.c <- ndop.ff.au.more.sf.POLE.pp %>% st_centroid()
saveRDS(ndop.ff.au.more.sf.POLE.pp.c, paste0(path.wd.prep, "ndop.ff.au.more.sf.POLE.pp.c.rds"))


# per pixel unikátní druhy
ndop.ff.au.more.sf.POLE.pp.sp <- ndop.ff.au.more.sf.POLE %>%
  group_by(DRUH, POLE) %>%
  slice_head(n = 1)
nrow(ndop.ff.au.more.sf.POLE.pp.sp) #
saveRDS(ndop.ff.au.more.sf.POLE.pp.sp, paste0(path.wd.prep, "ndop.ff.au.more.sf.POLE.pp.sp.rds"))



# filtrace problematických druhů (malý počet nálezů, název)
ndop.ff.au.more.sf.POLE.pp.sp.count <- ndop.ff.au.more.sf.POLE.pp.sp %>%
  st_drop_geometry() %>%
  group_by(DRUH) %>%
  count(DRUH) %>%
  arrange(desc(n))
saveRDS(ndop.ff.au.more.sf.POLE.pp.sp.count, paste0(path.wd.prep, "ndop.ff.au.more.sf.POLE.pp.sp.count.rds"))

ndop.sp.selected <- ndop.ff.au.more.sf.POLE.pp.sp.count %>%
  filter(n >= 10) %>%
  filter(!str_detect(DRUH, "×")) %>%
  filter(!str_detect(DRUH, "/")) %>%
  filter(!str_detect(DRUH, "sp\\."))
saveRDS(ndop.sp.selected, paste0(path.wd.prep, "ndop.sp.selected.rds"))


#
# vizualizace map TGOB (per pixel occs) a presencí
#
SPH_STAT.source.s <- st_simplify(SPH_STAT.source, preserveTopology = FALSE, dTolerance = 0.001)
# druhy <- unique(as.vector(unlist(ndop.ff.au.more.sf.POLE$DRUH)))
druhy <- unique(as.vector(unlist(ndop.sp.selected$DRUH)))
first <- TRUE
for (druh in druhy) {
  print(druh)
  pres <- ndop.ff.au.more.sf.POLE.pp.c %>% filter(DRUH == druh)
  pres.au <- unique(as.vector(unlist(pres$AUTOR)))
  tgob <- ndop.ff.au.more.sf.POLE.pp.c %>% filter(AUTOR %in% pres.au)
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
    geom_sf(data = SPH_STAT.source.s, fill = "LightBlue") +
    geom_sf(data = pl, aes(colour = pp), pch = 15, size = 2) +
    labs(color = "TGOB (log10)", caption = paste0("occs source: AOPK NDOP (years: ", min(ndop.fs$years), "-", max(ndop.fs$years), "; months: ", min(ndop.fs$months), "-", max(ndop.fs$months), ") | author: Petr Balej | WGS84 | grid: KFME (SitMap_2Rad) | generated: ", Sys.Date())) +
    scale_colour_gradient(low = "white", high = "red", trans = "log10") +
    theme_light() +
    geom_sf(data = pres.unique) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    ggtitle(
      label = (druh),
      subtitle = paste0("TGOB (same species observers occs per pixel) | ", nrow(pres.unique), " / ", nrow(tgob.unique))
    )

  png(paste0(path.wd.prep, "tgob/tgob-", str_replace_all(druh, "[\\/\\:]", "_"), ".png"), width = 1200, height = 700)

  print(ggpl)
  dev.off()
}


saveRDS(df.out, paste0(path.wd.prep, "df.out.rds"))
write.csv(df.out, paste0(path.wd.prep, "df.out.csv"))




# uložení rozdělených dat per species pro model


plot(sitmap_2rad.czechia)

# vygenerovat i mapy presencí druhu (počet) a backgroundu podle odpovídajících autorů (počty i per pixel)



############
# execution
############




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