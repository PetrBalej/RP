start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "lubridate") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/"
setwd(path.wd)
path.data <- paste0(path.wd, "../../projects-data/")
path.prep <- paste0(path.wd, "../../dataPrep/")
path.rgee <- paste0(path.wd, "../../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.ndop <- paste0(path.prep, "ndop/")


############
# inputs
############
ndop.csv.path <- paste0(path.data, "aopk/ndop/download20230307/") # download20230328 - 2014-2017: v řadě krajů (~RP) není garantována ani validována značná část nálezů - nepoužitelné!!!
sitmap_2rad.czechia <- readRDS(paste0(path.prep, "sitmap_2rad/sitmap_2rad-czechia.rds")) # 4326
SPH_STAT.source <- st_read(paste0(path.data, "cuzk/SPH_SHP_WGS84/WGS84/SPH_STAT.shp"))
lsd.pa.min <- readRDS(paste0(path.prep, "lsd/lsd.pa.min.rds"))
############
# settings
############
ndop.fs <- list(
  "months" = c(4:6), "years" = c(2019:2022), "precision" = 1000,
  "version" = "v1",
  "authorsFilter" = c(0),
  "filterProblematicSpecies" = c("×", "/", "sp\\.", "f\\. domestica", "sensu lato"),
  "minPres" = 10 # aplikován pouze pro generování mapek SSOS
)


st_crs(SPH_STAT.source) <- 4326



ndop_ugc <-
  function(years_range = list(from = "2017-01-01", to = "2019-12-31"),
           season_months_range = list(from = 4, to = 7),
           import_path_ndop = "/../ndop/csv",
           res_crs = 5514,
           presicion = 100,
           removeLSD = FALSE) {
    col_types <- cols_only(
      DATUM_DO = col_date("%Y%m%d"),
      DATUM_OD = col_date("%Y%m%d"),
      ZDROJ = "c", SITMAP = "i", X = "d", Y = "d", DAT_SADA = "c",
      CXPRESNOST = "d", PROJEKT = "c", KAT_TAX = "c", DRUH = "c",
      GARANCE = "c", VALIDACE = "c", ID_NALEZ = "n", NEGATIV = "i",
      VEROH = "i", AUTOR = "c",
      UMIST_NAL = "c"
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

    if (removeLSD) {
      # čtverec 2rad má cca 3x3 km (reálně trochu méně) polovina úhlopříčky je sqrt(3^2+3^2)/2 = 2.1 (bezpečnější je spodní hranice cca 1.8)
      # přesnost (CXPRESNOST) <1800m tedy odstraní LSD a Hodinovkové datasety AVIFu automaticky!!!
      # do budoucna na to ale nelze spoléhat, pokud bych z nějakého důvodu uznal za dostatečné přesnosti větší než 1.8 km!!!
      print("LSD dataset odstraněn")
      csv_ndop %<>% filter(str_detect(UMIST_NAL, "Liniové sčítání druhů") %in% c(FALSE, NA)) # str_detect nad NA vrací NA, ne T/F
      # lze dofiltrovat i Hodinovky z NDOP?
    }

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

############
# execution
############

# zvážit jiný filtr garance?
ndop <- ndop_ugc(
  years_range = list(from = paste0(min(ndop.fs$years), "-01-01"), to = paste0(max(ndop.fs$years), "-12-31")),
  season_months_range = list(from = min(ndop.fs$months), to = max(ndop.fs$months)),
  import_path_ndop = ndop.csv.path,
  res_crs = 5514,
  presicion = ndop.fs$precision
)

ndop.orig <- ndop
# ndop <- ndop.orig
# ndop$AUTOR %<>% as.factor
ndop %<>% filter(KAT_TAX == "Ptáci")


# sjednocuji (ať zbytečně neodstraňuji...)
ndop %<>% mutate(DRUH = ifelse(DRUH == "Acanthis flammea/cabaret", "Acanthis flammea", DRUH))


# záloha původních autorských jmen pro kontrolu
ndop %<>% mutate(AUTOR.orig = AUTOR)

# odstraním diakritiku a převedu na malá písmena
ndop %<>% mutate(AUTOR = tolower(stri_trans_general(str = AUTOR, id = "Latin-ASCII")))
ndop %<>% mutate(AUTOR = str_replace(AUTOR, "\\.", "\\. "))
ndop %<>% mutate(AUTOR = trimws(AUTOR))


au.base <- unique(trimws(as.vector(unlist(ndop$AUTOR))))
length(au.base)
## nebudu dělat explode (čárkou) více spoluautorů - dohromady mohou mít jako tým jiné určovací schopnosti, které může jednotlivý autor postrádat a tím bych mohl individuálnímu autoru přiřadit druhy, které sám nepozná
# AUTOR.unique <- unique(trimws(unlist(strsplit(as.vector(unlist(au.base)), "[,]"))))
# length(AUTOR.unique)
au.base.one <- au.base[stri_count_fixed(au.base, " ") == 0]
print(au.base.one)
write.csv(au.base.one, file = paste0(path.ndop, "autors-check-single-remove.csv"), row.names = FALSE)
# ponechat jen autory dvou- a víceslovné (u jednoslovných hrozí zgroupování nepříslušných autorů podle křestního jména a nemusí být důvěryhodné)
au.base.more <- au.base[stri_count_fixed(au.base, " ") > 0]

# dofiltr podezřelých a problematických autorů
au.bad <- au.base.more[str_detect(au.base.more, "anonym|krouzkovaci stanice")]
print(au.bad)
write.csv(au.bad, file = paste0(path.ndop, "autors-check-bad-remove.csv"), row.names = FALSE)

write.csv(au.base.more, file = paste0(path.ndop, "autors-check-good.csv"), row.names = FALSE)


ndop %<>% mutate(AUTOR0 = ifelse(stri_count_fixed(AUTOR, " ") == 0, 1, ifelse(str_detect(AUTOR, "anonym|krouzkovaci stanice"), 2, 0)))



# souřadnice na SF
ndop %<>% st_as_sf(coords = c("X", "Y"), crs = 5514)
ndop %<>% st_transform(st_crs(4326))
ndop %<>% mutate("X" = st_coordinates(geometry)[, 1], "Y" = st_coordinates(geometry)[, 2])


# namapování nálezů na POLE
ndop.POLE <- sitmap_2rad.czechia %>% st_join(ndop)
ndop.POLE %<>% filter(!is.na(ID_NALEZ))
ndop.POLE.join <- as_tibble(ndop.POLE %>% dplyr::select(POLE, ID_NALEZ)) %>% dplyr::select(-geometry)

ndop %<>% left_join(ndop.POLE.join, by = "ID_NALEZ")

ndop %<>% filter(AUTOR0 %in% ndop.fs$authorsFilter) # dálě pracuju jen s neproblematickými názvy autorů
ndop.POLE %<>% filter(AUTOR0 %in% ndop.fs$authorsFilter) # dálě pracuju jen s neproblematickými názvy autorů
print("NDOP základ:")
nrow(ndop) # 425933
nrow(ndop.POLE)
length(unique(ndop$DRUH))

#
# pro celkové TGOB úsilí (kernel smoothing) má význam mít vše
#
saveRDS(ndop, paste0(path.ndop, "ndop.rds"))
saveRDS(ndop.POLE, paste0(path.ndop, "ndop.POLE.rds"))

# záloha
ndop.tgob <- ndop
ndop.POLE.tgob <- ndop.POLE
#######################################################################################################################################################################################

#
# odstraňuji problematické druhy, pro vstup jako presence jednotlivých druhů do SDM a SSOS
#
ndop %<>% filter(DRUH != "Luscinia svecica svecica")
ndop.POLE %<>% filter(DRUH != "Luscinia svecica svecica")
if (length(ndop.fs$filterProblematicSpecies) > 0) {
  # filtrace problematických druhů (název)
  for (bad in ndop.fs$filterProblematicSpecies) {
    ndop %<>%
      filter(!str_detect(DRUH, bad))
    ndop.POLE %<>%
      filter(!str_detect(DRUH, bad))
  }
  print("NDOP po odstranění problematických druhů:")
  nrow(ndop)
  nrow(ndop.POLE)
  length(unique(ndop$DRUH))
}


#
# sjednocení poddruhů na druhy (nerozlišuju, nemělo by mít význam z hlediska nároku různých poddruhů?)
#
ndop %<>% rowwise() %>% mutate(DRUH = paste(strsplit(DRUH, " ")[[1]][1:2], collapse = " "))
ndop.POLE %<>% rowwise() %>% mutate(DRUH = paste(strsplit(DRUH, " ")[[1]][1:2], collapse = " "))
print("NDOP po downgrade poddruhů:")
nrow(ndop)
nrow(ndop.POLE)
length(unique(ndop$DRUH))

#
# synonymizace (sjednocení různé taxonomie)
#
ndop <- synonyms_unite(ndop, spCol = "DRUH") %>% st_as_sf()
ndop.POLE <- synonyms_unite(ndop.POLE, spCol = "DRUH") %>% st_as_sf()

print("NDOP po synonymizaci poddruhů:")
nrow(ndop)
nrow(ndop.POLE)
length(unique(ndop$DRUH))


# # omezení: modeluju jen druhy které můžu z LSD ověřit (v budoucnu za dalších podmínek klidně všechny...) - v základu nechávám všechny (abych měl top druhy), doomezuju až při modelování
# ndop %<>% filter(DRUH %in% unique(lsd.pa.min$TaxonNameLAT))
# ndop.POLE %<>% filter(DRUH %in% unique(lsd.pa.min$TaxonNameLAT))

saveRDS(ndop, paste0(path.ndop, "ndopP.rds"))
saveRDS(ndop.POLE, paste0(path.ndop, "ndopP.POLE.rds"))

#####################
# základní statistiky
#####################
ndop.stat <- as_tibble(ndop) %>% dplyr::select(-geometry) # pro výpočty je rychlejší udělat z toho normální tibble bez geometrie (sfc)
ID_NALEZ.total <- nrow(ndop.stat)
DRUH.total <- length(unique(unlist((ndop.stat$DRUH))))
POLE.total <- length(unique(unlist((ndop.stat$POLE))))
AUTOR.total <- length(unique(unlist((ndop.stat$AUTOR))))

### podle autora
ndop.stat.res <- ndop.stat %>%
  group_by(AUTOR) %>%
  summarise(
    ID_NALEZ.n = n_distinct(ID_NALEZ),
    ID_NALEZ.pct = n_distinct(ID_NALEZ) / ID_NALEZ.total * 100,
    DRUH.n = n_distinct(DRUH),
    DRUH.pct = n_distinct(DRUH) / DRUH.total * 100,
    POLE.n = n_distinct(POLE),
    POLE.pct = n_distinct(POLE) / POLE.total * 100
  ) %>%
  arrange(desc(ID_NALEZ.n))

saveRDS(ndop.stat.res, paste0(path.ndop, "ndop.stat.res.rds"))
write.csv(ndop.stat.res, file = paste0(path.ndop, "ndop.stat.res.csv"), row.names = FALSE)

### podle druhu
ndop.stat.res.sp <- ndop.stat %>%
  group_by(DRUH) %>%
  summarise(
    ID_NALEZ.n = n_distinct(ID_NALEZ),
    ID_NALEZ.pct = n_distinct(ID_NALEZ) / ID_NALEZ.total * 100,
    AUTOR.n = n_distinct(AUTOR),
    AUTOR.pct = n_distinct(AUTOR) / AUTOR.total * 100,
    POLE.n = n_distinct(POLE),
    POLE.pct = n_distinct(POLE) / POLE.total * 100
  ) %>%
  arrange(desc(ID_NALEZ.n))

saveRDS(ndop.stat.res.sp, paste0(path.ndop, "ndop.stat.res.sp.rds"))
write.csv(ndop.stat.res.sp, file = paste0(path.ndop, "ndop.stat.res.sp.csv"), row.names = FALSE)



#
# generování přehledových grafů (AUTOR)
#
options(scipen = 999)
ndop.stat.res.boxplots <- ndop.stat.res %>% dplyr::select(-AUTOR)

# nálezy
ndop.stat.res.boxplots %<>% arrange(desc(ID_NALEZ.n))
png(paste0(path.ndop, "AUTOR-observers-ID_NALEZ.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(ndop.stat.res.boxplots$ID_NALEZ.n, log = "x", main = "observers contribution to NDOP by ID_NALEZ", xlab = "observers ordered by activity (most -> less ID_NALEZ) [log10]", ylab = "sum ID_NALEZ")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(ndop.stat.res.boxplots$ID_NALEZ.pct, log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(ndop.stat.res.boxplots$ID_NALEZ.pct))) # sekundární y osa
mtext("sum ID_NALEZ [%]", side = 4, line = 3)
dev.off()

# nálezy cumsum
ndop.stat.res.boxplots %<>% arrange(desc(ID_NALEZ.n))
png(paste0(path.ndop, "AUTOR-observers-ID_NALEZ-cumsum.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(cumsum(ndop.stat.res.boxplots$ID_NALEZ.n), log = "x", main = "observers contribution to NDOP by ID_NALEZ", xlab = "observers ordered by activity (most -> less ID_NALEZ) [log10]", ylab = "cumsum ID_NALEZ")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(cumsum(ndop.stat.res.boxplots$ID_NALEZ.pct), log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(cumsum(ndop.stat.res.boxplots$ID_NALEZ.pct)))) # sekundární y osa
mtext("cumsum ID_NALEZ [%]", side = 4, line = 3)
dev.off()


# nej autoři
ndop.stat.res.topAUTOR <- ndop.stat.res %>%
  arrange(desc(ID_NALEZ.n)) %>%
  slice_head(n = 100)
saveRDS(ndop.stat.res.topAUTOR, paste0(path.ndop, "ndop.stat.res.topAUTOR.rds"))
write.csv(ndop.stat.res.topAUTOR, paste0(path.ndop, "ndop.stat.res.topAUTOR.csv"), row.names = FALSE)

# 100
ndop.topAUTOR100 <- ndop %>% filter(AUTOR %in% unique(ndop.stat.res.topAUTOR$AUTOR))
ndop.POLE.topAUTOR100 <- ndop.POLE %>% filter(AUTOR %in% unique(ndop.stat.res.topAUTOR$AUTOR))
saveRDS(ndop.topAUTOR100, paste0(path.ndop, "ndop.topAUTOR100.rds"))
saveRDS(ndop.POLE.topAUTOR100, paste0(path.ndop, "ndop.POLE.topAUTOR100.rds"))

# 10
ndop.stat.res.topAUTOR %<>% slice_head(n = 10)
ndop.topAUTOR10 <- ndop %>% filter(AUTOR %in% unique(ndop.stat.res.topAUTOR$AUTOR))
ndop.POLE.topAUTOR10 <- ndop.POLE %>% filter(AUTOR %in% unique(ndop.stat.res.topAUTOR$AUTOR))
saveRDS(ndop.topAUTOR10, paste0(path.ndop, "ndop.topAUTOR10.rds"))
saveRDS(ndop.POLE.topAUTOR10, paste0(path.ndop, "ndop.POLE.topAUTOR10.rds"))

# druhy (0)
png(paste0(path.ndop, "AUTOR-observers-DRUH-ordered-ID_NALEZ.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(ndop.stat.res.boxplots$DRUH.n, log = "x", main = "observers contribution to NDOP by DRUH", xlab = "observers ordered by activity (most -> less ID_NALEZ) [log10]", ylab = "sum DRUH")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(ndop.stat.res.boxplots$DRUH.pct, log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(ndop.stat.res.boxplots$DRUH.pct))) # sekundární y osa
mtext("sum DRUH [%]", side = 4, line = 3)
dev.off()

# druhy
ndop.stat.res.boxplots %<>% arrange(desc(DRUH.n))
png(paste0(path.ndop, "AUTOR-observers-DRUH.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(ndop.stat.res.boxplots$DRUH.n, log = "x", main = "observers contribution to NDOP by DRUH", xlab = "observers ordered by activity (most -> less DRUH) [log10]", ylab = "sum DRUH")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(ndop.stat.res.boxplots$DRUH.pct, log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(ndop.stat.res.boxplots$DRUH.pct))) # sekundární y osa
mtext("sum DRUH [%]", side = 4, line = 3)
dev.off()

# pole
ndop.stat.res.boxplots %<>% arrange(desc(POLE.n))
png(paste0(path.ndop, "AUTOR-observers-POLE.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(ndop.stat.res.boxplots$POLE.n, log = "x", main = "observers contribution to NDOP by POLE", xlab = "observers ordered by activity (most -> less POLE) [log10]", ylab = "sum POLE")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(ndop.stat.res.boxplots$POLE.pct, log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(ndop.stat.res.boxplots$POLE.pct))) # sekundární y osa
mtext("sum POLE [%]", side = 4, line = 3)
dev.off()



#
# generování přehledových grafů (DRUH) (! šlo by v cyklu hromadně...)
#
options(scipen = 999)
ndop.stat.res.boxplots.sp <- ndop.stat.res.sp %>% dplyr::select(-DRUH)

# nálezy
ndop.stat.res.boxplots.sp %<>% arrange(desc(ID_NALEZ.n))
png(paste0(path.ndop, "DRUH-observers-ID_NALEZ.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(ndop.stat.res.boxplots.sp$ID_NALEZ.n, log = "x", main = "species contribution to NDOP by ID_NALEZ", xlab = "species ordered by activity (most -> less ID_NALEZ) [log10]", ylab = "sum ID_NALEZ")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(ndop.stat.res.boxplots.sp$ID_NALEZ.pct, log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(ndop.stat.res.boxplots.sp$ID_NALEZ.pct))) # sekundární y osa
mtext("sum ID_NALEZ [%]", side = 4, line = 3)
dev.off()

# nálezy cumsum
ndop.stat.res.boxplots.sp %<>% arrange(desc(ID_NALEZ.n))
png(paste0(path.ndop, "DRUH-observers-ID_NALEZ-cumsum.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(cumsum(ndop.stat.res.boxplots.sp$ID_NALEZ.n), log = "x", main = "species contribution to NDOP by ID_NALEZ", xlab = "species ordered by activity (most -> less ID_NALEZ) [log10]", ylab = "cumsum ID_NALEZ")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(cumsum(ndop.stat.res.boxplots.sp$ID_NALEZ.pct), log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(cumsum(ndop.stat.res.boxplots.sp$ID_NALEZ.pct)))) # sekundární y osa
mtext("cumsum ID_NALEZ [%]", side = 4, line = 3)
dev.off()

# nej druhy
ndop.stat.res.sp.topDRUH <- ndop.stat.res.sp %>%
  arrange(desc(ID_NALEZ.n)) %>%
  slice_head(n = 50)
saveRDS(ndop.stat.res.sp.topDRUH, paste0(path.ndop, "ndop.stat.res.sp.topDRUH.rds"))
write.csv(ndop.stat.res.sp.topDRUH, paste0(path.ndop, "ndop.stat.res.sp.topDRUH.csv"), row.names = FALSE)

# 50
ndop.topDRUH50 <- ndop %>% filter(DRUH %in% unique(ndop.stat.res.sp.topDRUH$DRUH))
ndop.POLE.topDRUH50 <- ndop.POLE %>% filter(DRUH %in% unique(ndop.stat.res.sp.topDRUH$DRUH))
saveRDS(ndop.topDRUH50, paste0(path.ndop, "ndop.topDRUH50.rds"))
saveRDS(ndop.POLE.topDRUH50, paste0(path.ndop, "ndop.POLE.topDRUH50.rds"))

# 10
ndop.stat.res.sp.topDRUH %<>% slice_head(n = 10)
ndop.topDRUH10 <- ndop %>% filter(DRUH %in% unique(ndop.stat.res.sp.topDRUH$DRUH))
ndop.POLE.topDRUH10 <- ndop.POLE %>% filter(DRUH %in% unique(ndop.stat.res.sp.topDRUH$DRUH))
saveRDS(ndop.topDRUH10, paste0(path.ndop, "ndop.topDRUH10.rds"))
saveRDS(ndop.POLE.topDRUH10, paste0(path.ndop, "ndop.POLE.topDRUH10.rds"))

# autoři (0)
png(paste0(path.ndop, "DRUH-observers-AUTOR-ordered-ID_NALEZ.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(ndop.stat.res.boxplots.sp$AUTOR.n, log = "x", main = "species contribution to NDOP by AUTOR", xlab = "species ordered by activity (most -> less ID_NALEZ) [log10]", ylab = "sum AUTOR")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(ndop.stat.res.boxplots.sp$AUTOR.pct, log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(ndop.stat.res.boxplots.sp$AUTOR.pct))) # sekundární y osa
mtext("sum AUTOR [%]", side = 4, line = 3)
dev.off()

# autoři
ndop.stat.res.boxplots.sp %<>% arrange(desc(AUTOR.n))
png(paste0(path.ndop, "DRUH-observers-AUTOR.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(ndop.stat.res.boxplots.sp$AUTOR.n, log = "x", main = "species contribution to NDOP by DRUH", xlab = "species ordered by activity (most -> less AUTOR) [log10]", ylab = "sum AUTOR")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(ndop.stat.res.boxplots.sp$AUTOR.pct, log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(ndop.stat.res.boxplots.sp$AUTOR.pct))) # sekundární y osa
mtext("sum AUTOR [%]", side = 4, line = 3)
dev.off()

# pole
ndop.stat.res.boxplots.sp %<>% arrange(desc(POLE.n))
png(paste0(path.ndop, "DRUH-observers-POLE.png"), width = 520, height = 500)
par(mar = c(4, 4, 2, 4))
plot(ndop.stat.res.boxplots.sp$POLE.n, log = "x", main = "species contribution to NDOP by POLE", xlab = "species ordered by activity (most -> less POLE) [log10]", ylab = "sum POLE")
grid(NULL, NULL, lty = 6)
par(new = TRUE)
plot(ndop.stat.res.boxplots.sp$POLE.pct, log = "x", axes = FALSE, xlab = "", ylab = "") # negeneruju popisky os
grid(NULL, NULL, lty = 6)
axis(side = 4, at = pretty(range(ndop.stat.res.boxplots.sp$POLE.pct))) # sekundární y osa
mtext("sum POLE [%]", side = 4, line = 3)
dev.off()


#########################

# koreluje počet druhů s počtem nálezů (per autor?)
png(paste0(path.ndop, "species-occurrences.png"), width = 520, height = 500)
ndop.stat.res.boxplots %<>% arrange(desc(ID_NALEZ.n))
plot(ndop.stat.res.boxplots$DRUH.n ~ ndop.stat.res.boxplots$ID_NALEZ.n, log = "x", main = "species ~ occurrences (per observer)", xlab = "occurrences  [log10]", ylab = "species")
grid(NULL, NULL, lty = 6)
dev.off()


#
# přehledové mapy
#

# počet druhů na pole
ndop.POLE.join.sp <- as_tibble(ndop.POLE) %>% dplyr::select(POLE, DRUH, -geometry)
ndop.POLE.join.sp %<>%
  group_by(POLE) %>% summarise(
    DRUH.n = n_distinct(DRUH)
  )
ndop.POLE.join.sp %<>% left_join(sitmap_2rad.czechia, by = "POLE") %>% ungroup()
png(paste0(path.ndop, "map-species-per-POLE.png"), width = 1000)
plot(st_as_sf(ndop.POLE.join.sp) %>% dplyr::select(DRUH.n), breaks = unname(quantile(unique(as.vector(unlist(as_tibble(ndop.POLE.join.sp) %>% dplyr::select(DRUH.n, -geometry)))))))
dev.off()

# počet nálezů na pole
ndop.POLE.join.occs <- as_tibble(ndop.POLE) %>% dplyr::select(POLE, ID_NALEZ, -geometry)
ndop.POLE.join.occs %<>%
  group_by(POLE) %>% summarise(
    ID_NALEZ.n = n_distinct(ID_NALEZ)
  )
ndop.POLE.join.occs %<>% left_join(sitmap_2rad.czechia, by = "POLE") %>% ungroup()
png(paste0(path.ndop, "map-occs-per-POLE.png"), width = 1000)
plot(st_as_sf(ndop.POLE.join.occs) %>% dplyr::select(ID_NALEZ.n), breaks = unname(quantile(unique(as.vector(unlist(as_tibble(ndop.POLE.join.occs) %>% dplyr::select(ID_NALEZ.n, -geometry)))))))
dev.off()

# počet autorů na pole
ndop.POLE.join.au <- as_tibble(ndop.POLE) %>% dplyr::select(POLE, AUTOR, -geometry)
ndop.POLE.join.au %<>%
  group_by(POLE) %>% summarise(
    AUTOR.n = n_distinct(AUTOR)
  )
ndop.POLE.join.au %<>% left_join(sitmap_2rad.czechia, by = "POLE") %>% ungroup()
png(paste0(path.ndop, "map-observers-per-POLE.png"), width = 1000)
plot(st_as_sf(ndop.POLE.join.au) %>% dplyr::select(AUTOR.n), breaks = unname(quantile(unique(as.vector(unlist(as_tibble(ndop.POLE.join.au) %>% dplyr::select(AUTOR.n, -geometry)))))))
dev.off()







#
# vizualizace map SSOS-TGOB (per pixel occs) a presencí
#

dir.create(paste0(path.ndop, "ssos"), showWarnings = FALSE)

SPH_STAT.source.s <- st_simplify(SPH_STAT.source, preserveTopology = FALSE, dTolerance = 0.001)
# druhy <- unique(as.vector(unlist(ndop.sf.POLE$DRUH)))
druhy <- unique(as.vector(unlist(ndop$DRUH)))
first <- TRUE
first2 <- TRUE
for (druh in druhy) {
  pres <- ndop.POLE %>% filter(DRUH == druh)
  pres.au <- unique(as.vector(unlist(pres$AUTOR)))
  ssos <- ndop.POLE %>% filter(AUTOR %in% pres.au)
  pres.unique <- pres %>%
    group_by(POLE) %>%
    slice_head(n = 1) %>%
    st_centroid()

  if (nrow(pres.unique) < ndop.fs$minPres) {
    # negeneruju zbytečně mapky s druhy s malým počtem presencí...
    print("přeskakuju (limit minPres):")
    print(druh)
    next
  }
  print(druh)

  ssos.unique <- ssos %>%
    group_by(POLE) %>%
    slice_head(n = 1)

  # počet (unikátních, zrušit znásobení spoluautory) nálezů na pixel
  ssos.count <- ssos %>%
    group_by(POLE) %>%
    mutate(pp = length(unique(ID_NALEZ))) %>%
    slice_head(n = 1) %>%
    st_centroid()


  df.temp <- data.frame(
    "species" = druh,
    "p" = nrow(pres), "pu" = nrow(pres.unique),
    "ssos" = nrow(ssos), "ssosu" = nrow(ssos.unique),
    "authors" = length(pres.au)
  )
  if (first) {
    first <- FALSE
    df.out <- df.temp
  } else {
    df.out %<>% add_row(df.temp)
  }

  pl <- ssos.count %>% dplyr::select(pp)

  ggpl <- ggplot() +
    geom_sf(data = SPH_STAT.source.s, fill = "LightBlue") +
    geom_sf(data = pl, aes(colour = pp), pch = 15, size = 2) +
    labs(color = "SSOS-TGOB (log10)", caption = paste0("occs source: AOPK NDOP (years: ", min(ndop.fs$years), "-", max(ndop.fs$years), "; months: ", min(ndop.fs$months), "-", max(ndop.fs$months), ") | author: Petr Balej | WGS84 | grid: KFME (SitMap_2Rad) | generated: ", Sys.Date())) +
    scale_colour_gradient(low = "white", high = "red", trans = "log10") +
    theme_light() +
    geom_sf(data = pres.unique) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    ggtitle(
      label = (druh),
      subtitle = paste0("SSOS-TGOB (same species observers occs per pixel) | p:", df.temp$pu, " / bg:", df.temp$ssosu, " / autors:", df.temp$authors)
    )

  png(paste0(path.ndop, "ssos/", str_replace_all(druh, "[\\/\\:]", "_"), ".png"), width = 1200, height = 700)

  print(ggpl)
  dev.off()




  #
  # ssos - omezení autorů s "nepřirozeným"" poměrem sběru daného druhu
  #

  ### celkové
  ssos.DRUH <- as_tibble(ssos) %>%
    dplyr::select(-geometry) %>%
    group_by(DRUH) %>%
    summarise(POLE.n = n_distinct(POLE))
  ssos.DRUH.sum.all <- sum(ssos.DRUH$POLE.n)
  ssos.DRUH.sum <- as.numeric(ssos.DRUH %>% filter(DRUH == druh) %>% dplyr::select(POLE.n))
  # "přirozený" poměr
  ssos.DRUH.ratio <- (ssos.DRUH.sum / ssos.DRUH.sum.all) * 0.5 # snížím hranici na 30%, odstraní extrémy (lidi co druh značí jen občas)

  ### autorské
  ssos.DRUH.AUTOR <- as_tibble(ssos) %>%
    dplyr::select(-geometry) %>%
    group_by(DRUH, AUTOR) %>%
    summarise(POLE.n = n_distinct(POLE))
  ssos.DRUH.AUTOR.sum.all <- ssos.DRUH.AUTOR %>%
    ungroup() %>%
    group_by(AUTOR) %>%
    summarise(POLE.n.sum.all = sum(POLE.n))
  ssos.DRUH.AUTOR.sum <- ssos.DRUH.AUTOR %>%
    ungroup() %>%
    filter(DRUH == druh) %>%
    group_by(AUTOR) %>%
    summarise(POLE.n.sum = sum(POLE.n))

  ssos.DRUH.AUTOR.ratio <- ssos.DRUH.AUTOR.sum %>%
    left_join(ssos.DRUH.AUTOR.sum.all, by = "AUTOR") %>%
    mutate(ratio = POLE.n.sum / POLE.n.sum.all) %>%
    mutate(ratioLess = ifelse(ratio < ssos.DRUH.ratio, 1, 0))

  # unname(quantile(ssos.DRUH.AUTOR.ratio$ratio, probs=0.25)


  ssos.natural.au <- ssos.DRUH.AUTOR.ratio %>% filter(ratioLess == 0)

  pres.au <- unique(ssos.natural.au$AUTOR)

  ssos <- ndop.POLE %>% filter(AUTOR %in% pres.au)

  ### opakování

  ssos.unique <- ssos %>%
    group_by(POLE) %>%
    slice_head(n = 1)

  # počet (unikátních, zrušit znásobení spoluautory) nálezů na pixel
  ssos.count <- ssos %>%
    group_by(POLE) %>%
    mutate(pp = length(unique(ID_NALEZ))) %>%
    slice_head(n = 1) %>%
    st_centroid()


  df.temp <- data.frame(
    "species" = druh,
    "p" = nrow(pres), "pu" = nrow(pres.unique),
    "ssos" = nrow(ssos), "ssosu" = nrow(ssos.unique),
    "authors" = length(pres.au)
  )
  if (first2) {
    first2 <- FALSE
    df.out2 <- df.temp
  } else {
    df.out2 %<>% add_row(df.temp)
  }

  pl <- ssos.count %>% dplyr::select(pp)

  ggpl <- ggplot() +
    geom_sf(data = SPH_STAT.source.s, fill = "LightBlue") +
    geom_sf(data = pl, aes(colour = pp), pch = 15, size = 2) +
    labs(color = "SSOS-TGOB (log10)", caption = paste0("occs source: AOPK NDOP (years: ", min(ndop.fs$years), "-", max(ndop.fs$years), "; months: ", min(ndop.fs$months), "-", max(ndop.fs$months), ") | author: Petr Balej | WGS84 | grid: KFME (SitMap_2Rad) | generated: ", Sys.Date())) +
    scale_colour_gradient(low = "white", high = "red", trans = "log10") +
    theme_light() +
    geom_sf(data = pres.unique) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    ggtitle(
      label = (druh),
      subtitle = paste0("SSOS-TGOB (same species observers occs per pixel) | p:", df.temp$pu, " / bg:", df.temp$ssosu, " / autors:", df.temp$authors)
    )

  png(paste0(path.ndop, "ssos/", str_replace_all(druh, "[\\/\\:]", "_"), "-natural.png"), width = 1200, height = 700)

  print(ggpl)
  dev.off()
}


saveRDS(df.out, paste0(path.ndop, "df.out.rds"))
write.csv(df.out, paste0(path.ndop, "df.out.csv"), row.names = FALSE)

saveRDS(df.out2, paste0(path.ndop, "df.out2.rds"))
write.csv(df.out2, paste0(path.ndop, "df.out2.csv"), row.names = FALSE)


# vygenerovat i mapy presencí druhu (počet) a backgroundu podle odpovídajících autorů (počty i per pixel)


# TODO: cross check - odpovídají názvy kvadrátů ze SiteName reálně poloze (Lat, Lon)? - udělat i prostorový průnik se sítí
# lsd.filter %>% left_join(sitmap_2rad.source, by = "POLE")
# st.sitmap_2rad.centroids <- st_centroid(st.sitmap_2rad)

# minimum x presenci a y absencí - možné až po získání prediktorů a "děr"

# + remove spatialy autocorrelated squares - geographically/environmentally
# alespoň reprezentativnost vůči altitude a landcoveru?

# do protokolu časy i adresu skriptu, který to vygeneroval + název vygenerovaného souboru

# https://github.com/PetrBalej/rgee/blob/master/R/clean/export_raster.R
# https://github.com/PetrBalej/rgee/blob/master/R/igaD.R

end_time <- Sys.time()
print(end_time - start_time)