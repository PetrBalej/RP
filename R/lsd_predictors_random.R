#
# Are the LSD PA sites selected randomly across Czechia?
#
# Two Sample Kolmogorov-Smirnov Test:
# test if two samples (LSD PA sites vs. ramdom sites) came from the same distribution.
# If p-value is less than 0.05, we reject the null hypothesis.
#

start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "raster") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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
# path.rgee <- paste0(path.wd, "../../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
# source(paste0(path.rgee, "R/export_raster/functions.R"))
path.lsd <- paste0(path.prep, "../dataPrep/lsd/") # hodinovka / lsd


############
# inputs
############

predictors <- stack(readRDS(paste0(path.prep, "predictors/20192022/preds_final_vifs.rds"))) # vifstep2 pcamap6
predictors <- predictors / sum(as.vector(predictors), na.rm = TRUE)

# LSD PA
lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds")) %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci
lsd.pa.centroids %<>% group_by(POLE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(geometry)

############
# execution
############

n.points <- nrow(lsd.pa.centroids)
reps <- 100

results <- data.frame(rep = integer(), layer = character(), pvalue = numeric())

r.sampled.lsd <- extract(predictors, st_coordinates(lsd.pa.centroids), df = TRUE)[, -1]

for (rep in 1:reps) {
  print(rep)
  r.sampled.random <- as.data.frame(sampleRandom(predictors, n.points))
  for (cn in names(r.sampled.lsd)) {
    print(cn)
    ks.res <- ks.test(r.sampled.lsd[[cn]], r.sampled.random[[cn]])
    results %<>% add_row(rep = rep, layer = cn, pvalue = ks.res$p.value)
  }
}

results.mean <- results %>%
  group_by(layer) %>%
  mutate(pvalue.mean = mean(pvalue)) %>%
  mutate(rep = max(rep)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(-pvalue) 
  # %>% mutate(sameDistribution = ifelse(pvalue.mean > 0.05, 1, 0))

print(results.mean)
# # A tibble: 9 × 3
# rep layer            pvalue.mean
# <int> <chr>                  <dbl>
# 1   100 l8_4_mndwi_mean       0.600 
# 2   100 l8_4_msavi_mean       0.307 
# 3   100 l8_4_raw_mean.B6      0.272 
# 4   100 l8_6_ndvi_mean        0.328 
# 5   100 l8_6_raw_mean.B1      0.0600
# 6   100 wc_0_mean.bio02       0.153 
# 7   100 wc_0_mean.bio04       0.282 
# 8   100 wc_0_mean.bio06       0.0813
# 9   100 wc_0_mean.bio19       0.102 
