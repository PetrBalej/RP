#
# Are the LSD PA sites selected randomly across Czechia?
#
# Two Sample Kolmogorov-Smirnov Test:
# test if two samples (LSD PA sites vs. ramdom sites) came from the same distribution.
# If p-value is less than 0.05, we reject the null hypothesis.
#

start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "raster", "ape") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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

lsd.pa.centroids.t <- st_coordinates(lsd.pa.centroids)

r.sampled.lsd <- extract(predictors, lsd.pa.centroids.t, df = TRUE)[, -1]

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
write.csv(results.mean, paste0(path.lsd, "ks.test_LSDxRandom.csv"))



##################################################################################################################################################

#
# Moran.I
#

# distance matrix

lsd.dists <- as.matrix(dist(cbind(lsd.pa.centroids.t[, 1], lsd.pa.centroids.t[, 2])))
lsd.dists.inv <- 1 / lsd.dists
diag(lsd.dists.inv) <- 0

results.mi <- data.frame(layer = character(), observed = numeric(), expected = numeric(), sd = numeric(), pvalue = numeric())

for (cn in names(r.sampled.lsd)) {
  print(cn)
  mi.res <- ape::Moran.I(r.sampled.lsd[[cn]], lsd.dists.inv, scaled = TRUE)
  results.mi %<>% add_row(layer = cn, observed = mi.res$observed, expected = mi.res$expected, sd = mi.res$sd, pvalue = mi.res$p.value)
}
print(results.mi)
# layer   observed     expected         sd       pvalue
# 1  l8_4_mndwi_mean 0.02744470 -0.006535948 0.01288603 8.363908e-03
# 2  l8_4_msavi_mean 0.02393740 -0.006535948 0.01289155 1.808745e-02
# 3 l8_4_raw_mean.B6 0.07930399 -0.006535948 0.01284880 2.376810e-11
# 4   l8_6_ndvi_mean 0.08085675 -0.006535948 0.01236914 1.601830e-12
# 5 l8_6_raw_mean.B1 0.10347677 -0.006535948 0.01270486 0.000000e+00
# 6  wc_0_mean.bio02 0.26097223 -0.006535948 0.01290445 0.000000e+00
# 7  wc_0_mean.bio04 0.21438172 -0.006535948 0.01288428 0.000000e+00
# 8  wc_0_mean.bio06 0.25266352 -0.006535948 0.01290913 0.000000e+00
# 9  wc_0_mean.bio19 0.23652256 -0.006535948 0.01277942 0.000000e+00
write.csv(results.mi, paste0(path.lsd, "moransI_4326.csv"))


# projekce např. do Křováka pro správný výpočet vzdáleností - vychází skoro stejně...

lsd.pa.centroids.t <- as.data.frame(lsd.pa.centroids %>% st_transform(5514) %>% st_coordinates())

# distance matrix
lsd.dists <- as.matrix(dist(cbind(lsd.pa.centroids.t[, 1], lsd.pa.centroids.t[, 2])))
lsd.dists.inv <- 1 / lsd.dists
diag(lsd.dists.inv) <- 0

results.mi <- data.frame(layer = character(), observed = numeric(), expected = numeric(), sd = numeric(), pvalue = numeric())

for (cn in names(r.sampled.lsd)) {
  print(cn)
  mi.res <- ape::Moran.I(r.sampled.lsd[[cn]], lsd.dists.inv, scaled = TRUE)
  results.mi %<>% add_row(layer = cn, observed = mi.res$observed, expected = mi.res$expected, sd = mi.res$sd, pvalue = mi.res$p.value)
}
print(results.mi)
# layer   observed     expected         sd       pvalue
# 1  l8_4_mndwi_mean 0.03093599 -0.006535948 0.01234945 2.410940e-03
# 2  l8_4_msavi_mean 0.02979441 -0.006535948 0.01235474 3.275769e-03
# 3 l8_4_raw_mean.B6 0.07313161 -0.006535948 0.01231377 9.813728e-11
# 4   l8_6_ndvi_mean 0.08583455 -0.006535948 0.01185415 6.661338e-15
# 5 l8_6_raw_mean.B1 0.10171464 -0.006535948 0.01217584 0.000000e+00
# 6  wc_0_mean.bio02 0.28386071 -0.006535948 0.01236710 0.000000e+00
# 7  wc_0_mean.bio04 0.21316936 -0.006535948 0.01234777 0.000000e+00
# 8  wc_0_mean.bio06 0.24331980 -0.006535948 0.01237159 0.000000e+00
# 9  wc_0_mean.bio19 0.24162623 -0.006535948 0.01224729 0.000000e+00
write.csv(results.mi, paste0(path.lsd, "moransI_5514.csv"))
