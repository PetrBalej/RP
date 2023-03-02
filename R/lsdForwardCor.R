cmd_arg <- commandArgs(trailingOnly = TRUE)
if (is.na(cmd_arg[1])) {
  print("*********************************************************************************************************************")
  print("nezadán parametr pro část druhů")
  print("*********************************************************************************************************************")
  cmd_arg <- 1
} else {
  cmd_arg <- cmd_arg[1]
}

print(cmd_arg)

start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "raster") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############

# nově: pc02/rgee20230227/; dataPrep/lsdForwardCor

# nastavit working directory
path.wd <- "C:/Users/petr/Documents/pc02/projects/"
setwd(path.wd)
path.data <- "C:/Users/petr/Documents/pc02/projects-data/"
path.rgee <- "C:/Users/petr/Documents/pc02/rgee20230227/" # samsung500ntfs # paste0(path.expand("~"), "/Downloads/rgee2/rgee")
path.wd.prep <- paste0(path.wd, "dataPrep/lsdForwardCor/")
# source(paste0(path.rgee, "R/export_raster/functions.R"))

############
# inputs (reálně už zprocesované různě dříve - nahradit orig. jasně zpracovanými daty později)
############
df <- readRDS(paste0(path.rgee, "R/kostelec2023/NB/clean/lsd-PA-FINAL-57preds-118sp.df.rds"))
remove <- readRDS(paste0(path.rgee, "R/kostelec2023/NB/k2023/preds07-remove.rds"))[-1] # jedna záměna ve VIFu
rasterStack100na23 <- readRDS(paste0(path.rgee, "R/kostelec2023/NB/clean/raster_stack_100_na_23.rds")) # problematické DM - CV - inf NoData

# všechny rastery normalizovat/standardizovat (nutné pro RF?)


############
# settings
############
lsd.fs <- list("minSpecies" = 30, "version" = "v1", "speciesPerGroup" = 40, "speciesPart" = cmd_arg)


############
# execution
############

# preventivní znáhodnění, ale asi netřeba, dělá sdm
# rozdělení druhů do skupin
set.seed(85)
df <- df[sample(1:nrow(df)), ]
rownames(df) <- NULL
df <- df[, !(names(df) %in% remove)]
speciesNames <- names(df)[24:141]
speciesParts <- split(speciesNames, ceiling(seq_along(speciesNames) / lsd.fs$speciesPerGroup))


# redukce na 16 prediktorů
predsNames <- names(df)[1:23] # 1:15 L8;  16:23 Bio
remove2 <- c("l8_4_raw_cv_B5", "l8_6_raw_cv_B7", "l8_6_raw_cv_B1", "wc_mean_bio04", "l8_6_ndvi_cv", "l8_5_ndvi_cv", "wc_cv_bio04")
positive2 <- c("l8_4_evi_mean", "l8_4_mndwi_cv", "l8_4_ndvi_cv", "l8_4_raw_cv_B1", "l8_5_mndwi_mean", "l8_5_ndvi_mean", "l8_5_raw_cv_B3", "l8_5_raw_mean_B5", "l8_6_raw_cv_B4", "l8_6_raw_mean_B1", "wc_cv_bio06", "wc_cv_bio11", "wc_cv_bio15", "wc_mean_bio02", "wc_mean_bio06", "wc_mean_bio15")
predsNames <- predsNames[predsNames %in% positive2]
rasterStack16 <- rasterStack100na23[[predsNames]]



# freq(rasterStack16[[1]])
# nasamplovat 10 BG variant - nebo to nechat udělat sdm - jako u LG (10% bude stačit?)

sample.random <- list()
for (i in 1:10) {
  sample.random[[i]] <- sampleRandom(rasterStack16, size = 1000, cells = FALSE, sp = FALSE)
}

# plot(rasterStack16[[1]])
# par(bg=NA)
# plot(sample.random, add=TRUE)



# nemůžu, prediktory stejně vybírám druhově specificky...
# speciesParts.f <- paste(speciesParts[[lsd.fs["speciesPart"]]], collapse = "+")



d <- sdm::sdmData(as.formula(paste0(sp, "~", pred.names.all.f, "+coords(x+y)")), train = df)

m <- sdm::sdm(as.formula(paste0(sp, "~", pred.names.all.f)), data = d, methods = c("glm"), replication = c("cv"), cv.folds = 3, n = 5, seed = TRUE)


# https://rdrr.io/cran/sdm/man/predict.html
# w 	numeric, specifies which model(s) should be used if the object contains several models; with NULL all models are used



sdm::predict(m, sample.random,
  # filename=paste0(export_path, "export/",ct,"/",sps,".tif"),
  mean = TRUE, overwrite = TRUE
)






# do protokolu časy i adresu skriptu, který to vygeneroval + název vygenerovaného souboru

# https://github.com/PetrBalej/rgee/blob/master/R/clean/export_raster.R
# https://github.com/PetrBalej/rgee/blob/master/R/igaD.R

end_time <- Sys.time()
print(end_time - start_time)