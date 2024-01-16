start_time <- Sys.time()

# kontrola (do)instalace všech dodatečně potřebných balíčků
required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat", "ENMeval", "sdmATM", "ggplot2", "grid", "gridExtra", "ggtext", "rmarkdown", "grDevices") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

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

# "C:\Program Files\R\R-4.2.1\bin\x64\Rscript.exe" "C:\Users\petr\Documents\2023-03-20\RP\RP\R\ndopTGOBeval.R"

# nastavit working directory
# path.wd <- "/mnt/2AA56BAE3BB1EC2E/Downloads/rgee2/RP/RP/" # "D:/PERSONAL_DATA/pb/RP20230313/RP/RP/"
setwd(path.wd)
path.data <- paste0(path.wd, "../../projects-data/") # "D:/PERSONAL_DATA/pb/RP20230313/RP/projects-data/"
path.prep <- paste0(path.wd, "../../dataPrep/")
source(paste0(path.wd, "shared.R"))
path.rgee <- paste0(path.wd, "../../../rgee/") # "D:/PERSONAL_DATA/pb/kostelec2023/RP-fw/rgee20230303/"
source(paste0(path.rgee, "R/export_raster/functions.R"))
path.lsd <- paste0(path.prep, "lsd/")
path.eval <- paste0(path.prep, "ndopTGOBevalimp999mapy/")
path.tgob <- paste0(path.prep, "ndopTGOBimp999m/")
path.tgob.tgbg <- paste0(path.prep, "ndopTGOBtgbg999m/") # path.wd.prep.ndop
path.PP <- paste0(path.eval, "PP/")

############
# inputs
############
# uložit si results, rastery a udělat evaluate; ukládat i počty presencí a bg (z modelů) a WKT presencí (NDOP i LSD)
lsd.pa.centroids <- readRDS(paste0(path.lsd, "lsd.pa.centroids.rds")) %>% filter(POLE != "5260ca") # krkonoše, nejsou v prediktorech (automatizovat odebrání!!) - až v evaluaci

predictors <- stack(readRDS(paste0(path.prep, "predictors/20192022/preds_final_vifs.rds"))) # vifstep2 pcamap6
predictors <- predictors / sum(as.vector(predictors), na.rm = TRUE)
crs <- raster::crs(predictors[[1]])


############
# settings
############

ndop.fs <- list("aucTresholds" = c(0.00, 0.70), "version" = "v1")

############
# execution
############

se <- function(x) sqrt(var(x) / length(x))

sumNormal <- function(x) {
  return(x / sum(as.vector(x), na.rm = TRUE))
}
# tofo: doplnit možnost "un" jako zobrazované samostatné varianty!!!. jen nebude figurovat v diff grafech, tam je to nesmysl

perform_ttest <- TRUE
if (perform_ttest) {
  breakNum <- 3
} else {
  breakNum <- 333
}

# select null if performs better
withNull <- FALSE
# remove thin to get only TGOB derived versions
withThin <- TRUE

# tgobXostatní
anyAlt <- FALSE


"%notin%" <- Negate("%in%")
dir.create(path.PP, showWarnings = FALSE)
dir.create(paste0(path.PP, "predictions01"), showWarnings = FALSE)
dir.create(paste0(path.PP, "predictions"), showWarnings = FALSE)
dir.create(paste0(path.PP, "predictions_individual"), showWarnings = FALSE)


# z ndopTGOBevalPP_wilcox_groupKernelTgob.R

modelsResults.avg <- paste0(path.eval, "tbl.final.rds")
if (file.exists(modelsResults.avg)) {
  print("načítám existující tbl.final")
  tbl <- readRDS(modelsResults.avg)
} else {
  stop()
}



tbl.f <- tbl %>%
  group_by(version.old, species) %>%
  slice_max(AUC, with_ties = FALSE) %>%
  ungroup() %>%
  filter(version.old != "un")

versions <- unique(tbl.f$version.old)
print(versions)
sps <- unique(tbl.f$species)
reps <- 1:5


# pdf(filename=paste0(path.PP,"predictions.pdf"), width = 4, height = 3)

additionalInfo <- list()

for (sp in sps) {
  print(sp)

  tbl.f.sp <- tbl.f %>% filter(species == sp)

  predictions.mean <- list()
  predictions.mean.names <- c()
  for (ver in versions) {
    print(ver)

    tbl.f.sp.ver <- tbl.f.sp %>% filter(version.old == ver)

    ver.orig <- tbl.f.sp.ver$version.group
    if (tbl.f.sp.ver$version.group == "random" | tbl.f.sp.ver$version.group == "rthin") {
      ver.orig <- "un"
    }


    predictions <- list()
    for (rep in reps) {
      print(rep)
      predictions[[rep]] <- readRDS(paste0(path.tgob, "34ssos.new_", sp, "_", rep, ".rds"))[[sp]][[ver.orig]][[tbl.f.sp.ver$adjust]][[tbl.f.sp.ver$fc]][[tbl.f.sp.ver$rm]]

      # `34ssos.new_Accipiter gentilis_5`[["Accipiter gentilis"]][["kernel0112000"]][["1"]][["LQP"]][["1"]]
    }

    predictions.mean.names <- c(predictions.mean.names, paste0(ver, " (AUC=", round(tbl.f.sp.ver$AUC, 3), ")"))
    predictions.mean[[tbl.f.sp.ver$version.old]] <- raster::calc(stack(predictions), median)


    additionalInfo[[sp]][[paste0("AUC_", tbl.f.sp.ver$version.old)]] <- tbl.f.sp.ver$AUC


    # obrázky predikcí samostatně
    png(filename = paste0(path.PP, "predictions_individual/", sp, "_", ver, ".png"), width = 1000, height = 800)
    plot(predictions.mean[[tbl.f.sp.ver$version.old]], main = paste0(sp, " / ", ver, " / (AUC=", round(tbl.f.sp.ver$AUC, 3), ")"), sub = paste0("adj=", tbl.f.sp.ver$adjust, ", FC=", tbl.f.sp.ver$fc, ", rm=", tbl.f.sp.ver$rm))
    dev.off()


    # bias rastery
    if (tbl.f.sp.ver$version.old == "OOA") {
      nw <- readRDS(paste0(path.tgob.tgbg, str_replace(ver.orig, "kernel", "kdr"), "_", sp, ".rds"))[[tbl.f.sp.ver$adjust]]


      br.temp <- raster::resample(raster::raster(sumNormal(nw)), predictors[[1]], method = "bilinear") # sumNormal

      # normalize to 0-1
      r.min <- raster::minValue(br.temp)
      r.max <- raster::maxValue(br.temp)
      br.temp <- ((br.temp - r.min) / (r.max - r.min))
      br.temp <- raster::setMinMax(br.temp)

      crs <- raster::crs(predictors[[1]])
      raster::crs(br.temp) <- crs
      br.temp <- raster::mask(raster::crop(br.temp, raster::extent(predictors[[1]])), predictors[[1]])
      br.temp <- raster::setMinMax(br.temp)

      predictions.mean[["OOA.bias"]] <- br.temp
      predictions.mean.names <- c(predictions.mean.names, "OOA.bias")
      # `kdr0112000_Accipiter gentilis`[["0.1"]]
    }


    if (tbl.f.sp.ver$version.old == "TGOB") {
      # nw <- readRDS(paste0(path.tgob, "50ssos_collector.subset_", sp, "_1.rds"))[[ver.orig]][["kd"]][[tbl.f.sp.ver$adjust]]
      # readRDS("D:/PersonalWork/Balej/v2/RP/dataPrep/ndopTGOBimp999m/50ssos_collector.subset_Vanellus vanellus_4.rds")
      nw <- readRDS(paste0(path.tgob.tgbg, "kdr90000999", "_", sp, ".rds"))[[tbl.f.sp.ver$adjust]]

      br.temp <- raster::resample(raster::raster(sumNormal(nw)), predictors[[1]], method = "bilinear") # sumNormal

      # normalize to 0-1
      r.min <- raster::minValue(br.temp)
      r.max <- raster::maxValue(br.temp)
      br.temp <- ((br.temp - r.min) / (r.max - r.min))
      br.temp <- raster::setMinMax(br.temp)

      crs <- raster::crs(predictors[[1]])
      raster::crs(br.temp) <- crs
      br.temp <- raster::mask(raster::crop(br.temp, raster::extent(predictors[[1]])), predictors[[1]])
      br.temp <- raster::setMinMax(br.temp)

      predictions.mean[["TGOB.bias"]] <- br.temp
      predictions.mean.names <- c(predictions.mean.names, "TGOB.bias")
      # `50ssos_collector.subset_Vanellus vanellus_4`[["TGOB01"]][["kd"]][["0.5"]]
    }
  }
  no <- c(5, 4, 2, 6, 3, 1)




  additionalInfo[[sp]]$cor_br <- cor(na.omit(as.vector(predictions.mean[["OOA.bias"]])), na.omit(as.vector(predictions.mean[["TGOB.bias"]])))
  additionalInfo[[sp]]$cor_tgobOoa <- cor(na.omit(as.vector(predictions.mean[["OOA"]])), na.omit(as.vector(predictions.mean[["TGOB"]])))

  predictions.mean <- subset(stack(predictions.mean), no)

  # predikce normalizované 01
  png(filename = paste0(path.PP, "predictions01/", sp, ".png"), width = 1400, height = 600)

  # par(mfrow = c(3, 3))
  # par(mar = c(0,0,5,0))
  # par(mar = c(2, 2, 20, 2))
  # plot(predictions.mean,zlim=c(0,1),main=predictions.mean.names[no])
  plot(predictions.mean, zlim = c(0, 1), main = predictions.mean.names[no], cex.main = 1.5)
  mtext(sp, side = 3, line = -1, outer = TRUE, adj = 0, cex = 1)
  mtext(paste0("cor: tgobOoa=", round(additionalInfo[[sp]]$cor_tgobOoa, 3), " / br=", round(additionalInfo[[sp]]$cor_br, 3)), side = 3, line = -1, outer = TRUE, adj = 1, cex = 1)
  dev.off()



  # predikce NEnormalizované
  png(filename = paste0(path.PP, "predictions/", sp, ".png"), width = 1400, height = 600)
  plot(predictions.mean, main = predictions.mean.names[no], cex.main = 1.5)
  mtext(sp, side = 3, line = -1, outer = TRUE, adj = 0, cex = 1)
  mtext(paste0("cor: tgobOoa=", round(additionalInfo[[sp]]$cor_tgobOoa, 3), " / br=", round(additionalInfo[[sp]]$cor_br, 3)), side = 3, line = -1, outer = TRUE, adj = 1, cex = 1)
  dev.off()
}

saveRDS(additionalInfo, paste0(path.PP, "additionalInfo.rds"))
# dev.off()