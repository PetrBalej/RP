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
required_packages <- c("tidyverse", "sf", "magrittr", "raster", "sdm", "usdmPB") # c("sp", "rgdal", "mapview", "raster", "geojsonio", "stars", "httpuv", "tidyverse", "sf", "lubridate", "magrittr", "dplyr", "readxl", "abind", "stringr")

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)

############
# paths
############
"C:/Users/petr/ownloads/RP/RP"
# nově: pc02/rgee20230227/; dataPrep/lsdForwardCor

# nastavit working directory
path.wd <- "C:/Users/petr/Downloads/RP/RP/"
setwd(path.wd)
path.data <- "C:/Users/petr/Downloads/RP/projects-data/"
path.rgee <- "C:/Users/petr/Documents/pc02/rgee20230227/" # samsung500ntfs # paste0(path.expand("~"), "/Downloads/rgee2/rgee")
path.wd.prep <- paste0(path.wd, "dataPrep/lsdForwardCor/")
source(paste0(path.rgee, "R/export_raster/functions.R"))

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
lsd.fs <- list("minSpecies" = 30, "version" = "v1", "speciesPerGroup" = 5, "speciesPart" = cmd_arg)


############
# execution
############

"%notin%" <- Negate("%in%")
# preventivní znáhodnění, ale asi netřeba, dělá sdm
# rozdělení druhů do skupin
set.seed(85)
# df <- df[sample(1:nrow(df)), ]
# rownames(df) <- NULL
df <- df[, !(names(df) %in% remove)]
speciesNames <- names(df)[24:141]
speciesParts <- split(speciesNames, ceiling(seq_along(speciesNames) / lsd.fs$speciesPerGroup))


# redukce na 16 prediktorů
predsNames <- names(df)[1:23] # 1:15 L8;  16:23 Bio
remove2 <- c("l8_4_raw_cv_B5", "l8_6_raw_cv_B7", "l8_6_raw_cv_B1", "wc_mean_bio04", "l8_6_ndvi_cv", "l8_5_ndvi_cv", "wc_cv_bio04")
positive2 <- c("l8_4_evi_mean", "l8_4_mndwi_cv", "l8_4_ndvi_cv", "l8_4_raw_cv_B1", "l8_5_mndwi_mean", "l8_5_ndvi_mean", "l8_5_raw_cv_B3", "l8_5_raw_mean_B5", "l8_6_raw_cv_B4", "l8_6_raw_mean_B1", "wc_cv_bio06", "wc_cv_bio11", "wc_cv_bio15", "wc_mean_bio02", "wc_mean_bio06", "wc_mean_bio15")
predsNames <- predsNames[predsNames %in% positive2]
# rasterStack16 <- rasterStack100na23[[predsNames]]
predsNames <- predsNames[1:5]




varCoCo <- function(df, sp, preds = NULL, preds.orig, predsK = 2, topN = 3, dAucMax = 0.05, speciesPart = NULL, first = TRUE) {
  print("**********************************************************************")
  if (is.null(preds)) {
    predsNamesComb <- comb_k(preds.orig, predsK)
    preds.n <- length(preds.orig)
  } else {
    # dělám couples - doplňuji do all
    predsNamesComb <- list()

    for (pcs in preds) {
      print("pcs------")
      print(pcs)
      pcs.vector <- sort(unlist(strsplit(pcs, "\\+")))

      print(pcs.vector)
      preds.orig.notin <- preds.orig[preds.orig %notin% pcs.vector]
      print(preds.orig.notin)

      if (is_empty(preds.orig.notin)) {
        return(NULL)
      }
      print("ccc")
      # print(comb_k(c(pcs, preds.orig.notin ), predsK))
      # print(length(comb_k(c(pcs, preds.orig.notin ), predsK)))
      # rovnou dávat připravené kombinace prediktorů do vstupu????,
      print("ccca")
      # predsNamesComb <- append(predsNamesComb, comb_k(c(pcs, preds.orig.notin ), predsK))


      predsNamesComb.temp <- list()
      for (po in preds.orig.notin) {
        predsNamesComb.temp[[po]] <- c(pcs, po)
      }


      # počítá jen s dvojkombinacema!!!
      predsNamesComb <- append(predsNamesComb, predsNamesComb.temp)
      print(predsNamesComb)
      print(length(predsNamesComb))
    }

    print("ccc-total")
    print(predsNamesComb)
    print(length(predsNamesComb))


    print("ccc-total")


    stop()
    preds.n <- length(preds.orig.notin) + 1
  }

  #  predsNamesComb <- append(predsNamesComb, lapply(comb_k(preds.orig.notin, predsK), function(x, y) x <- c(x, y), y = pcs))


  print(preds)
  print(preds.n)


  print("počet všech kombinací:")
  print(length(predsNamesComb))
  print(predsNamesComb)

  # predsCouple <-  preds[str_detect(preds, "\\+")] # odstranit z orig. seznamu prediktorů
  # predsCouple.vector <- sort(unlist(strsplit(predsCouple, "\\+")))
  # preds.real.vector <- sort(unlist(strsplit(preds, "\\+")))
  # preds.real.vector.n <- length(preds.real.vector)



  # if (preds.n > 1) {
  # predsNamesComb <- comb_k(preds, predsK) # dělám kombinace už se spojenou couple zdánlivě jako jedním prediktorem (pevnou kombinací)
  # }else{
  #   predsNamesComb <- preds
  # }

  if (preds.n > 0) {
    predsCouples <- list()
    for (predsNameComb in predsNamesComb) {
      print(predsNameComb)

      predsNameComb <- sort(unique(unlist(strsplit(predsNameComb, "\\+"))))


      pred.names.all.f <- paste(predsNameComb, collapse = "+")

      d <- sdm::sdmData(as.formula(paste0(sp, "~", pred.names.all.f, "+coords(x+y)")), train = df)

      m <- sdm::sdm(as.formula(paste0(sp, "~", pred.names.all.f)), data = d, methods = c("glm"), replication = c("cv"), cv.folds = 3, n = 3, seed = TRUE)

      ### train
      df.merge.train.temp <- getEvaluation(m,
        wtest = "training",
        stat = c("AUC", "COR", "Deviance", "obs.prevalence", "threshold", "sensitivity", "specificity", "TSS", "Kappa", "NMI", "phi", "ppv", "npv", "ccr", "prevalence")
      )
      ### test
      df.merge.test.temp <- getEvaluation(m,
        wtest = "test.dep",
        stat = c("AUC", "COR", "Deviance", "obs.prevalence", "threshold", "sensitivity", "specificity", "TSS", "Kappa", "NMI", "phi", "ppv", "npv", "ccr", "prevalence")
      )





      df.merge.train.temp$preds <- df.merge.test.temp$preds <- pred.names.all.f
      # df.merge.train$predsCouple <- df.merge.test$predsCouple <- predsCouple

      auc.test <- median(df.merge.test.temp$AUC)
      dAuc <- auc.test - median(df.merge.train.temp$AUC)
      df.merge.train.temp$dAuc <- df.merge.test.temp$dAuc <- dAuc


      df.merge.train.temp$predsCouples <- df.merge.test.temp$predsCouples <- 0
      # k dalšímu postupu jen dobré modely - ale dělám to celkově (ne individuální porovnání - jen pro průchod - postprocesem ještě asi později dofiltruju další...)
      if (abs(dAuc) < dAucMax) {
        predsCouples[[pred.names.all.f]] <- auc.test
        df.merge.train.temp$predsCouples <- df.merge.test.temp$predsCouples <- 1
      }



      # merge train a test s info
      df.merge.test <- merge(df.merge.test.temp,
        getModelInfo(m),
        by = "modelID"
      )
      df.merge.train <- merge(df.merge.train.temp,
        getModelInfo(m),
        by = "modelID"
      )
      print("path......................................................")
      print(paste0(path.wd.prep, speciesPart, "-train.csv"))
      write.table(df.merge.train, file = paste0(path.wd.prep, speciesPart, "-train.csv"), append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = first)
      write.table(df.merge.test, file = paste0(path.wd.prep, speciesPart, "-test.csv"), append = TRUE, quote = TRUE, sep = ",", row.names = FALSE, col.names = first)

      if (first) {
        first <- FALSE
      }

      gc()
    }

    if (length(predsCouples) > 0 & preds.n > 1 & length(predsNameComb) != length(preds.orig)) {
      print("rekurze!!!!!!!!!!!!!!!")


      # + výběr rozpětí test-train
      # vyhodnotit (nutný zásobník AUC+predsCouple) a vybrat nej X kombinací do dalšího volání
      # ve výchozím stavu ale žádná couple není - netřeba - couple je předchozí iterace,
      # já musím získat jen 3 nej kombinace pro aktuální K a z něj vytvořit couple do další iterace


      predsCoupleTop <- names(na.omit(sort(unlist(predsCouples), decreasing = TRUE)[1:topN]))

      varCoCo(df, sp, preds = predsCoupleTop, preds.orig, speciesPart = speciesPart, first = FALSE)
    }
  }
}







first <- TRUE
for (sp in speciesParts[[lsd.fs$speciesPart]]) {
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print(sp)
  varCoCo(df, sp, preds = NULL, predsNames, speciesPart = lsd.fs$speciesPart, first = first)
  if (first) {
    first <- FALSE
  }
}











# k6_top <- k6 %>% group_by(species) %>% slice_max(AUC, n = 3)
# k6_top  %>% group_by(species) %>% summarize(allPreds = length(unique(strsplit(paste(preds, collapse = "+"), "\\+")[[1]])))







# do protokolu časy i adresu skriptu, který to vygeneroval + název vygenerovaného souboru

# https://github.com/PetrBalej/rgee/blob/master/R/clean/export_raster.R
# https://github.com/PetrBalej/rgee/blob/master/R/igaD.R

end_time <- Sys.time()
print(end_time - start_time)