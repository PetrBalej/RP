# testovací příklad podle https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html

# Load packages -- the order here is important because some pkg functions overwrite others.
library(ENMeval)
library(raster)
library(dplyr)

# Set a random seed in order to be able to reproduce this analysis.
set.seed(48)

# You can search online databases like GBIF using the spocc package (commented below),
# but here we will load in some pre-downloaded data.
bv <- spocc::occ("Bradypus variegatus", "gbif", limit = 300, has_coords = TRUE)
occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[, 2:3])
# occs <- readRDS("bvariegatus.rds")

# Removing occurrences that have the same coordinates is good practice to
# avoid pseudoreplication.
occs <- occs[!duplicated(occs), ]

# Locate the predictor raster files from the dismo folder.
envs.files <- list.files(
  path = paste(system.file(package = "dismo"), "/ex", sep = ""),
  pattern = "grd", full.names = TRUE
)


# Read the raster files into a RasterStack.
# These variables represent 8 bioclimatic variables and one categorical variable "biome".
# Find the descriptions of the bioclimatic variables here:
# https://www.worldclim.org/data/bioclim.html
envs <- raster::stack(envs.files)
# The biome raster has some NAs for cells that have values in the other rasters.
# Let's mask all rasters to biome to change the value of these cells to NA for all rasters.
# ENMeval will do this automatically, but let's do it here to avoid the warning message later.
# We change back from a RasterBrick to RasterStack because of issues with assigning
# factor rasters for RasterBricks.
envs <- raster::mask(envs, envs[[9]]) %>% raster::stack()
# Make sure to declare the categorical variable as a factor
envs$biome <- raster::as.factor(envs$biome)
# Let's now remove occurrences that are cell duplicates -- these are
# occurrences that share a grid cell in the predictor variable rasters.
# Although Maxent does this by default, keep in mind that for other algorithms you may
# or may not want to do this based on the aims of your study.
# Another way to space occurrence records a defined distance from each other to avoid
# spatial autocorrelation is with spatial thinning (Aiello-Lammens et al. 2015).
occs.cells <- raster::extract(envs[[1]], occs, cellnumbers = TRUE)
occs.cellDups <- duplicated(occs.cells[, 1])
occs <- occs[!occs.cellDups, ]



# There are some points east of the Amazon River.
# Suppose we know this is a population that we don't want to include in the model.
# We can remove these points from the analysis by subsetting the occurrences by
# latitude and longitude.
occs <- filter(occs, latitude > -20, longitude < -45)

# First we extract the climatic variable values at the occurrence points -- these values are
# our "reference".
# We remove the categorical variable for these operations because the math only makes sense
# for continuous variables -- the function will not work with categorical variables.
occs.z <- raster::extract(envs[[-9]], occs)
# Now we use the similarity() function (borrowed from the rmaxent package) to calculate
# environmental similarity metrics of our predictor variable extent compared to the reference
# points.
library(rmaxent)
occs.sim <- similarity(envs[[-9]], occs.z)
occs.mess <- occs.sim$similarity_min
# This is the MESS plot -- increasingly negative values represent increasingly different
# climatic conditions from the reference (our occurrences), while increasingly positive
# values are more similar. First, we'll make a SpatialPoints object for our occurrences
# for plotting with levelplot() from the rasterVis package (Lamigueiro & Hijmans 2021).
# This package has great plotting functionality for rasters, and by default bins values for
# display when data is continuous.
occs.sp <- sp::SpatialPoints(occs)


# We'll now experiment with a different spatial R package called sf (simple features).
# Let's make our occs into a sf object -- as the coordinate reference system (crs) for these
# points is WGS84, a geographic crs (lat/lon) and the same as our envs rasters, we specify it
# as the RasterStack's crs.
occs.sf <- sf::st_as_sf(occs, coords = c("longitude", "latitude"), crs = raster::crs(envs))

# Now, we project our point data to an equal-area projection, which converts our
# degrees to meters, which is ideal for buffering (the next step).
# We use the typical Eckert IV projection.
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)

# Buffer all occurrences by 500 km, union the polygons together
# (for visualization), and convert back to a form that the raster package
# can use. Finally, we reproject the buffers back to WGS84 (lat/lon).
# We choose 500 km here to avoid sampling the Caribbean islands.
occs.buf <- sf::st_buffer(occs.sf, dist = 500000) %>%
  sf::st_union() %>%
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))


# First, let's specify a fake testing occurrence dataset and plot the testing points with
# the rest of our data
occs.testing <- data.frame(longitude = -runif(10, 55, 65), latitude = runif(10, -10, 0))





# Crop environmental rasters to match the study extent
envs.bg <- raster::crop(envs, occs.buf)
# Next, mask the rasters to the shape of the buffers
envs.bg <- raster::mask(envs.bg, occs.buf)



### vlastní random body
bg <- dismo::randomPoints(envs.bg[[1]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(occs)

# rand <- get.randomkfold(occs, bg, k = 5)
# evalplot.grps(pts = occs, pts.grp = rand$occs.grp, envs = envs.bg)



e.mx <- ENMevaluate(
  occs = occs, envs = envs.bg[[1:8]],
  bg.coords = bg,
  algorithm = "maxnet", partitions = "randomkfold",
  # kfolds = 5,
  tune.args = list(fc = c("L", "LQ"), rm = 1:5)
)

# data pro evaluaci + null model
res <- eval.results(e.mx)
mod.null <- ENMnulls(e.mx, mod.settings = list(fc = "LQ", rm = 5), no.iter = 100)
