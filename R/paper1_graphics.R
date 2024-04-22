required_packages <- c("tidyverse", "sf", "magrittr", "stringi", "raster", "spatstat","spatstat.explore", "spatstat.geom","spatstat.random") 

install.packages(setdiff(required_packages, rownames(installed.packages())))

# načte všechny požadované knihovny jako dělá jednotlivě library()
lapply(required_packages, require, character.only = TRUE)


wd <- "C:/Users/petr/Documents/2024-04-11/"
setwd(wd)

kernels <- c(0.2, 0.3, 0.4, 0.5, 0.6) # c(0.25, 0.5, 1, 2, 3, 4)
px.sizes <- # c(0.25, 0.5, 1, 2, 4, 8)

  
  minMaxNormalize <- function(d, r, w = 1) {
    
    br.temp <- raster::resample(raster::raster(d), r, method = "bilinear")
    
    # normalize to 0-1
    r.min <- raster::minValue(br.temp)
    r.max <- raster::maxValue(br.temp)
    br.temp <- ((br.temp - r.min) / (r.max - r.min))
    br.temp <- raster::setMinMax(br.temp)
    return(br.temp/sum(as.vector(br.temp)) * w)
  
  }
  
  
  

# vzorové gridy
r.unit <- raster(nrows=2, ncols=3, xmn=0, xmx=3, ymn=0, ymx=2, vals=1:6)
plot(r.unit)
grid.unit <- rasterToPolygons(r.unit)
plot(grid.unit)

r.unit2 <- disaggregate(r.unit, fact = 2) 
grid.unit2 <- rasterToPolygons(r.unit2)
plot(grid.unit2)


r.high <- disaggregate(r.unit, fact = 100) 



u4 <- 0.25/2
points.a1 <- data.frame(
  x = c(0.25,  0.75, 0.75),
  y = c(1.25,  1.25,  1.75)
  )


points.a2 <- data.frame(
  x = c(1.75, 2.25, 2.25),
  y = c(1.25, 1.25, 0.75)
)

points.all <- rbind(points.a1, points.a2)

points.a1.sf <- sf::st_as_sf(points.a1, coords = c("x","y"))
points.a2.sf <- sf::st_as_sf(points.a2, coords = c("x","y"))
points.all.sf <- sf::st_as_sf(points.all, coords = c("x","y"))

plot(points.a1.sf, add=TRUE, col="red")
plot(points.a2.sf, add=TRUE, col="blue")


#
# ppp
#

# tgob
ext <- raster::extent(r.unit)
ow <- owin(xrange = c(ext@xmin, ext@xmax), yrange = c(ext@ymin, ext@ymax))
points.ppp <- ppp(points.all[, 1], points.all[, 2], window = ow)

bw <- bw.scott.iso(points.ppp) # sigma 0.3871331
#bw.diggle(points.ppp) # sigma 0.1247554
# bw.CvL(points.ppp) # sigma 0.6461961

# d <- density.ppp(points.ppp,  sigma = bw, positive = TRUE)



for(k in kernels){
  print(k)
  d <- density.ppp(points.ppp,  sigma = k, positive = TRUE)
  # plot(d, main=k)
  
  
  png(paste0(wd,"outputs/tgob-", k, ".png"), width=200, height=131)
  par(mar=c(2,2,2,0), cex=0.8, bg=NA)
  plot(minMaxNormalize(d, r.high), cex.lab=0.9, cex.axis=0.9, cex.main=0.9, cex.sub=0.9, legend=FALSE)
  plot(points.all.sf, add=TRUE)
  mtext(paste0("sigma=", k), side=3)
  dev.off()
  
  
}




# ooa

# o1, o2
points.ppp.a1 <- ppp(points.a1[, 1], points.a1[, 2], window = ow)
points.ppp.a2 <- ppp(points.a2[, 1], points.a2[, 2], window = ow)



# o1
d <- density.ppp(points.ppp.a1,  sigma = 0.3, positive = TRUE)
# plot(d, main=k)


png(paste0(wd,"outputs/ooa-a1-1-0.25.png"), width=200, height=131)
par(mar=c(2,2,2,0), cex=0.8, bg=NA)
r.a1.w1 <- minMaxNormalize(d, r.high, 0.25)
plot(r.a1.w1, cex.lab=0.9, cex.axis=0.9, cex.main=0.9, cex.sub=0.9, legend=FALSE)
plot(points.a1.sf, add=TRUE, col="red")
plot(grid.unit, add=TRUE)
mtext(paste0("sq=1; ","w=", 0.25), side=3)
dev.off()

png(paste0(wd,"outputs/ooa-a1-2-0.5.png"), width=200, height=131)
par(mar=c(2,2,2,0), cex=0.8, bg=NA)
r.a1.w2 <- minMaxNormalize(d, r.high, 0.5)
plot(r.a1.w2, cex.lab=0.9, cex.axis=0.9, cex.main=0.9, cex.sub=0.9, legend=FALSE)
plot(points.a1.sf, add=TRUE, col="red")
plot(grid.unit2, add=TRUE)
mtext(paste0("sq=3; ","w=", 0.5), side=3)
dev.off()


# o2
d <- density.ppp(points.ppp.a2,  sigma = 0.3, positive = TRUE)
# plot(d, main=k)


png(paste0(wd,"outputs/ooa-a2-1-0.75.png"), width=200, height=131)
par(mar=c(2,2,2,0), cex=0.8, bg=NA)
r.a2.w1 <- minMaxNormalize(d, r.high, 0.75)
plot(r.a2.w1, cex.lab=0.9, cex.axis=0.9, cex.main=0.9, cex.sub=0.9, legend=FALSE)
plot(points.a2.sf, add=TRUE, col="blue")
plot(grid.unit, add=TRUE)
mtext(paste0("sq=3; ","w=", 0.75), side=3)
dev.off()

png(paste0(wd,"outputs/ooa-a2-2-0.5.png"), width=200, height=131)
par(mar=c(2,2,2,0), cex=0.8, bg=NA)
r.a2.w2 <- minMaxNormalize(d, r.high, 0.5)
plot(r.a2.w2, cex.lab=0.9, cex.axis=0.9, cex.main=0.9, cex.sub=0.9, legend=FALSE)
plot(points.a2.sf, add=TRUE, col="blue")
plot(grid.unit2, add=TRUE)
mtext(paste0("sq=3; ","w=", 0.5), side=3)
dev.off()




# o1+o2


png(paste0(wd,"outputs/ooa-a-1.png"), width=200, height=131)
par(mar=c(2,2,2,0), cex=0.8, bg=NA)
plot(r.a1.w1+r.a2.w1, cex.lab=0.9, cex.axis=0.9, cex.main=0.9, cex.sub=0.9, legend=FALSE)
plot(points.a1.sf, add=TRUE, col="red")
plot(points.a2.sf, add=TRUE, col="blue")
plot(grid.unit, add=TRUE)
mtext(paste0("sq=4; ","w=", 1.0), side=3)
dev.off()

png(paste0(wd,"outputs/ooa-a-2.png"), width=200, height=131)
par(mar=c(2,2,2,0), cex=0.8, bg=NA)
plot(r.a1.w2+r.a2.w2, cex.lab=0.9, cex.axis=0.9, cex.main=0.9, cex.sub=0.9, legend=FALSE)
plot(points.a1.sf, add=TRUE, col="red")
plot(points.a2.sf, add=TRUE, col="blue")
plot(grid.unit2, add=TRUE)
mtext(paste0("sq=6; ","w=", 1.0), side=3)
dev.off()




