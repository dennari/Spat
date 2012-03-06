offsetX <- 3478085
offsetY <- 7085263

cut.val <- 0 ### Just to force it.
theme.novpadding <-
    list(layout.heights =
         list(top.padding = cut.val,
              main.key.padding = cut.val,
              key.axis.padding = cut.val,
              axis.xlab.padding = cut.val,
              xlab.key.padding = cut.val,
              key.sub.padding = cut.val,
              bottom.padding = cut.val),
         layout.widths =
         list(left.padding = cut.val,
              key.ylab.padding = cut.val,
              ylab.axis.padding = cut.val,
              axis.key.padding = cut.val,
              right.padding = cut.val))
exportFigs <- 1
displayFigs <- 0

require("maptools")
require("spdep")
require("RColorBrewer")
require("spatstat")

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  jet.colorf <-
  colorRamp(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

mar_lab <- c(2.5,2.5,1.5,1.0)
mar_tight <- c(0.1,0.1,0.1,0.1)
myplot <- function(...,plotfn=plot,width=6,height=6,mar=mar_lab,file=FALSE,nodevoff=FALSE,afterfn=NULL,k=NULL) {
	if(displayFigs) {
		quartz()
		par(mar=mar)
		p <- plotfn(...)
		if(is.function(afterfn)) {
			afterfn(p,k)
		}
		if(!exportFigs) {
			return(p)
		}
	}
	if(exportFigs && file != FALSE) {
		pdf(file=file,width=width,height=height)
		par(mar=mar)
		p <- plotfn(...)
		if(is.function(afterfn)) {
			afterfn(p,k)
		}
		if(!nodevoff) {
			dev.off()
		}
		return(p)
	}
}
listplot <- function(k,v,file=FALSE,formula=FALSE,main="",...) {
	if(formula != FALSE) {
		p <- myplot(v,formula,file=sprintf(file,k),main=sprintf(main,k),k=k,...)
	} else {
		p <- myplot(v,file=sprintf(file,k),main=sprintf(main,k),k=k,...)
	}
	return(p)
}

#raj <- readShapeSpatial("map/HallintoalueRaja1Milj_ETRS.shp")
if(0) {
	map <- readShapeSpatial("map/HallintoAlue1Milj_ETRS.shp")
	hcont <- readShapeSpatial("map/hcont_l.shp")
	forest <- readShapeSpatial("map/forest.shp")
	
	population <- read.csv(
		"asukastiheys.csv",
		sep=",",
		header=FALSE,
		skip=0,
		col.names=c("Kunta","land_area","water_area","water_area2","total_area","population","density_area","density_land_area"))
	population <- population[,c("Kunta","land_area","density_land_area")]

	obs2009 <- read.csv(	
		"Obs2009XYS_ETRS.csv",
		sep=",",
		header=FALSE,
		skip=0,
		col.names=c("x","y","laji"))
	
	valtraj <- hcont[hcont@data$Z==0|hcont@data$Z==9999,]


	map@data <- map@data[,c("Kunta","Maakunta")]
	map@data$Kunta <- as.numeric(as.vector(map@data$Kunta))
	map@data$Maakunta <- as.numeric(as.vector(map@data$Maakunta))
	
	centroids <- coordinates(map)
}


comp <- function(k) {
	if(sum(population$Kunta==k) > 0) {
		tmp <- population[population$Kunta==k,c("density_land_area")]
		return(tmp[1])
	}
	print(k)
	return(NA)
}
addAlpha <- function(c) {
	op <- 0.7 
	return(paste(c,as.hexmode(as.integer(op*255)),sep=""))
}


oranges <- sapply(brewer.pal(8,"Oranges"),addAlpha)
findcolor <- function(v) {
	return(grey(v))
	#return(rgb2hex(jet.colorf(v)))
}
densn <- sapply(map@data$Kunta,comp)
dens <- cut(densn,c(0,1,5,10,20,50,100,1000,Inf))
# ndens <- log(dens)
# ndens <- ndens - min(ndens)
# ndens <- ndens/max(ndens)
# map@data <- cbind(map@data,dens=dens,logdens=ndens)
map@data <- cbind(map@data,dens=dens)
#map@data <- cbind(map@data,logdens=log(dens))

obs <-read.csv("BearObs2010XY_ETRS.csv",header=FALSE,sep=" ",skip=0,col.names=c("x","y"))
# coordinates(obs) <- ~x+y
# lx <- max(obs@coords[,1])-min(obs@coords[,1])
# ly <- max(obs@coords[,2])-min(obs@coords[,2])
# asp <- lx/ly
# obs@coords[,1] <- (obs@coords[,1]-min(obs@coords[,1]))/lx 
# obs@coords[,2] <- (obs@coords[,2]-min(obs@coords[,2]))/lx 


rajj <- list("sp.lines", valtraj, col = "#0000ff50")

#cols <- sapply(ndens,findcolor)

# don't display the frame

pt <- list(
		"b2010"=list(
			"sp.points", 
			obs, 
			pch = 16, 
			col = "#00000088", 
			cex=0.5),
		"b2009"=list("sp.points", 
			obs2009[obs2009$laji=="karhu",c("x","y")], 
			pch = 16, 
			col = "#00000088", 
			cex=0.5),
		"l2009"=list("sp.points", 
			obs2009[obs2009$laji=="ilves",c("x","y")], 
			pch = 16, 
			col = "#00000088", 
			cex=0.5),
		"w2009"=list("sp.points", 
				obs2009[obs2009$laji=="susi",c("x","y")], 
				pch = 16, 
				col = "#00000088", 
				cex=0.5)
		)


pl_b10 <- mapply(listplot,
		names(pt),
		rep(list(map),4),
		sp.layout=mapply(list,pt,rep(list(rajj),4),SIMPLIFY=FALSE),
		MoreArgs=list(
			plotfn=spplot,
			file="report/%s.pdf",
			width=3,
			height=5.5,
			afterfn=function(p,k) {
				trellis.par.set(axis.line=list(col=NA))
				print(p)
			},
			par.settings=theme.novpadding,
			mar=mar_tight,
			main="",
			zcol="dens",
			col.regions=oranges,
			col="transparent",
			colorkey=FALSE,
			pretty=TRUE
		))


# pl_b09 <- spplot(
# 			map,
# 			main="Bear - 2009",
# 			zcol="dens",
# 			sp.layout=list(pts2,rajj),
# 			col.regions=oranges,
# 			col="transparent",
# 			colorkey=FALSE,
# 			pretty=TRUE)

# pt_w09 <- list("sp.points", 
# 			obs2009[obs2009$laji=="susi",c("x","y")], 
# 			pch = 16, 
# 			col = "#00000088", 
# 			cex=0.5)

# pl_w09 <- spplot(
# 			map,
# 			main="Wolf - 2009",
# 			zcol="dens",
# 			sp.layout=list(pt_w09,rajj),
# 			col.regions=oranges,
# 			colorkey=FALSE,
# 			col="transparent",
# 			pretty=TRUE)
# pt_l09 <- list("sp.points", 
# 			obs2009[obs2009$laji=="ilves",c("x","y")], 
# 			pch = 16, 
# 			col = "#00000088", 
# 			cex=0.5)

# pl_l09 <- spplot(
# 			map,
# 			main="Lynx - 2009",
# 			zcol="dens",
# 			sp.layout=list(pt_l09,rajj),
# 			col.regions=oranges,
# 			colorkey=FALSE,
# 			col="transparent",
# 			pretty=TRUE)

#print(pl_b09)
#quartz()
print(pl_b10)
#quartz()
#print(pl_w09)
#quartz()
#print(pl_l09)
#plot(obs, pch = 4, col = "black", cex=0.5,add=TRUE)

