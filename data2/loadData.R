kunnat <- readShapeSpatial("map45/HallintoAlue45Milj_ETRS.shp")
#kraj <- readShapeSpatial("map45/HallintoalueRaja45Milj_ETRS.shp",verbose=TRUE)
#kraj <- kraj[kraj@data$Kohdeluokk==82100|kraj@data$Kohdeluokk==84111,]
#city <- readShapeSpatial("map/cityp.shp")
# turn the factor to numeric vector
kunnat@data$Kunta <- as.numeric(as.vector(kunnat@data$Kunta))
hcont <- readShapeSpatial("map45/hcont_l.shp")
map <- hcont[row.names(hcont@data)=="761"|hcont@data$Z==9999,]
#forest <- readShapeSpatial("map/forest.shp")
population <- read.csv(
	"asukastiheys.csv",
	sep=",",
	header=FALSE,
	skip=0,
	col.names=c("Kunta",
				"land_area",
				"water_area",
				"water_area2",
				"total_area",
				"population",
				"pdens_area",
				"pdens_land_area"))


pdens <- sapply(kunnat@data$Kunta, function(k) {
	if(sum(population$Kunta==k) > 0) {
		tmp <- population[population$Kunta==k,c("pdens_land_area")]
		return(tmp[1])
	}
	print(k)
	return(NA)
})

kunnat_cnt <- coordinates(kunnat)

cov.pdens <- as.data.frame(cbind(
				kunta=kunnat@data$Kunta,
				x=kunnat_cnt[,1],
				y=kunnat_cnt[,1],
				pdens=pdens))

pdensbins <- c(0,1,5,10,20,50,100,1000,Inf)
kunnat@data <- cbind(kunnat@data,pdens=pdens,pdens_binned=cut(pdens,pdensbins))

obs.y2010 <-read.csv("BearObs2010XY_ETRS.csv",
	header=FALSE,
	sep=" ",
	skip=0,
	col.names=c("x","y"))
obs.y2009 <- read.csv(	
	"Obs2009XYS_ETRS.csv",
	sep=",",
	header=FALSE,
	skip=0,
	col.names=c("x","y","laji"))



