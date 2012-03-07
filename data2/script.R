
exportFigs <- 0
displayFigs <- 0
plots.obs <- FALSE

source("util.R",local=TRUE)
source("dist.R",local=TRUE)

#raj <- readShapeSpatial("map/HallintoalueRaja1Milj_ETRS.shp")
if(!exists("kunnat"))
	source("loadData.R",local=TRUE)	



# ## Build triangular mesh:
# mesh = (inla.mesh.create(
#                   ## Data locations:
#                   obs.y2010[1:300,],
#                   ## Where to put the mesh/graph files:
#                   ## Set to >=0 for visual (not on Windows):
#                   plot.delay=0,
#                   ## Encapsulate data region with a relative margin:
#                   extend=list(n=8, offset=-0.15),
#                   ## Refined triangulation,
#                   ## minimal angles >=26 degrees,
#                   ## interior maximal edge lengths 0.08,
#                   ## exterior maximal edge lengths 0.2:
#                   refine=(list(min.angle=26,
#                                max.edge.data=0.08,
#                                max.edge.extra=0.2))
#                   ))## Build triangular mesh:
#hel <- kunnat[6,]
N <- nrow(coords)
loc.bnd <- coords[N:1,]
#loc.bnd <- loc.bnd[1:50,]
#loc.bnd <- rbind(loc.bnd,loc.bnd[1,])
#loc.bnd[,1] <- (loc.bnd[,1]-min(loc.bnd[,1]))/(max(loc.bnd[,1])-min(loc.bnd[,1]))
#loc.bnd[,2] <- (loc.bnd[,2]-min(loc.bnd[,2]))/(max(loc.bnd[,2])-min(loc.bnd[,2]))
#loc.bnd[3,1] <- -1*loc.bnd[3,1]
#maxedge <- sqrt(hel@polygons[[1]]@Polygons[[1]]@area)/10
segm.bnd <- inla.mesh.segment(as.matrix(loc.bnd))

# loc.bnd <- matrix(c(0,0, 1,0, 0.5,0.5, 0.4,0.7, 1,1, 0,1),6,2,byrow=TRUE)
# loc.bnd <- matrix(c(0,0, 0.415,1, 0.3,1),3,2,byrow=TRUE)
# # loc.bnd <- rbind(loc.bnd,loc.bnd[1,])
# segm.bnd <- inla.mesh.segment(loc.bnd)

#maxedge <- 0.4 
maxedge <- (map@bbox[1,2]-map@bbox[1,1])/50
#loc <- city[sample.int(1874,size=1500),]@coords


mesh <- inla.mesh.create(
                  ## Data locations:
                 boundary=segm.bnd,
                  ## Where to put the mesh/graph files:
                  ## Set to >=0 for visual (not on Windows):
                  plot.delay=NULL,
                  keep=FALSE,
                  refine=(list(min.angle=26,
                               max.edge.data=maxedge,
                               max.edge.extra=2*maxedge))
        )

## Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
spde = inla.spde.create(mesh, model="matern", param=list(alpha=2))

#Convenient definitions - number of vertices and number of observations
nV=mesh$n
nData <- dim(loc)[1]

#See paper - we need two A matrices, one for the location of the points and one
#for the integration points (here taken to be the identity matrix)
LocationMatrix = inla.mesh.project(mesh, loc)$A
IntegrationMatrix = sparseMatrix(i=1:nV,j=1:nV,x=rep(1,nV)) #simple integration scheeme
ObservationMatrix=rBind(IntegrationMatrix,LocationMatrix) #Glue matrices together





#############################
### PLOTS ###################
#############################

if(displayFigs | exportFigs) {

ms <-summary(mesh)

wh <- (ms$xlim[2]-ms$xlim[1])/(ms$ylim[2]-ms$ylim[1])
myplot(mesh,
	file="report/mesh.pdf",
	width=5,
	height=5/wh,
	afterfn=function(p,k) {
		lines(map,col="#0000ffff")
		print(p)
	},
	mar=mar_tight,
	main="",
	col="#00000020"
)




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

oranges <- sapply(brewer.pal(8,"Oranges"),function(c) {addAlpha(c,0.7)})

rajj <- list("sp.lines", map, col = "#0000ff50")
pt <- list(
		"b2010"=list(
			"sp.points", 
			obs.y2010, 
			pch = 16, 
			col = "#00000088", 
			cex=0.5),
		"b2009"=list("sp.points", 
			obs.y2009[obs.y2009$laji=="karhu",c("x","y")], 
			pch = 16, 
			col = "#00000088", 
			cex=0.5),
		"l2009"=list("sp.points", 
			obs.y2009[obs.y2009$laji=="ilves",c("x","y")], 
			pch = 16, 
			col = "#00000088", 
			cex=0.5),
		"w2009"=list("sp.points", 
				obs.y2009[obs.y2009$laji=="susi",c("x","y")], 
				pch = 16, 
				col = "#00000088", 
				cex=0.5)
		)

if(plots.obs) {
pl_b10 <- mapply(listplot,
		names(pt),
		rep(list(kunnat),4),
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
			zcol="pdens_binned",
			col.regions=oranges,
			col="transparent",
			colorkey=FALSE,
			pretty=TRUE
		))
}
}

