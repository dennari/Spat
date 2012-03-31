
exportFigs <- 0
displayFigs <- 1

plots.obs <- F
plots.mesh <- F
plots.results <- T

run.mesh <- T
run.INLA <- T



source("util.R",local=TRUE)
if(1 | !exists("kunnat")) {
	source("loadData.R",local=TRUE)	
	source("dist.R",local=TRUE)
	ocoords <- coords
}
bb <- map@bbox
#corners <- matrix(c(bb[1,1],bb[2,1], bb[1,2],bb[2,1], bb[1,2],bb[2,2], bb[1,1],bb[2,2] ),4,2,byrow=TRUE)
#ocoords <- corners[dim(corners)[1]:1,] 
#coords <- ocoords
lx <- max(ocoords[,1])-min(ocoords[,1])
ly <- max(ocoords[,2])-min(ocoords[,2])

offX <- mean(ocoords[,1])
offY <- mean(ocoords[,2])

# coords[,1] <- (ocoords[,1]-offX)/ly
# coords[,2] <- (ocoords[,2]-offY)/ly

if(run.mesh) {

	N <- nrow(coords)
	loc.bnd <- coords[N:1,]
	segm.bnd <- inla.mesh.segment(as.matrix(loc.bnd))

	maxedge <- ly/50
	minedge <- maxedge/1.11  # 20544.2 works for both
	mesh <- inla.mesh.create(
	                  ## Data locations:
	                 bound=segm.bnd,
	                  ## Where to put the mesh/graph files:
	                  ## Set to >=0 for visual (not on Windows):
	                  plot.delay=NULL,
	                  keep=FALSE,
	                  cutoff=minedge,
	                  refine=(list(min.angle=26,
	                               max.edge.data=maxedge,
	                               max.edge.extra=2*maxedge))
	        )
	spde <- inla.spde.create(mesh, model="matern")
}

##########################
#### OBSERVATIONS ########
##########################

loc <- as.matrix(obs.y2010)
#loc <- loc[!duplicated(loc[,c("x","y")]),]
#loc <- loc[sample.int(min(2000,nrow(loc))),]
# # discard points that fall outside the window
bnd.polygon <- mesh$loc[mesh$segm$bnd$idx[,1],1:2]
bnd.polygon <- bnd.polygon[!duplicated(bnd.polygon[,1:2]),]

loc <- pip(loc,bnd.polygon,bound=FALSE)





myplot(bnd.polygon,
	plotfn=lines,
	width=4,
	height=ly/lx*4,
	xrange=range(mesh$loc[,1]),
	yrange=range(mesh$loc[,2])
)
myplot(mesh$segm$bnd, mesh$loc,
	plotfn=lines,
	width=4,
	height=ly/lx*4,
	xrange=range(mesh$loc[,1]),
	yrange=range(mesh$loc[,2])
)
#pp <- ppp(mesh$loc[,1],mesh$loc[,2],xrange=range(mesh$loc[,1]),yrange=range(mesh$loc[,2]))
#loc[t(cbind(segm$idx[idx,, drop=FALSE], NA)), 1]
##############################
#### COVARIATES ##############
##############################
#distmat <- crossdist(mesh$loc[,1],mesh$loc[,2],cov.pdens[,"x"],cov.pdens[,"y"])
#pp <- ppp(mesh$loc[,1],mesh$loc[,2],poly=bnd.polygon,marks=cov.pdens)
#pp1 <- ppp(loc[,1],loc[,2],poly=bnd.polygon)
pp2 <- ppp(cov.pdens[,"x"],cov.pdens[,"y"],poly=bnd.polygon,marks=log(cov.pdens[,"pdens"]))
contdens <- markmean(pp2,dimyx=c(400,400))
msk <- as.mask(pp2,dimyx=c(400,400))
pxl <- nearest.raster.point(cov.pdens[,"x"],cov.pdens[,"y"],msk,indices=T)
vv <- contdens$v[pxl$row,pxl$col]

#distmat <- crossdist(mesh$loc[,1],mesh$loc[,2],contdens$xcol,contdens$yrow)

distmat <- crossdist(mesh$loc[,1],mesh$loc[,2],cov.pdens[,"x"],cov.pdens[,"y"])


cov.mesh_pdens <- cov.pdens[apply(distmat,1,which.min),"pdens"]
#loc <- loc[loc[,1] < 0.17,]

# kappa =0.5
# tau =0.5

############################
#### RESCALE ###############
############################

# loc[,1] <- (loc[,1]-offX)/ly
# loc[,2] <- (loc[,2]-offY)/ly

# ##image.plot(inla.mesh.project(proj, variances^0.5), col=cp(256))


# sample = inla.spde.query(spde, sample=list(tau=tau, kappa2=kappa^2))$sample

# max_points = 10000
# intercept =log(max_points) - max(sample)
# source("Daniel/simulate_lgcp.R",local=TRUE)
# loc <- simulate_lgcp(mesh,intercept+sample)

# print(dim(loc)[1])

#plot.new()
myplot(mesh)
points(loc,col="#ff000080",pch=16)
nV <- mesh$n
nData <- dim(loc)[1]
print(nData)
print(nV)
# print(dim(loc)[1])

if(run.INLA) {
	## Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
	#Convenient definitions - number of vertices and number of observations

	#See paper - we need two A matrices, one for the location of the points and one
	#for the integration points (here taken to be the identity matrix)
	LocationMatrix = inla.mesh.project(mesh, loc)$A
	IntegrationMatrix = sparseMatrix(i=1:nV,j=1:nV,x=rep(1,nV)) #simple integration scheeme
	ObservationMatrix=rBind(IntegrationMatrix,LocationMatrix) #Glue matrices together

	#The integration weights (alpha in the paper)
	IntegrationWeights = diag(spde$internal$c0)
	E_point_process =c(IntegrationWeights,rep(0,nData)) #Glue everything together



	fake_data <- c(rep(0,nV),rep(1,nData)) #Put a zero where there aren't observations and a 1 where there is a point

	formula <- y ~ 1 + f(idx, model=spde) #Basic latent model - feel free to add covariates etc
	formulac <- y ~ 1 + pdens + f(idx, model=spde) #Basic latent model - feel free to add covariates etc


	data <- list(y=fake_data,idx = c(1:nV),pdens=cov.mesh_pdens) #put hte data in

	#The INLA call.  Likelihood is Poisson with Observation Matrix and appropriate value fo E.
	#result = inla(formula, data=data, family="poisson",
	#	control.predictor=list(A=ObservationMatrix),E=E_point_process,verbose=TRUE)

	sigma2.approx = 0.5^2
	range.approx = 0.4
	kappa = sqrt(8)/range.approx
	tau = 1/(4*pi*kappa^2*sigma2.approx)^0.5

	#init.mode = c(-5.5, 3) # works
	#init.mode = c(-8, 5.5) # works
	#init.mode = c(-4.905784,6.317846)
	#init.mode = c(log(tau), log(kappa^2))
	#init.mode = c(-1.304829,2.717723)

	#control.mode = list(theta = init.mode, restart=TRUE)

#The INLA call.  Likelihood is Poisson with Observation Matrix and appropriate value fo E.
result1 =
    inla(formula, data=data, family="poisson",
         control.predictor=list(A=ObservationMatrix,compute=TRUE),
         E=E_point_process,
         #control.mode = list(theta = init.mode, restart=TRUE),
         verbose=TRUE,
          ## We don't need the marginals:
          control.compute = list(return.marginals=FALSE),
#control.mode = control.mode,
          ## We don't need to overoptimise:
          control.inla=list(tolerance=1e-4)
         )

#  result2 =
#     inla(formulac , data=data, family="poisson",
#          control.predictor=list(A=ObservationMatrix,compute=TRUE),
#          E=E_point_process,
#          #control.mode = list(theta = init.mode, restart=TRUE),
#          verbose=TRUE,
#           ## We don't need the marginals:
#           control.compute = list(return.marginals=FALSE),
# #control.mode = control.mode,
#           ## We don't need to overoptimise:
#           control.inla=list(tolerance=1e-4)
#          )


}




#############################
### PLOTS ###################
#############################

if(displayFigs | exportFigs) {

	ms <-summary(mesh)

	wh <- (ms$xlim[2]-ms$xlim[1])/(ms$ylim[2]-ms$ylim[1])

if(plots.results) {
	#########################
	#Some fancy plotting code stolen from Finn
	###########################
	## Prepare for plotting:

	## Calculate mapping between triangulation vertices and grid points:

	proj = inla.mesh.projector(mesh, dims=c(200,200))

	## Construct greyscale palette function:
	my.grey.palette = function (n,...) { return (grey.colors(n,0.05,0.95,...))}
	## Use it:
	my.palette = my.grey.palette

  	cp = colorRampPalette(c("darkblue", "blue", "cyan", "yellow", "red", "darkred"))
  	my.palette = cp

	## Construct map data appropriate for easy plotting:
	#mm = calc.map(map)

	levelplotmap = function(..., mm) {
	    panel.levelplot(...)
	    panel.lines(mm$x, mm$y, col="red")
	}


	#####################
	## Plot results:


	## Map resulting posterior mean field to a grid:
	plotdata = inla.mesh.project(proj, result$summary.random$idx[,"mean"])
	## Plot PM contours:
	if(exportFigs) {
		pdf(file="report/postint_b2010.pdf",width=5,height=5/wh)
		par(mar=mar_tight)
	} else {

		dev.new()
	}
	bbb = (levelplot(row.values=proj$x, column.values=proj$y, x=plotdata,
	                 mm=NULL, panel=levelplotmap,
	                 col.regions=my.palette,
	                 xlim=range(proj$x), ylim=range(proj$y), aspect="iso",
	                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
	                 xlab="Easting",ylab="Northing"))
	print(bbb)
	if(exportFigs) {
		dev.off()
	}

	# pp<-myplot(plotdata,
	# 	file="report/postint_b2010.pdf",
	# 	plotfn=levelplot,
	# 	width=5,
	# 	height=5/wh,
	# 	afterfn=function(p,k) {
	# 		print(p)
	# 	},
	# 	mar=mar_tight,
	# 	main="",
	# 	col="#00000020",
	# 	row.values=proj$x, 
	# 	column.values=proj$y,
	#     mm=NULL, 
	#     panel=levelplotmap,
	#     col.regions=my.palette,
	#     xlim=range(proj$x), ylim=range(proj$y), aspect="iso",
	#     contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
	#   	lab="Easting",ylab="Northing"
	# )



	#Plot the interecept and the hyperparameter posteriors
	#plot(result, plot.random.effects=FALSE)

	#Plot the posterior for the RANGE of the random field
	range=inla.tmarginal(function(x) sqrt(8)/exp(x/2),(result$marginals.hyperpar$"K.0 for idx-basisK" ))
	attr(range,"inla.tag")="Range"
	#plot(range)



}



if(plots.mesh) {

	ms <-summary(mesh)

	wh <- (ms$xlim[2]-ms$xlim[1])/(ms$ylim[2]-ms$ylim[1])
	myplot(mesh,
		file="report/mesh.pdf",
		width=5,
		height=5/wh,
		mar=mar_tight,
		main="",
		sub="",
		col="#00000020"
	)

}

if(plots.obs) {

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

	rajj <- list("sp.lines", map, col = "#0000ff70")
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
trellis.par.set(theme = theme.novpadding)

lw <- list(left.padding = list(x = 0, units = "inches"))
lw$right.padding <- list(x = -0.1, units = "inches")
lh <- list(bottom.padding = list(x = 0, units = "inches"))
lh$top.padding <- list(x = -0.2, units = "inches")

lattice.options(layout.widths = lw, layout.heights = lh)
	# points and window
	myplot(mesh$segm$bnd, mesh$loc,
		plotfn=lines,
		width=4,
		height=ly/lx*4,
		xrange=range(mesh$loc[,1]),
		yrange=range(mesh$loc[,2]),
		mar=mar_tight,
		file="report/b2010.pdf",
		xlab=NULL,ylab=NULL,
		afterfn=function(p,k) {
			points(loc,col="#00000080",pch=16)
		}
	)

	pp2 <- ppp(cov.pdens[,"x"],cov.pdens[,"y"],poly=bnd.polygon,marks=log(cov.pdens[,"pdens"]))
	contdens <- markmean(pp2,dimyx=c(400,400))
	pl_b10 <- myplot(contdens,
			col=jet.colors(512),
			main="",
			mar=mar_tight,
			ribbon=FALSE,
			file="report/pdens.pdf",
			width=4,
			height=ly/lx*4,
			afterfn=function(p,k) {
				lines(mesh$segm$bnd, mesh$loc)
				print(p)
			}
		)

}
}

