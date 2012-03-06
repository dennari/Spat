###########################################################################
## oceans.R
## Runs the point process on oceans example from Section 8.3 of Simpson et al. (2011)
## Author: Daniel Simpson
## This code comes with no warranty or guarantee of any kind.
###########################################################################

source("simulate_LGCP_ocean.R")

#Convenient definitions - number of vertices and number of observations
nV=mesh$n
nData <- dim(pts)[1]
loc =pts

#Make the spde - simple non-intrinsic Matern
spde = inla.spde.create(mesh, model="matern")


#See paper - we need two A matrices, one for the location of the points and one for the integration points (here taken to be the identity matrix)
LocationMatrix = inla.mesh.project(mesh, loc)$A
IntegrationMatrix = sparseMatrix(i=1:nV,j=1:nV,x=rep(1,nV)) #simple integration scheeme
ObservationMatrix=rBind(IntegrationMatrix,LocationMatrix) #Glue matrices together

#The integration weights (alpha in the paper)
IntegrationWeights = diag(spde$internal$c0)
E_point_process =c(IntegrationWeights,rep(0,nData)) #Glue everything together





fake_data = c(rep(0,nV),rep(1,nData), rep(NA,nV)) #Put a zero where there aren't observations and a 1 where there is a point

formula = y ~ 1 + f(idx, model=spde) #Basic latent model - feel free to add covariates etc


data = list(y=fake_data,idx = c(1:nV)) #put hte data in

init.mode = c(log(tau), log(kappa^2))
control.mode = list(theta = init.mode, restart=TRUE)

#The INLA call.  Likelihood is Poisson with Observation Matrix and appropriate value fo E.
result =
    inla(formula, data=data, family="poisson",
         control.predictor=list(A=rBind(ObservationMatrix,IntegrationMatrix), compute=TRUE),
         E=c(E_point_process,IntegrationWeights),
         control.mode = list(theta = init.mode, restart=TRUE),
         verbose=TRUE
         )





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

## Construct map data appropriate for easy plotting:
#mm = calc.map(map)

levelplotmap = function(..., mm) {
    panel.levelplot(...)
    panel.lines(mm$x, mm$y, col="red")
}


#####################
## Plot results:

if (FALSE){
## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj, result$summary.random$idx[,"mean"]+result$summary.fixed[1,"mean"])
## Plot PM contours:
dev.new()
bbb = (levelplot(row.values=proj$x, column.values=proj$y, x=plotdata,
                 mm=NULL, panel=levelplotmap,
                 col.regions=my.palette,
                 xlim=range(proj$x), ylim=range(proj$y), aspect="iso",
                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                 xlab="Easting",ylab="Northing"))
print(bbb)

#Plot the interecept and the hyperparameter posteriors
plot(result, plot.random.effects=FALSE)

#Plot the posterior for the RANGE of the random field
range=inla.tmarginal(function(x) sqrt(8)/exp(x/2),(result$marginals.hyperpar$"K.0 for idx-basisK" ))
attr(range,"inla.tag")="Range"
plot(range)
}


source("ocean_figs.R")
