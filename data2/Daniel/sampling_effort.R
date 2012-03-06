###########################################################################
## sampling_effort.R
## This code runs the example for a point process with a hole in the sampling
## effort from Section 8.2 of Simpson et al. (2011)
## Author: Daniel Simpson
## This code comes with no warranty or guarantee of any kind.
###########################################################################

require(INLA)
source("simulate_lgcp.R")
source("utils.R")

set.seed(3216)
DEGRADED=TRUE
EXACT=TRUE
BAD=TRUE
do.print=FALSE
#make a mesh
boundary = inla.mesh.segment(loc=matrix(c(-1,-1,-1,1,1,1,1,-1),byrow=T,ncol=2)[4:1,])
mesh = inla.mesh.create(boundary=boundary,refine=list(max.edge=0.05))

#sample the point process
max_points = 1000 #Maximum number of pts you want to simulate

kappa = 15
tau =0.1


spde=inla.spde.create(mesh, model="matern")
sample = inla.spde.query(spde, sample=list(tau=tau, kappa2=kappa^2))$sample
intercept =log(max_points) - log(4) - max(sample)
pts = simulate_lgcp(mesh, intercept + sample)

#degrade sample - remove anything in the rectange [xL,xR]x[yL,yR]
xL =-0.5
xR =0.4
yL =-0.1
yR =0.4


remove <- (pts[,1] < xR)& (pts[,1]>xL)&(pts[,2]<yR)&(pts[,2]>yL)
loc = pts[!remove,]
dev.new()
plot(pts)
points(pts[remove,],col="red")
lines(matrix(c(xL,yL,xL,yR,xR,yR,xR,yL,xL,yL),ncol=2,byrow=T),col="red")

if(do.print){

dev.copy2pdf(file="../Point Processs with SPDEs/hole_data.pdf")


}


########################
#Fit true pattern
#########################

if(EXACT){
nV=mesh$n
nData <- dim(pts)[1]

#See paper - we need two A matrices, one for the location of the points and one for the integration points (here taken to be the identity matrix)
LocationMatrix = inla.mesh.project(mesh, pts)$A
IntegrationMatrix = sparseMatrix(i=1:nV,j=1:nV,x=rep(1,nV)) #simple integration scheeme
ObservationMatrix=rBind(IntegrationMatrix,LocationMatrix) #Glue matrices together

#The integration weights (alpha in the paper)
IntegrationWeights = diag(spde$internal$c0)
E_point_process =c(IntegrationWeights,rep(0,nData)) #Glue everything together



fake_data = c(rep(0,nV),rep(1,nData)) #Put a zero where there aren't observations and a 1 where there is a point

formula = y ~ 1 + f(idx, model=spde) #Basic latent model - feel free to add covariates etc


data = list(y=fake_data,idx = c(1:nV)) #put hte data in

#The INLA call.  Likelihood is Poisson with Observation Matrix and appropriate value fo E.
result.true = inla(formula, data=data, family="poisson",control.predictor=list(A=ObservationMatrix),E=E_point_process,verbose=TRUE)


}




########################
#Fit degraded
########################
if(DEGRADED) {
bndIN = inla.mesh.segment(loc=matrix(c(xL,yL,xL,yR,xR,yR,xR,yL,xL,yL),ncol=2,byrow=T))
bndOUT = inla.mesh.segment(loc=matrix(c(-1,-1,-1,1,1,1,1,-1),byrow=T,ncol=2)[4:1,])
mesh0 = inla.mesh.create(boundary=list(bndIN,bndOUT),refine=list(max.edge=0.05))

mesh.degraded = inla.mesh.create(loc=mesh0$loc[,1:2],,boundary=bndOUT,refine=list(max.edge=0.1))



if (do.print) {
  dev.new()
  plot(mesh.degraded)

  dev.copy2pdf(file="../Point Processs with SPDEs/hole_mesh.pdf")

}


spde.degraded = inla.spde.create(mesh.degraded, model="matern")


nV=mesh.degraded$n
nData <- dim(loc)[1]


#The integration weights (alpha in the paper)
###################
### THIS IS THE DIFFERENT PART
### Work out which bits are actually sampled and modify the E parameter appropriately
###################
meshLoc <- mesh.degraded$loc
sampled  <- !( ((meshLoc[,1] >xL)&(meshLoc[,1]<xR))&((meshLoc[,2]>yL)&(meshLoc[,2]<yR)))

if (FALSE) {
## THE OLD WAY - THIS WORKS! and it's slightly faster than what's in the paper (the extra points hole doesn't turn up in the likelihood at all!  But it's more complicated to explain and only works when things are unobserved (not when there is variable sampling effort), so we haven't done it this way.


#See paper - we need two A matrices, one for the location of the points and one for the integration points (here taken to be the identity matrix)
LocationMatrix = inla.mesh.project(mesh.degraded, loc)$A
IntegrationMatrix = sparseMatrix(i=1:nV,j=1:nV,x=rep(1,nV)) #simple integration scheeme
ObservationMatrix=rBind(IntegrationMatrix[sampled,],LocationMatrix) #Glue matrices together


IntegrationWeights = diag(spde.degraded$internal$c0)[sampled]
E_point_process =c(IntegrationWeights,rep(0,nData)) #Glue everything together
}

if (TRUE) {
sampled = as.numeric(sampled)
LocationMatrix = inla.mesh.project(mesh.degraded, loc)$A
IntegrationMatrix = sparseMatrix(i=1:nV,j=1:nV,x=rep(1,nV)) #simple integration scheeme
ObservationMatrix=rBind(IntegrationMatrix,LocationMatrix) #Glue matrices together
IntegrationWeights = diag(spde.degraded$internal$c0)
E_point_process =c(sampled*IntegrationWeights,rep(0,nData)) #Glue everything together
}


fake_data = c(rep(0,nV),rep(1,nData)) #Put a zero where there aren't observations and a 1 where there is a point

formula = y ~ 1 + f(idx, model=spde.degraded) #Basic latent model - feel free to add covariates etc


data = list(y=fake_data,idx = c(1:nV)) #put hte data in

#The INLA call.  Likelihood is Poisson with Observation Matrix and appropriate value fo E.
result.degraded = inla(formula, data=data, family="poisson",control.predictor=list(A=ObservationMatrix),E=E_point_process,verbose=TRUE)
#,control.mode= list(theta= result.true$mode$theta,restart=TRUE),control.inla = list(cmin=0.05) )

#result.degraded2= inla(formula, data=data, family="poisson",control.predictor=list(A=ObservationMatrix),E=E_point_process,verbose=TRUE,control.mode= list(result=result.degraded,restart=TRUE),control.inla = list(cmin=0.005) )

}



########################
#Fit degraded model with normal (i.e. unadapted) mesh
#########################

if(BAD){
nV=mesh$n
nData <- dim(loc)[1]


#The integration weights (alpha in the paper)
###################
### THIS IS THE DIFFERENT PART
### Work out which bits are actually sampled and modify the E parameter appropriately
###################
meshLoc <- mesh$loc
sampled  <- !( ((meshLoc[,1] >xL)&(meshLoc[,1]<xR))&((meshLoc[,2]>yL)&(meshLoc[,2]<yR)))

#sampled[sampled==0] = 1e-6

#See paper - we need two A matrices, one for the location of the points and one for the integration points (here taken to be the identity matrix)
LocationMatrix = inla.mesh.project(mesh, loc)$A
IntegrationMatrix = sparseMatrix(i=1:nV,j=1:nV,x=rep(1,nV)) #simple integration scheeme
ObservationMatrix=rBind(IntegrationMatrix[sampled,],LocationMatrix) #Glue matrices together


IntegrationWeights = diag(spde$internal$c0)[sampled]
E_point_process =c(IntegrationWeights,rep(0,nData)) #Glue everything together


fake_data = c(rep(0,length(which(sampled))),rep(1,nData)) #Put a zero where there aren't observations and a 1 where there is a point

formula = y ~ 1 + f(idx, model=spde) #Basic latent model - feel free to add covariates etc


data = list(y=fake_data,idx = c(1:nV)) #put hte data in

#The INLA call.  Likelihood is Poisson with Observation Matrix and appropriate value fo E.
result.bad = inla(formula, data=data, family="poisson",control.predictor=list(A=ObservationMatrix),E=E_point_process,verbose=TRUE)


}




#########################
#Some fancy plotting code stolen from Finn



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

###########################
## Prepare for plotting:

## Calculate mapping between triangulation vertices and grid points:

proj.degraded = inla.mesh.projector(mesh.degraded, dims=c(200,200))


#####################
## Plot results:

if (DEGRADED) {

## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj.degraded, result.degraded$summary.random$idx[,"mean"])
## Plot PM contours:
dev.new()
bbb = (levelplot(row.values=proj.degraded$x, column.values=proj.degraded$y, x=plotdata,
                 mm=NULL,
                 col.regions=my.palette,
                 xlim=range(proj.degraded$x), ylim=range(proj.degraded$y), aspect="iso",
                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                 xlab="Easting",ylab="Northing"))
print(bbb)

if (do.print) {
dev.copy2pdf(file="../Point Processs with SPDEs/degraded_mean.pdf")
}
}

if (EXACT) { 
#########################
#Some fancy plotting code stolen from Finn
###########################
## Prepare for plotting:

## Calculate mapping between triangulation vertices and grid points:

proj = inla.mesh.projector(mesh, dims=c(200,200))


#####################
## Plot results:


## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj, result.true$summary.random$idx[,"mean"])
## Plot PM contours:
dev.new()
bbb = (levelplot(row.values=proj$x, column.values=proj$y, x=plotdata,
                 mm=NULL, panel=levelplotmap,
                 col.regions=my.palette,
                 xlim=range(proj$x), ylim=range(proj$y), aspect="iso",
                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                 xlab="Easting",ylab="Northing"))
print(bbb)
if (do.print)
  {
dev.copy2pdf(file="../Point Processs with SPDEs/true_mean.pdf")


  }
}


if (do.print) {
dev.new()
plot(result.degraded$marginals.hyperpar$"T.0 for idx-basisT",lty=2,type="l")
lines(result.true$marginals.hyperpar$"T.0 for idx-basisT",lty=1)
lines(result.bad$marginals.hyperpar$"T.0 for idx-basisT",lty=3)
dev.copy2pdf(file="../Point Processs with SPDEs/T0_compare.pdf")

dev.new()
plot(result.degraded$marginals.hyperpar$"K.0 for idx-basisK",lty=2,type="l")
lines(result.true$marginals.hyperpar$"K.0 for idx-basisK",lty=1)
lines(result.bad$marginals.hyperpar$"K.0 for idx-basisK",lty=3)
dev.copy2pdf(file="../Point Processs with SPDEs/K0_compare.pdf")

}
