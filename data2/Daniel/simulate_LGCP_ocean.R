###########################################################################
## simulate_LGCP_ocean
## Simulates a LGCP  over the oceans. See Section 8.3 of Simpson et al. (2011)
## Author: Daniel Simpson
## This code comes with no warranty or guarantee of any kind.
###########################################################################

require(INLA)
require(rgl)
require(fields)

source("ocean_mesh.R")
source("simulate_lgcp.R")


max_points = 10000 #Maximum number of pts you want to simulate


spde=inla.spde.create(mesh, model="matern")

omesh = old.mesh.class(mesh)
cp=colorRampPalette(c("blue","cyan","yellow","red"))
proj = inla.mesh.projector(mesh,dims=c(180,90),projection="longsinlat")

##hemi = as.numeric(mesh$loc[,3]>0)

sigma2.approx = 0.5^2
range.approx = 0.4
kappa = sqrt(8)/range.approx
tau = 1/(4*pi*kappa^2*sigma2.approx)^0.5


##image.plot(inla.mesh.project(proj, variances^0.5), col=cp(256))


sample = inla.spde.query(spde, sample=list(tau=tau, kappa2=kappa^2))$sample
##sample = sample-sqrt(1-mesh$loc[,3]^2)*0.25
##plotdata = inla.mesh.project(proj, sample)

##plot(omesh, sample, color.palette=cp, draw.edges=FALSE, draw.vertices=FALSE)

##image.plot(plotdata, col=cp(256))

intercept =log(max_points) - log(4*pi) - 2*log(radius.of.earth)- max(sample)

pts = simulate_lgcp(mesh, intercept + sample)

##points3d(pts)

