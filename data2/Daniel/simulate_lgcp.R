###########################################################################
## simulate_lgcp.R
## Code for simulating LGCPs on triangulated domains (planar or spherical)
## used throughout Simpson et al. (2011)
## Author: Daniel Simpson
## This code comes with no warranty or guarantee of any kind.
###########################################################################

require(INLA)

simulate_lgcp <- function(mesh,weights) {

if(mesh$manifold == "R2") {
  print("It's flat!")
#Construct bounding rectangle
loc <- mesh$loc
xmin = min(loc[,1])
xmax = max(loc[,1])
ymin = min(loc[,2])
ymax = max(loc[,2])
area =(xmax- xmin)*(ymax- ymin)

#Simulate number of points
lambda_max <- max(weights)
Npoints <- rpois(n=1,lambda=area*exp(lambda_max))

  if (Npoints > 1e5) {
    print(Npoints)
    print("You've got to be joking")
    return(NULL)
  }

#Simulate uniform points on the bounding rectangle
x <- runif(n=Npoints, min=xmin, max=xmax)
y <- runif(n=Npoints, min=ymin, max=ymax)

points <-cbind(x,y)

#Do some thinning
A <- inla.mesh.project(mesh,points)$A

weights =exp( weights-lambda_max)

pointValues = A%*%weights

keep =which(runif(Npoints) < pointValues)

return(points[keep,])
}

if (mesh$manifold=="S2") {
  print("Ground control to Major Tom")
area =   4*pi*radius.of.earth^2 #we're on the earth

#Simulate number of points
lambda_max <- max(weights)
Npoints <- rpois(n=1,lambda=area*exp(lambda_max))
print(area*exp(lambda_max))
print(Npoints)

  if (Npoints > 1e5) {
    print(Npoints)
    print("You've got to be joking")
    return(NULL)
  }
  
#Simulate random points on the sphere
#this could be faster
points= c()
for (i in 1:Npoints)
  {
    tmp = rnorm(3)
    points = rbind(points,tmp/sqrt(t(tmp)%*%tmp))
  }

#Do some thinning
A <- inla.mesh.project(mesh,points)$A

weights =exp( weights-lambda_max)

pointValues = A%*%weights

keep =which(runif(Npoints) < pointValues)

return(points[keep,])

}

}
