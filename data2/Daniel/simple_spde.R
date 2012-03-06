require(INLA)
require(rgl) #For plotting!
require(boot) #for inverse logit

# Make some fake data
N=500

Ntrials = rep(100,N)
points = matrix(runif(2*N),ncol=2)

#The 'true' surface

fcn = function(loc){
	return( cos(loc[,1])*sin(4*loc[,2]))
	}
	
noisy_data = rbinom(N, Ntrials,inv.logit(fcn(points)))

## Build a mesh

bnd = inla.mesh.segment(matrix(c(0,0,1,0,1,1,0,1),ncol=2,byrow=TRUE))

mesh = inla.mesh.create(points,boundary=bnd,refine=list(max.edge=0.1))

plot(mesh)

##Make the SPDE

spde = inla.spde.create(mesh,model="imatern")

## Fit the model

fake_data = list(y=noisy_data, data_points = mesh$idx$loc)

formula = y ~ f(data_points, model=spde)-1

r = inla(formula, family="binomial",Ntrials= Ntrials, data=fake_data,verbose=TRUE)

## And some output
summary(r)
plot(r, plot.random.effects=FALSE)

plot(old.mesh.class(mesh),r$summary.random$data_points$mean)

#Calculate the MSE

print(sqrt(var(r$summary.random$data_points$mean - fcn(mesh$loc))))