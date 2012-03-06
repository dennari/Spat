require(INLA)
require(maps) #install.packages(c("maps"))
load("anom1962.RData")

#Work out the size of the rectangel
xrange_data = range(loc[,1]) 
yrange_data = range(loc[,2])

xlength = diff(xrange_data)
ylength = diff(yrange_data)


#Add a buffer to account for edge effects
xrange = xrange_data + c(-xlength/5,xlength/5)
yrange = yrange_data + c(-ylength/5,ylength/5)


# approximately one knot every two degrees.
Nx = 100 
Ny = 50

# Findthe x and y knot points
xpoints = seq(xrange[1], xrange[2], length.out=Nx)
ypoints = seq(yrange[1], yrange[2], length.out=Ny)


#Tell fmesher that it's a lattice! (convenient if you have both lattice and non-lattice points and you want to separate them in the output)
lattice = inla.mesh.lattice(xpoints,ypoints)

#Make a mesh
mesh = inla.mesh.create(lattice=lattice, extend=FALSE)

#mesh = inla.mesh.create(loc, extend= TRUE,cutoff=2, refine = list(max.edge=2, min.angle=24)) #mesh to data points. Extend the domain, treat points that are closer than two units as a single point, make the largest edge 2 units, and make the smallest angle at least 24 degrees.



#Make an spde model object
spde = inla.spde.create(mesh,model="matern")


#Organise the data
data = list(z=z, loc = seq(1,mesh$n))

#The data points and the observation locations don't line up---construct an observation matrix
ObservationMatrix= inla.mesh.project(mesh,loc)$A


#make a formula
formula = z ~ 1 + f(loc, model =spde)

#do the inference - Observation matrix incorporated into control.predictor
result = inla(formula,family = "gaussian", data = data, control.predictor = list(A = ObservationMatrix), verbose= TRUE)


#The posterior mean
z = result$summary.random$loc$mean

#Do an easy plot.
plot(old.mesh.class(mesh), z)



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
    panel.lines(mm$x, mm$y, col="black")
}



#####################
## Plot results:


## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj, result$summary.random$loc[,"mean"])
## Plot PM contours:
dev.new()
bbb = (levelplot(row.values=proj$x, column.values=proj$y, x=plotdata,
                 mm=map.poly("usa"), panel=levelplotmap,
                 col.regions=my.palette,
                 xlim=range(proj$x), ylim=range(proj$y), aspect="iso",
                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                 xlab="Easting",ylab="Northing"))
print(bbb)

#Plot the interecept and the hyperparameter posteriors
plot(result, plot.random.effects=FALSE)

#Plot the posterior for the RANGE of the random field
range=inla.tmarginal(function(x) sqrt(8)/exp(x/2),(result$marginals.hyperpar$"K.0 for loc-basisK" ))
attr(range,"inla.tag")="Range"
plot(range)
