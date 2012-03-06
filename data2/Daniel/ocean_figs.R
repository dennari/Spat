###########################################################################
## ocean_figs.R
## Code for plotting the figures from Section 8.3 of Simpson et al. (2011)
## Author: Daniel Simpson
## This code comes with no warranty or guarantee of any kind.
###########################################################################

graphics.off()
rgl.quit()
require(rgl)
COLOUR = TRUE


windowRect.globe = c(50,50,50+840,50+400)

if (COLOUR) { 
  cp = colorRampPalette(c("darkblue", "blue", "cyan", "yellow", "red", "darkred"))
  my.palette = cp
} else
{
  ## Construct greyscale palette function:
  my.grey.palette = function (n,...) { return (grey.colors(n,0.05,0.95,...))}
  ## Use it:
  my.palette = my.grey.palette
}
cp = my.palette

draw.globe = function(mesh, sample, draw.edges=FALSE)
{
    mesh0 = inla.mesh.create(loc=cbind(0,0), extend=list(offset=1.1,n=4))

    mesh01 = old.mesh.class(mesh0)
    mesh02 = old.mesh.class(mesh0)
    mesh1 = old.mesh.class(mesh)
    mesh2 = old.mesh.class(mesh)
    mesh02$mesh$s[,1] = mesh02$mesh$s[,1]*(-1)
    mesh02$mesh$s[,3] = mesh02$mesh$s[,3]*(-1)
    mesh2$mesh$s[,1] = mesh2$mesh$s[,1]*(-1)
    mesh2$mesh$s[,3] = mesh2$mesh$s[,3]*(-1)

    mesh01$mesh$s[,1] = mesh01$mesh$s[,1]-1.1
    mesh02$mesh$s[,1] = mesh02$mesh$s[,1]+1.1
    mesh1$mesh$s[,1] = mesh1$mesh$s[,1]-1.1
    mesh2$mesh$s[,1] = mesh2$mesh$s[,1]+1.1

    open3d(windowRect=windowRect.globe)
    plot(mesh01, col="white", color.palette=cp,
         draw.vertices=FALSE, draw.edges=FALSE, add=TRUE)
    plot(mesh02, col="white", color.palette=cp,
         draw.vertices=FALSE, draw.edges=FALSE, add=TRUE)
    plot(mesh1, col=sample, color.palette=cp,
         draw.vertices=FALSE, draw.edges=draw.edges, add=TRUE)
    plot(mesh2, col=sample, color.palette=cp,
         draw.vertices=FALSE, draw.edges=draw.edges, add=TRUE)

    view3d(0,0,fov=0,zoom=0.4)
    rgl.bringtotop()
}


finn.do.print=FALSE


draw.globe(mesh,result$summary.random$idx[,"mean"]*NA, draw.edges=TRUE)
if (finn.do.print) rgl.snapshot("/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_mesh.png", top=TRUE)

##draw.globe(mesh,result$summary.random$idx[,"sd"])

source("utils.R")

proj = inla.mesh.projector(mesh, dims=c(360,180), projection='longsinlat')
proja = proj
proja$x = proja$x*pi/180



if (FALSE) {
## Note: the y-axis labels on this plot are wrong!!!
## Need my.levelplot-like method to fix it.
dev.new()
ptsll = cbind(atan2(pts[,2],pts[,1])*180/pi,pts[,3]*90)
plot(ptsll[,1],ptsll[,2],xlab='Longitude',ylab='Latitude',pch=".")
if (finn.do.print)
    dev.copy2pdf(width=4.5,height=4,file="/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_pts.pdf")
}

#make sure the axes are equal
at.mean.and.sample =  pretty(c(min(c(result$summary.linear.predictor[nV+nData+(1:nV),"mean"],intercept + sample),na.rm=TRUE), max(c(result$summary.linear.predictor[nV+nData+(1:nV),"mean"],intercept + sample),na.rm=TRUE)),11)


## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj, result$summary.random$idx[,"mean"])
##ptsll = cbind(atan2(pts[,2],pts[,1])*180/pi,pts[,3]*90)
ptsll = cbind(atan2(pts[,2],pts[,1]),pts[,3]*2)
## Plot PM contours:
dev.new()
my.levelplot(proja, plotdata*NA, map=globemapflat, color.palette=cp,
             at=at.mean.and.sample,
             contour=TRUE,
             colorkey=FALSE,
             main='')
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(ptsll[,1], ptsll[,2],col=1,cex=2,pch=46)
trellis.unfocus()
if (finn.do.print)
    dev.copy2pdf(width=4.5,height=4,file="/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_points.pdf")


att.range=range(c(range(result$summary.linear.predictor[nV+nData+(1:nV),"0.025quant"]),
range(result$summary.linear.predictor[nV+nData+(1:nV),"0.975quant"])))
att=pretty(att.range,21)

## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj, result$summary.linear.predictor[nV+nData+(1:nV),"mean"])
##ptsll = cbind(atan2(pts[,2],pts[,1])*180/pi,pts[,3]*90)
ptsll = cbind(atan2(pts[,2],pts[,1]),pts[,3]*2)
## Plot PM contours:
dev.new()
my.levelplot(proja, plotdata, map=globemapflat, color.palette=cp,
             at=at.mean.and.sample,
             contour=TRUE,
             colorkey=list(space='top'),
             main='')
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(ptsll[,1], ptsll[,2],col=1,cex=2,pch=46)
trellis.unfocus()
if (finn.do.print)
    dev.copy2pdf(width=4.5,height=4,file="/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_mean.pdf")


if (FALSE) {
## Map resulting posterior IQR field to a grid:
plotdata = inla.mesh.project(proj, result$summary.linear.predictor[nV+nData+(1:nV),"0.975quant"]-result$summary.linear.predictor[nV+nData+(1:nV),"0.025quant"])
##ptsll = cbind(atan2(pts[,2],pts[,1])*180/pi,pts[,3]*90)
ptsll = cbind(atan2(pts[,2],pts[,1]),pts[,3]*2)
## Plot PM contours:
dev.new()
my.levelplot(proja, plotdata, map=globemapflat, color.palette=cp,
             at=att,
             contour=TRUE,
             colorkey=list(space='top'),
             main='')
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(ptsll[,1], ptsll[,2],col=1,cex=2,pch=46)
trellis.unfocus()
if (finn.do.print)
    dev.copy2pdf(width=4.5,height=4,file="/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_iqr.pdf")
}



if (FALSE){
  #Plot risk map using Gaussian Approximation - very fast!)
u.level=5.5
## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj,
   pnorm((result$summary.linear.predictor[nV+nData+(1:nV),"mean"]-u.level)/
         result$summary.linear.predictor[nV+nData+(1:nV),"sd"]))
##ptsll = cbind(atan2(pts[,2],pts[,1])*180/pi,pts[,3]*90)
ptsll = cbind(atan2(pts[,2],pts[,1]),pts[,3]*2)
## Plot PM contours:
dev.new()
my.levelplot(proja, plotdata, map=globemapflat, color.palette=cp,
             at=pretty(c(0,1),11),
             contour=TRUE,
             colorkey=list(space='top'),
             main='')
##trellis.focus("panel", 1, 1, highlight=FALSE)
##lpoints(ptsll[,1], ptsll[,2],col=1,cex=2,pch=46)
##trellis.unfocus()
if (finn.do.print)
    dev.copy2pdf(width=4.5,height=4,file="/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_prob.pdf")

}





if (TRUE) {

## Map resulting posterior sd to a grid:
## Overestimation of the std.dev., since the intercept and spatial effect are
## negatively correlated in the posterior.  Most prominent for large ranges.
plotdata = inla.mesh.project(proj, (result$summary.linear.predictor[nV+nData+(1:nV),"sd"]))
## Plot PM contours:
dev.new()
my.levelplot(proja, plotdata, map=globemapflat, color.palette=cp,
             at=pretty(c(min(plotdata,na.rm=TRUE),
                         max(plotdata,na.rm=TRUE)),11),
             contour=TRUE,
             colorkey=list(space='top'),
             main='')
if (finn.do.print)
    dev.copy2pdf(width=4.5,height=4,file="/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_sd.pdf")



}



## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj, intercept +  sample)
##ptsll = cbind(atan2(pts[,2],pts[,1])*180/pi,pts[,3]*90)
ptsll = cbind(atan2(pts[,2],pts[,1]),pts[,3]*2)
## Plot PM contours:
dev.new()
my.levelplot(proja, plotdata, map=globemapflat, color.palette=cp,
             at=at.mean.and.sample,
             contour=TRUE,
              colorkey=list(space='top'),
             main='')
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(ptsll[,1], ptsll[,2],col=1,cex=2,pch=46)
trellis.unfocus()
if (finn.do.print)
    dev.copy2pdf(width=4.5,height=4,file="/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_sample.pdf")


#Get the proper risk map using inla.pmarginal (slower)
u.level=5.5
risk.level= rep(-100,nV)
for (i in 1:nV) {
  N = as.character( nV + nData + i)
  risk.level[i] = inla.pmarginal(x=u.level,marginal=result$marginals.linear.predictor[[paste("Apredictor.",N,sep="")]])
}


## Map resulting posterior mean field to a grid:
plotdata = inla.mesh.project(proj,1-risk.level)
##ptsll = cbind(atan2(pts[,2],pts[,1])*180/pi,pts[,3]*90)
ptsll = cbind(atan2(pts[,2],pts[,1]),pts[,3]*2)
## Plot PM contours:
dev.new()
my.levelplot(proja, plotdata, map=globemapflat, color.palette=cp,
             at=pretty(c(0,1),11),
             contour=TRUE,
             colorkey=list(space='top'),
             main='')
##trellis.focus("panel", 1, 1, highlight=FALSE)
##lpoints(ptsll[,1], ptsll[,2],col=1,cex=2,pch=46)
##trellis.unfocus()
if (finn.do.print)
    dev.copy2pdf(width=4.5,height=4,file="/Users/danielsimpson/Dropbox/Janine\ and\ Dan/Point\ Processs\ with\ SPDEs/ocean_prob_correct.pdf")
