###########################################################################
## ocean_mesh.R
## Constructs the mesh over the oceans Section 8.3 of Simpson et al. (2011)
## Author: Daniel Simpson
## This code comes with no warranty or guarantee of any kind.
###########################################################################

library(INLA)

true.radius.of.earth = 6371
radius.of.earth = 1

aus <- read.table("australia.txt")
aus = rbind(aus,aus[1,])
ant <- read.table("antarctica.txt")
ant <- rbind(ant,ant[1,])
#ams <-  read.table("americas.txt")
#ams <- rbind(ams,ams[1,])
eur <-  read.table("eurasia_africa_americas.txt")
eur <- rbind(eur,eur[1,])

lonlat3D=function(lon,lat){
  cbind(cos((lon/180)*pi)*cos((lat/180)*pi),
          sin((lon/180)*pi)*cos((lat/180)*pi),
          sin((lat/180)*pi))
}

#Make everything 3D
aus3 <- lonlat3D(aus[,1],aus[,2])
ant3 <- lonlat3D(ant[,1],ant[,2])
#ams3 <- lonlat3D(ams[,1],ams[,2])
eur3 <- lonlat3D(eur[,1],eur[,2])

segm_aus = inla.mesh.segment(loc=aus3)
segm_ant = inla.mesh.segment(loc=ant3)
#segm_ams = inla.mesh.segment(loc=ams3)
segm_eur = inla.mesh.segment(loc=eur3)

segm_world = list(segm_aus,segm_ant,segm_eur)

#mesh1=inla.mesh.create(boundary=segm_eur,cutoff=10/true.radius.of.earth,refine=list(max.edge=10/true.radius.of.earth))

mesh=inla.mesh.create(boundary=segm_world,cutoff=100/true.radius.of.earth,refine=list(max.edge=500/true.radius.of.earth))




