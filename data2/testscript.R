source("util.R")
library("RandomFields")

d <- 0.5
s1 <- seq(0+d/2, 20-d/2, by=d); s2 <- seq(0+d/2, 20-d/2, by=d)
z <- GaussRF(x=s1, y=s2, grid=T, model="stable",
param=c(0,1,0,2,1.5))
# image(s1, s2, z)

# meshgrid <- function(a,b) {
#   list(
#        x=outer(b*0,a,FUN="+"),
#        y=outer(b,a*0,FUN="+")
#        )
# } 
# xy <- meshgrid(s1,s2)
# x <- as.vector(xy[[1]])
# y <- as.vector(xy[[2]])

v <- 10
if(0) {
nx <- length(s1)
ny <- length(s2)
cov <- matrix(nrow=nx*ny,ncol=nx*ny)

for(i1 in 1:ny) {
	for(j1 in 1:nx) {
		p1k <- (i1-1)*ny+j1
		p1 <- c(s1[j1], s2[i1])
		for(i2 in 1:ny) {
			for(j2 in 1:nx) {
				p2k <- (i2-1)*ny+j2
				p2 <- c(s1[j2], s2[i2])
				#r <- (p1[1]-p2[1])^2+(p1[2]-p2[2])^2
				d <- p2-p1
				cov[p1k,p2k] <- exp(-sum(d*d)/v)
			}
		}
	}
	print(i1)
}
}
if(0) {
	cov.pd <- nearPD(cov)
}
if(0) {
	std <- chol(as.matrix(cov.pd$mat))
}

field1 <- rmvnorm(n=1,mean=rep(0,nx*ny),sigma=diag(nx*ny))
field <- matrix(std%*%t(field1),nrow=ny)
persp(s1,s2,field,col="lightblue",shade=0.25,expand=0.5,theta=30,phi=30)









