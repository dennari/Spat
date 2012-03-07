sp1 <- sapply(coordinates(map), function(x) do.call("rbind", x))
N <- length(sp1)
#sp1 <- sp1[c(1:14,17:N)]
N <- length(sp1)
lps <- matrix(nrow=N,ncol=2)
fps <- matrix(nrow=N,ncol=2)
NNN <- 0
for(i in 1:N) {
	line <- sp1[[i]]
	NN <- nrow(line)
	NNN <- NNN+NN
	lps[i,] <- line[NN,]	
	fps[i,] <- line[1,]	
	
	# for(j in 1:N) {
	#  	mi <- colMeans(sp1[[i]])
	#  	mj <- colMeans(sp1[[j]])
	#  	dist[i,j] = sqrt((mi[1]-mj[1])^2+(mi[1]-mj[1])^2)
	# }
}
d <- as.matrix(dist(mns))
d[d==0]=NA
c <- vector("integer",N)
c[1] <- 1
coords <- sp1[[c[1]]][nrow(sp1[[c[1]]]):1,]

for(i in 1:(N-1)) {
	dmin <- 99999999
	jmin <- 0
	fp <- fps[c[i],]
	for(j in 1:N) {
		lp <- lps[j,]
		d <- dist(rbind(fp,lp))[1]

		if(!is.element(j,c) & d < dmin) {
			jmin <- j
			dmin <- d 
		}
	}
	if(jmin == 16)
		break
	# print(c[i])
	# print(jmin)
	# print(dmin)
	c[i+1]<-jmin
	nline <- sp1[[jmin]]
	coords <- rbind(coords,nline[nrow(nline):1,])
	#coords <- rbind(coords,nline)
}
coords <- rbind(coords,getAllCoords(map,c(16,3,4,5)))

# coords <- matrix(nrow=0,ncol=2)
# #c <- c(1,unique(c))
# for(i in 1:3) {
# 	line <- sp1[c(i)]
# 	NN <- nrow(line)
# 	coords <- rbind(coords,line[NN:1,])
# }

# coords2 <- matrix(nrow=0,ncol=2)
# for(i in 1:N) {
# 	coords2 <- rbind(coords2,fps[i,])
# 	coords2 <- rbind(coords2,lps[c(i),])
# }