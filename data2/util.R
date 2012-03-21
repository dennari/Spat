require("maptools")
require("spdep")
require("RColorBrewer")
require("spatstat")
require("INLA")
require("splancs")

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  jet.colorf <-
  colorRamp(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

addAlpha <- function(c,op) {
	return(paste(c,as.hexmode(as.integer(op*255)),sep=""))
}

mar_lab <- c(2.5,2.5,1.5,1.0)
mar_tight <- c(0.1,0.1,0.1,0.1)
myplot <- function(...,plotfn=plot,width=6,height=6,mar=mar_lab,file=FALSE,nodevoff=FALSE,afterfn=NULL,k=NULL) {
	if(displayFigs) {
		quartz()
		par(mar=mar)
		p <- plotfn(...)
		if(is.function(afterfn)) {
			afterfn(p,k)
		}
		if(!exportFigs) {
			return(p)
		}
	}
	if(exportFigs && file != FALSE) {
		pdf(file=file,width=width,height=height)
		par(mar=mar)
		p <- plotfn(...)
		if(is.function(afterfn)) {
			afterfn(p,k)
		}
		if(!nodevoff) {
			dev.off()
		}
		return(p)
	}
}
listplot <- function(k,v,file=FALSE,formula=FALSE,main="",...) {
	if(formula != FALSE) {
		p <- myplot(v,formula,file=sprintf(file,k),main=sprintf(main,k),k=k,...)
	} else {
		p <- myplot(v,file=sprintf(file,k),main=sprintf(main,k),k=k,...)
	}
	return(p)
}


getLines <- function(map) {
	lines <- list()
	for(Lines in map@lines) {
		for(line in Lines@Lines) {
			lines <- c(lines,line)
		} 
	}
}

getAllCoords <- function(sp,ic) {
	
	sp1 <- sapply(coordinates(sp), function(x) do.call("rbind", x))
	N <- length(ic)
	nline <- sp1[[ic[1]]]
	coords <- nline[1:nrow(nline),]
	for(i in 2:N) {
		#print(ic[i])
		nline <- sp1[[ic[i]]]
		coords <- rbind(coords,nline[1:nrow(nline),])
	}
	return(coords)
}

