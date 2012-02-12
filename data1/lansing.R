# Lansing wood data analysis

exportFigs <- 1
displayFigs <- 0
interaction <- 1
speciesinteraction <- 1
intensity <- 1
dimyx <- ifelse(exportFigs,c(500,500),c(100,100))

require("spatstat");
require("RColorBrewer")
#require("playwith");
data(lansing)
lansingm <- lansing

#unitname(lansingm) <- c("metre","metres",round(924/3.2808399))
unitname(lansingm) <- list("metre","metres",1)
ft2m <- round(lansing$window$units$multiplier/3.2808399)
lansingm <- affine(lansingm,diag(c(ft2m,ft2m)))
range <- lansingm$window$xrange

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

mar_lab <- c(2.5,2.5,1.5,1.0)
mar_tight <- c(0.1,0.1,0.1,0.1)
myplot <- function(...,width=6,height=6,mar=mar_lab,file=FALSE,nodevoff=FALSE) {
	if(displayFigs) {
		quartz()
		par(mar=mar)
		p <- plot(...)
		if(!exportFigs) {
			return(p)
		}
	}
	if(exportFigs && file != FALSE) {
		pdf(file=file,width=width,height=height)
		par(mar=mar)
		p <- plot(...)
		if(!nodevoff) {
			dev.off()
		}
		return(p)
	}
}
listplot <- function(k,v,file=FALSE,formula=FALSE,...) {
	if(formula != FALSE) {
		return(myplot(v,formula,file=sprintf(file,k),...))
	} else {
		return(myplot(v,file=sprintf(file,k),...))
	}
}

sigma <- bw.relrisk(lansing);


nlansing <- lansingm[lansingm$marks!="misc"];
levels(nlansing$marks) <- c("oak","hickory","maple",NA,"oak","oak")

hm <- lansingm[lansingm$marks=="maple" | lansingm$marks=="hickory"];
levels(hm$marks) <- c(NA,"hickory","maple",NA,NA,NA)

oaks <- lansing[grep("oak",lansing$marks)]
levels(oaks$marks) <- c("blackoak",NA,NA,NA,"redoak","whiteoak")

oakhm <- nlansing
levels(oakhm$marks) <- c("oak","hm","hm")

bw <- bw.diggle(nlansing)

if(intensity) {
	dens1 <- density(split(nlansing),
		sigma=2.5*sigma,
		dimyx=dimyx)
	
	rl=relrisk(nlansing,sigma=2.5*sigma)
	
	# smoothed intensities
	mapply(listplot,names(rl),rl,
		MoreArgs=list(
			zlim = c(0, 0.7),
			col=jet.colors(512),
			main="",
			mar=mar_tight,
			ribbon=FALSE,
			file="intensity_relative_%s.pdf",
			width=3,
			height=3
		))
	
	# original point patterns
	mapply(listplot,names(split(lansing)),split(lansing),
		MoreArgs=list(
			main="",
			mar=mar_tight,
			file="lansing_%s.pdf",
			width=3,
			height=3
		))
	
	# combined point patterns
	mapply(listplot,list("oaks","hm","all"),list(oaks,hm,lansing),
		MoreArgs=list(
			use.marks=FALSE,
			pch=21,
			main="",
			mar=mar_tight,
			file="lansing_%s_combined.pdf",
			width=3,
			height=3
		))
}

if(interaction) {
	
	patrns <- c(split(nlansing),list(hm=hm,all=nlansing))
	Ls <- mapply(envelope,patrns,list(Lest,Linhom,Linhom,Lest,Lest),
		MoreArgs=list(
			nsim=99,
			correction="Ripley",
			normpower=2,
			sigma=2.5*sigma
		),SIMPLIFY=FALSE)
	mapply(listplot,names(patrns),Ls,
		MoreArgs=list(
			main="",
			formula=.-r~r,
			file="l_%s.pdf",
			legend=FALSE,
			width=3,
			height=3,
			mar=c(2.0,0.3,0.1,0.3),
			yaxt="n"
		))
}

if(speciesinteraction) {
	

	# mark connection functions, pairwise
	bw <- 2*bw.stoyan(nlansing)
	i <- c("hickory","hickory","maple","hickory","maple","oak")
	j <- c("oak","maple","oak","hickory","maple","oak")
	markcs <- mapply(
			markconnect,
			rep(list(nlansing),length(i)),
			i,
			j,
			MoreArgs=list(
				r=seq.int(range[1],range[2],(range[2]-range[1])/500),
				correction="Ripley",
				bw=bw,
				normalise=FALSE
			),SIMPLIFY=FALSE)
	markc <- markcs[[1]]
	markc <- markc[,c("r","iso")]
	for(m in markcs[2:length(markcs)]) {
		markc <- cbind(markc,m[,c("r","iso")])
	}
	col <- sapply(brewer.pal(length(markc)-1,"Dark2"),function(c) {
			return(paste(c,as.hexmode(round(0.7*255)),sep=''))
		},USE.NAMES=FALSE)

	v <- myplot(
			markc,
			legend=FALSE,
			col=col,
			lty=1,
			lwd=4,
			ylab="mark-connection",
			main="",
			file="markc.pdf",nodevoff=TRUE,
			xlim=c(range[1],range[2]),
			ylim=c(0,0.3),
			width=5,
			height=5)
	legend('topright',
		mapply(function(i,j){
			return(sprintf("%s-%s",i,j))
		},rev(i),rev(j)),
		col=col,
		lwd=4,
		lty=v$lty)
	if(exportFigs) {
		dev.off()
	}

}

	
