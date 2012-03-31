# Lansing wood data analysis
options(error=dump.frames)
exportFigs <- 1
displayFigs <- 0
interaction <- 1
speciesinteraction <- 1
intensity <- 1
ppcf <- 1
dimyx <- ifelse(exportFigs,c(500,500),c(100,100))
nsim <- 3

require("spatstat");
require("RColorBrewer")
data(lansing)
lansingm <- lansing

unitname(lansingm) <- list("metre","metres",1)
ft2m <- round(lansing$window$units$multiplier/3.2808399)
lansingm <- affine(lansingm,diag(c(ft2m,ft2m)))
range <- c(0,150)#lansingm$window$xrange

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

mar_lab <- c(2.5,2.5,1.5,1.0)
mar_tight <- c(0.1,0.1,0.1,0.1)
myplot <- function(...,width=6,height=6,mar=mar_lab,file=FALSE,nodevoff=FALSE,afterfn=NULL,k=NULL) {
	if(displayFigs) {
		quartz()
		par(mar=mar)
		p <- plot(...)
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
		p <- plot(...)
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

nlansing <- lansingm[lansingm$marks!="misc"];
sigma <- 2.5*bw.relrisk(nlansing);


levels(nlansing$marks) <- c("oak","hickory","maple",NA,"oak","oak")

hm <- lansingm[lansingm$marks=="maple" | lansingm$marks=="hickory"];
levels(hm$marks) <- c(NA,"hickory","maple",NA,NA,NA)

oaks <- lansing[grep("oak",lansing$marks)]
levels(oaks$marks) <- c("blackoak",NA,NA,NA,"redoak","whiteoak")

oakhm <- nlansing
levels(oakhm$marks) <- c("oak","hm","hm")


if(intensity) {
	
	rl=relrisk(nlansing,sigma=sigma,dimyx=dimyx)
	
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
print("INTERACTION")	

	snlansing <- split(nlansing)	

	Lss <- list(oak=snlansing$oak,hm=hm,all=nlansing)
	Ls <- mapply(envelope,Lss,list(Lest,Lest,Lest),
		MoreArgs=list(
			nsim=nsim,savefuns=TRUE,
			correction="Ripley",
			r=seq.int(range[1],range[2],(range[2]-range[1])/500)
		),SIMPLIFY=FALSE)
	
	dens <- density(split(nlansing),
		sigma=sigma)
	
	Lssi <- list(maple=snlansing$maple,hickory=snlansing$hickory)

	Lsi <- mapply(
		envelope,
		Lssi,
		list(Linhom,Linhom),
		simulate=list(expression(rpoispp(dens$maple)),expression(rpoispp(dens$hickory))),
		MoreArgs=list(
			nsim=nsim,savefuns=TRUE,
			correction="Ripley",
			normpower=2,
			sigma=sigma,
			r=seq.int(range[1],range[2],(range[2]-range[1])/500)
		),SIMPLIFY=FALSE)
	

	nms <- names(c(Lssi,Lss))


	mapply(listplot,nms,c(Lsi,Ls),
		MoreArgs=list(
			main="",
			formula=.-r~r,
			file="l_%s.pdf",
			legend=FALSE,
			width=3,
			height=3,
			mar=c(2.0,0.3,0.1,0.3),
			yaxt="n",
			xlim=c(0,150),
			lty=1,
			lwd=2
		))	
}

if(speciesinteraction) {
print("SPECIESINTERACTION")	
	legendfn <- function(p,k) {
				legend(
					'topright',
					c(k,"theoretical"),
					col=p$col[1:2],
					lty=1,
					lwd=2
				)
			}
	legendfn <- NULL

	# CSRI
	i <- c("hickory","hickory","maple")
	j <- c("oak","maple","oak")
	fns <- mapply(function(i,j){
			return(sprintf("%s_%s",i,j))
		},i,j,USE.NAMES=FALSE)
	
	Ls1 <- envelope(
			nlansing,
			Lcross,
			r=seq.int(range[1],range[2],(range[2]-range[1])/500),
			i=i[1],
			j=j[1],
			nsim=nsim,savefuns=TRUE,
			correction="Ripley",
			savepatterns=TRUE)	
			

	Ls <- mapply(
			envelope,
			rep(list(nlansing),2),
			rep(list(Lcross),2),
			i=i[2:3],
			j=j[2:3],
			MoreArgs=list(
				r=seq.int(range[1],range[2],(range[2]-range[1])/500),
				nsim=nsim,savefuns=TRUE,
				simulate=Ls1
			),SIMPLIFY=FALSE)	
	csrd <- c(list(Ls1),Ls)
	
	csrp <- mapply(listplot,fns,csrd,
		MoreArgs=list(
			lwd=2,
			lty=1,
			main="",
			formula=.-r~r,
			file="csri_%s.pdf",
			legend=FALSE,
			width=3,
			height=3,
			mar=c(2.0,0.3,0.1,0.3),
			yaxt="n",
			afterfn=legendfn
		),SIMPLIFY=FALSE)	
	
	# independence of components
print("IOC")
	i <- c("hickory","hickory","maple")
	j <- c("oak","maple","oak")
	fns <- mapply(function(i,j){
			return(sprintf("%s_%s",i,j))
		},i,j,USE.NAMES=FALSE)
	
	
	Ls1 <- envelope(
			nlansing,
			Lcross,
			i=i[1],
			j=j[1],
			r=seq.int(range[1],range[2],(range[2]-range[1])/500),
			nsim=nsim,savefuns=TRUE,
			correction="Ripley",
			simulate = expression(rshift(nlansing)),
			savepatterns=TRUE)	

	Ls <- mapply(
			envelope,
			rep(list(nlansing),2),
			rep(list(Lcross),2),
			i=i[2:3],
			j=j[2:3],
			MoreArgs=list(
				simulate = Ls1,
				r=seq.int(range[1],range[2],(range[2]-range[1])/500),
				nsim=nsim,savefuns=TRUE
			),SIMPLIFY=FALSE)	
	
	iocd <- c(list(Ls1),Ls)

	iocp <- mapply(listplot,fns,iocd,
		MoreArgs=list(
			lwd=2,
			lty=1,
			main="",
			formula=.-r~r,
			file="ioc_%s.pdf",
			legend=FALSE,
			width=3,
			height=3,
			mar=c(2.0,0.3,0.1,0.3),
			yaxt="n",
			afterfn=legendfn
		),SIMPLIFY=FALSE)	

	# random labeling
print("RANDOMLABELING")
	Ldif <- function(X, ..., i) { 
		Lidot <- Ldot(X, ..., i = i) 
		L <- Lest(X, ...)
		return(eval.fv(Lidot - L))
	}

	Ls1 <- envelope(
			nlansing,
			Ldif,
			i="hickory",
			r=seq.int(range[1],range[2],(range[2]-range[1])/500),
			nsim=nsim,savefuns=TRUE,
			correction="Ripley",
			simulate = expression(rlabel(nlansing)),
			savepatterns=TRUE)	

	Ls <- mapply(
			envelope,
			rep(list(nlansing),2),
			rep(list(Ldif),2),
			i=c("oak","maple"),
			MoreArgs=list(
				r=seq.int(range[1],range[2],(range[2]-range[1])/500),
				simulate = Ls1,
				nsim=nsim,savefuns=TRUE
			),SIMPLIFY=FALSE)
	
	rld <- c(list(Ls1),Ls)	
	
	rlp <- mapply(listplot,fns,rld,
		MoreArgs=list(
			main="",
			formula=.~r,
			file="rl_%s.pdf",
			legend=FALSE,
			width=3,
			height=3,
			mar=c(2.0,0.3,0.1,0.3),
			yaxt="n",
		 	lwd=2,
		 	lty=1,
			afterfn=legendfn
		),SIMPLIFY=FALSE)
	
}

if(ppcf) {
print("PPCF")
	dens <- density(split(nlansing))
	# ppcf inhomog
	bw <- 2.5*bw.stoyan(nlansing)
	i <- c("hickory","hickory","maple","hickory","maple","oak")
	j <- c("oak","maple","oak","hickory","maple","oak")
	fns <- mapply(function(i,j){
		return(sprintf("%s_%s",i,j))
	},i,j,USE.NAMES=FALSE)
	
	ppcfdi <- mapply(
			envelope,
			rep(list(nlansing),6),
			rep(list(pcfcross.inhom),6),
			i=i,
			j=j,
			lambdaI=list(dens[[i[1]]],dens[[i[2]]],dens[[i[3]]],dens[[i[6]]],dens[[i[4]]],dens[[i[5]]]),
			lambdaJ=list(dens[[j[1]]],dens[[j[2]]],dens[[j[3]]],dens[[j[6]]],dens[[j[4]]],dens[[j[5]]]),
			MoreArgs=list(
				#r=seq.int(range[1],range[2],(range[2]-range[1])/500),
				simulate=expression(
					rmpoispp(
						dens,
						types=names(dens)
					)
				),
				correction="Ripley",
				bw=bw,
				nsim=nsim,savefuns=TRUE
			),SIMPLIFY=FALSE)


	v <- mapply(listplot,fns[1:3],ppcfdi[1:3],
		MoreArgs=list(
			main="",
			formula=.~r,
			file="ppcfi_%s.pdf",
			legend=FALSE,
			width=3,
			height=3,
			mar=c(2.0,2.0,0.2,0.3),
		 	lwd=2,
		 	lty=1,
			ylim=c(0.3,1.3)
		))	
	v <- mapply(listplot,fns[4:6],ppcfdi[4:6],
		MoreArgs=list(
			main="",
			formula=.~r,
			file="ppcfi_%s.pdf",
			legend=FALSE,
			width=3,
			height=3,
			mar=c(2.0,2.0,0.2,0.3),
		 	lwd=2,
		 	lty=1,
		 	xlim=c(2,70.5)
		))
	

}

