# Lansing wood data analysis

exportFigs <- 1
displayFigs <- 0
interaction <- 0
speciesinteraction <- 0
intensity <- 1
dimyx <- ifelse(exportFigs,c(500,500),c(100,100))

require("spatstat");
require("RColorBrewer")
#require("playwith");
data(lansing)
lansingm <- lansing

#unitname(lansingm) <- c("metre","metres",round(924/3.2808399))
unitname(lansingm) <- list("metre","metres",round(924/3.2808399))

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
		return(p)
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

sigma <- bw.relrisk(lansing);
# # plot everything


# plot oaks combined
# oaks = lansing[
# 	lansing$marks=="whiteoak" | 
# 	lansing$marks=="blackoak" | 
# 	lansing$marks=="redoak"]
# myplot(oaks,
# 	use.marks=FALSE,
# 	pch=21,
# 	main="Oaks",
# 	file="lansing_oaks.pdf")
# myplot(density(oaks,sigma=2.5*sigma,dimyx=dimyx),
# 	file="lansing_intensity_oaks.pdf",
# 	main="Intensity - oaks",
# 	ribbon=FALSE,
# 	col=jet.colors(256))
# myplot(density(lansing,sigma=2.5*sigma,dimyx=dimyx),
# 	file="lansing_intensity_combined.pdf",
# 	main="Intensity - trees",
# 	ribbon=FALSE,
# 	col=jet.colors(256))

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
	
	myplot(dens1,
		file="intensity_separate.pdf",
		main="",
		sub="Isotropic Gaussian kernel",
		ribbon=FALSE,
		mar=mar_tight,
		col=jet.colors(256))	

	rl=relrisk(nlansing,sigma=2.5*sigma)
	myplot(rl$oak,
		zlim = c(0, 0.7),
		col=jet.colors(512),
		main="",
		mar=mar_tight,
		ribbon=FALSE,
		file="intensity_relative_oak.pdf",
		width=3,
		height=3
		)
	myplot(rl$hickory,
		zlim = c(0, 0.7),
		col=jet.colors(512),
		main="",
		mar=mar_tight,
		ribbon=FALSE,
		file="intensity_relative_hickory.pdf",
		width=3,
		height=3
		)
	myplot(rl$maple,
		zlim = c(0, 0.7),
		col=jet.colors(512),
		main="",
		mar=mar_tight,
		ribbon=FALSE,
		file="intensity_relative_maple.pdf",
		width=3,
		height=3
		)

	myplot(split(nlansing)$hickory,
		main="",
		mar=mar_tight,
		file="lansing_hickory.pdf",
		width=3,
		height=3)	
	myplot(split(nlansing)$maple,
		main="",
		mar=mar_tight,
		file="lansing_maple.pdf",
		width=3,
		height=3)
	myplot(split(nlansing)$oak,
		main="",
		mar=mar_tight,
		file="lansing_oak.pdf",
		width=3,
		height=3)

	myplot(split(lansing),
		main="",
		mar=mar_tight,
		file="lansing.pdf")

	myplot(oaks,
		use.marks=FALSE,
		pch=21,
		main="",
		mar=mar_tight,
		file="lansing_oaks_combined.pdf",
		width=3,
		height=3)
	
	myplot(hm,
			use.marks=FALSE,
			pch=21,
			main="",
			mar=mar_tight,
			file="lansing_hm_combined.pdf",
			width=3,
			height=3)
	
	myplot(lansing,
			use.marks=FALSE,
			pch=21,
			main="",
			mar=mar_tight,
			file="lansing_combined.pdf",
			width=3,
			height=3)


}

if(interaction) {
	myplot(envelope(split(nlansing)$oak,Lest,nsim=20,correction="best",normpower=2,sigma=2.5*sigma),
		.-r~r,
		main="",
		file="l_oak.pdf",
		legend=FALSE,
		width=3,
		height=3)
	myplot(envelope(split(nlansing)$hickory,Linhom,nsim=20,correction="best",normpower=2,sigma=2.5*sigma),
		.-r~r,
		main="",
		file="l_hickory.pdf",
		legend=FALSE,
		width=3,
		height=3)
	myplot(envelope(split(nlansing)$maple,Linhom,nsim=20,correction="best",normpower=2,sigma=2.5*sigma),
		.-r~r,
		main="",
		file="l_maple.pdf",
		legend=FALSE,
		width=3,
		height=3)	
	#myplot(envelope(hm,Lest,nsim=20,correction="best",normpower=2,sigma=2.5*sigma),.-r~r)
}

if(speciesinteraction) {
	col <- brewer.pal(7,"Dark2")
	bw <- 2*bw.stoyan(nlansing)
	markchm <- markconnect(nlansing,"hickory","maple",correction="Ripley",bw=bw,normalise=TRUE)
	markcho <- markconnect(nlansing,"hickory","oak",correction="Ripley",bw=bw,normalise=TRUE)
	markchh <- markconnect(nlansing,"hickory","hickory",correction="Ripley",bw=bw,normalise=TRUE)
	markcmm <- markconnect(nlansing,"maple","maple",correction="Ripley",bw=bw,normalise=TRUE)
	markcmo <- markconnect(nlansing,"maple","oak",correction="Ripley",bw=bw,normalise=TRUE)
	markcoo <- markconnect(nlansing,"oak","oak",correction="Ripley",bw=bw,normalise=TRUE)
	markc <- cbind(
			markchm,
			markcho[,cbind("r","iso")],
			markchh[,cbind("r","iso")],
			markcmm[,cbind("r","iso")],
			markcmo[,cbind("r","iso")],
			markcoo[,cbind("r","iso")])

	v <- myplot(
			markc,
			legend=FALSE,
			col=col,
			lty=1,
			lwd=2,
			ylab="mark-connection",
			main="",
			file="markc.pdf",nodevoff=TRUE,ylim=c(0.3,2.0),
			width=5,
			height=5)
	legend('topright',
		c('oo','mo','mm','hh','ho','hm','theo'),
		col=col,
		lwd=2,
		lty=v$lty)
	if(exportFigs) {
		dev.off()
	}

	pcfhm <- pcfcross(nlansing,"hickory","maple",correction="Ripley",bw=bw)
	pcfho <- pcfcross(nlansing,"hickory","oak",correction="Ripley",bw=bw)
	pcfhh <- pcfcross(nlansing,"hickory","hickory",correction="Ripley",bw=bw)
	pcfmm <- pcfcross(nlansing,"maple","maple",correction="Ripley",bw=bw)
	pcfmo <- pcfcross(nlansing,"maple","oak",correction="Ripley",bw=bw)
	pcfoo <- pcfcross(nlansing,"oak","oak",correction="Ripley",bw=bw)
	pcf <- cbind(
			pcfhm,
			pcfho[,cbind("r","iso")],
			pcfhh[,cbind("r","iso")],
			pcfmm[,cbind("r","iso")],
			pcfmo[,cbind("r","iso")],
			pcfoo[,cbind("r","iso")])
	vv <- myplot(
			pcf,
			legend=FALSE,
			col=col,
			lty=1,
			lwd=2,
			ylab="pcf",
			main="",
			ylim=c(0.4,1.8),
			xlim=c(0.01,0.23),
			file="pcf.pdf",nodevoff=TRUE,width=5,
			height=5)

	legend('topright',
		c('oo','mo','mm','hh','ho','hm','theo'),
		col=col,
		lwd=2,
		lty=vv$lty)

	if(exportFigs) {
		dev.off()
	}

	pcfhm <- pcfcross.inhom(nlansing,"hickory","maple",lambdaI=dens1$hickory,lambdaJ=dens1$maple,correction="Ripley",stoyan=0.3,bw=bw)
	pcfho <- pcfcross.inhom(nlansing,"hickory","oak",lambdaI=dens1$hickory,lambdaJ=dens1$oak,correction="Ripley",stoyan=0.3,bw=bw)
	pcfhh <- pcfcross.inhom(nlansing,"hickory","hickory",lambdaI=dens1$hickory,lambdaJ=dens1$hickory,correction="Ripley",stoyan=0.3,bw=bw)
	pcfmm <- pcfcross.inhom(nlansing,"maple","maple",lambdaI=dens1$maple,lambdaJ=dens1$maple,correction="Ripley",stoyan=0.3,bw=bw)
	pcfmo <- pcfcross.inhom(nlansing,"maple","oak",lambdaI=dens1$maple,lambdaJ=dens1$oak,correction="Ripley",stoyan=0.3,bw=bw)
	pcfoo <- pcfcross.inhom(nlansing,"oak","oak",lambdaI=dens1$oak,lambdaJ=dens1$oak,correction="Ripley",stoyan=0.3,bw=bw)
	pcf2 <- cbind(
			pcfhm,
			pcfho[,cbind("r","iso")],
			pcfhh[,cbind("r","iso")],
			pcfmm[,cbind("r","iso")],
			pcfmo[,cbind("r","iso")],
			pcfoo[,cbind("r","iso")])

	vv <- myplot(
			pcf2,
			.~r,
			legend=FALSE,
			col=col,
			lty=1,
			lwd=2,
			ylab="pcf",
			main="",
			ylim=c(0.5,1.4),
			xlim=c(0.01,0.23),
			file="pcf_inhom.pdf",nodevoff=TRUE,	width=5,
			height=5)
	legend('topright',
		c('oo','mo','mm','hh','ho','hm','theo'),
		col=col,
		lwd=2,
		lty=vv$lty)
	if(exportFigs) {
		dev.off()
	}
}
# maple vs hickory

# myplot(nlansing[nlansing$marks!="oak"],
# 	use.marks=TRUE,
# 	chars=c(19,19,19),
# 	cols=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
# 	main="Hickory vs. Maple",
# 	file="lansing_oaks.pdf")
# myplot(nlansing,
# 	use.marks=TRUE,
# 	chars=c(19,19,19),
# 	cols=c(rgb(0,0.8,0,0.7),rgb(0,0,0.8,0.7),rgb(0.8,0,0,0.7)),
# 	main="All",
# 	file="lansing_oaks.pdf")
# myplot(oaks,
# 	use.marks=FALSE,
# 	chars=c(19),
# 	cols=rgb(0,0.6,0,0.5),
# 	main="Oak",
# 	file="lansing_oaks.pdf")
# myplot(oaks,
# 	use.marks=TRUE,
# 	chars=c(19,19,19),
# 	cols=c(rgb(0,0.2,0,0.5),rgb(0,0.6,0,0.5),rgb(0,1.0,0,0.5)),
# 	main="Different oaks",
# 	file="lansing_oaks.pdf")
# rl=relrisk(nlansing,sigma=2.5*sigma)
# myplot(rl,zlim = c(0, 0.7),col=jet.colors(512))
# rl2=relrisk(oaks,sigma=2.5*sigma)
# myplot(rl2,zlim = c(0, 0.7),col=jet.colors(512))
# rl3=relrisk(oakhm,sigma=2.5*sigma,casecontrol=FALSE)
# myplot(rl3,zlim = c(0, 0.7),col=jet.colors(512))



# m <- split(lansing)$maple
# myplot(envelope(m,Lest,nsim=19,correction="Ripley",global=TRUE));
# myplot(envelope(m,Linhom,nsim=19,correction="Ripley",global=TRUE));
# myplot(envelope(m,Linhom,nsim=19,correction="Ripley",global=TRUE));

#L <- alltypes(lansing,Lcross.inhom,envelope=TRUE,sigma=sigma)

# crosses <- list()
# levels <- levels(lansing$marks)
# for(l in 1:(length(levels)-1)) {
# 	for(k in (l+1):length(levels)) {
# 		t <- list()
# 		t$pcf <- pcfcross.inhom(lansingm,levels[l],levels[k],correction="Ripley")
# 		t$L <- Lcross.inhom(lansingm,levels[l],levels[k],correction="Ripley")
# 		t$main <- sprintf("%s - %s",levels[l],levels[k])
# 		crosses[[length(crosses)+1]] <- t
# 	}
# }

# quartz(width=14,height=12)
# par(mfrow=c(5,2))
# for(c in crosses[6:10]) {
# 	plot(c$pcf,main=c$main,xlim=c(0.01,0.24),ylab="",xlab="",legend=FALSE)
# 	plot(c$L,cbind(iso-theo,theo-theo)~r,main=c$main,ylab="",xlab="",legend=FALSE)
# }




	

# dens2 <- relrisk(lansing,
# 	sigma=sigma,
# 	dimyx=dimyx)
# myplot(dens2, 
# 	file="lansing_intensity_relative.pdf",
# 	zlim = c(0, 1),
# 	main="Relative intensity",
# 	ribbon=FALSE,
# 	col=jet.colors(256)
# 	)
	
