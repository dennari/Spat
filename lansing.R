# Lansing wood data analysis

exportFigs <- 0
displayFigs <- 1

dimyx <- ifelse(exportFigs,c(500,500),c(100,100))

require("spatstat");
#require("playwith");
data(lansing)
lansingm <- lansing

#unitname(lansingm) <- c("metre","metres",round(924/3.2808399))
unitname(lansingm) <- list("metre","metres",round(924/3.2808399))

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

myplot <- function(...,file=FALSE) {
	if(displayFigs) {
		quartz()
		plot(...)
	}
	if(exportFigs && file) {
		pdf(file=file)
		plot(...)
		dev.off()	
	}
}

sigma <- bw.relrisk(lansing);
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



dens1 <- density(split(lansing),
	sigma=sigma,
	dimyx=dimyx)
myplot(dens1,
	file="lansing_intensity_separate.pdf",
	main="Independent intensity",
	sub="Isotropic Gaussian kernel",
	ribbon=FALSE,
	col=jet.colors(256))
	

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
	
