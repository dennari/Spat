
listplot <- function(k,v,file=FALSE,...) {
	return(myplot(v,file=sprintf(file,k),...))
}

a <- mapply(listplot,names(rl),rl,
		MoreArgs=list(
		zlim = c(0, 0.7),
		main="",
		mar=mar_tight,
		ribbon=FALSE,
		file="intensity_relative_%s.pdf",
		width=3,
		height=3
		))
print(a)

	for(m in markcs[2:3]) {
		print(m)
	}