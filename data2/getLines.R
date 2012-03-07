lines <- list()
for(Lines in map@lines) {
	for(line in Lines@Lines) {
		lines <- c(lines,line)
	} 
}