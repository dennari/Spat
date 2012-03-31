envs.typeIerror = function(x, k=1, rmin=min(attributes(x)$simfuns$r), rmax=max(attributes(x)$simfuns$r)) {
    # (c) Mari Myllymäki & Pavel Grabarnik
    # Ref. Grabarnik, P., Myllymäki, M. and Stoyan, D. (2011). 
    # Correct testing of mark independence for marked point patterns. Ecological Modelling 222, 3888-3894. 
    # Calculates an approximation of the type I error for the envelope test
	# x = an object returned by envelopes using argument savefuns=TRUE, i.e. including the functions from simulations
    # k = The kth lowest and highest simulation will form the envelopes.
    # rmin = Minimum value of r to be used for testing mark hypothesis by deviation and envelope tests. Default is the minimum value of r for which envelopes have been calculated.
    # rmax = Maximum value of r to be used for testing mark hypothesis by deviation and envelope tests. Default is the maximum value of r for which envelopes have been calculated.

	if(is.null(attributes(x)$simfuns)) stop("Error: Functions are not saved in x!\n")
		
	# r, chosen indices
	if(rmin > rmax) stop("Cannot have rmin > rmax!\n")
	#if(rmin < min(x$r)) cat("Notice: your rmin is smaller than x$r, only those r-values >= min(x$r) will be used. \n")
	#if(rmax > max(x$r)) cat("Notice: your rmax is larger than x$r, only those r-values <= max(x$r) will be used. \n")
	rranks = which(rmin<=x$r & x$r <=rmax)
	# r, chosen indices in descreasing order
	rranksinv = sort(rranks, decreasing=TRUE)
		
	# Take the simulations from "x" into a matrix so that it is easier to handle with them
	NofSimulations = ncol(attributes(x)$simfuns)-1
	sims = matrix(ncol=NofSimulations, nrow=nrow(attributes(x)$simfuns))
	for(i in 1:NofSimulations) {
		sims[,i] = attributes(x)$simfuns[[i+1]]
	}	

	if(k < 1) stop("Unreasonable value of k. \n")
	RANKLO = k
	RANKUP = (NofSimulations+1) - k
	
	# initialization of vectors
	sims.type = rep("normal", times=NofSimulations)	  # "normal" (inside envelopes), "strange" (contributing or outside), "false strange" (see below sims.type.false)
	sims.NofStrangeR = rep(0, times=NofSimulations)   # "strange" (contributing or outside)
	sims.type.false = rep(0, times=NofSimulations)    # a simulation that is contributing together with other simulations for some distance, but it not taken as contributing simulation for this distance, and is not strange or contributing for any other distance  
	
	# ranks for each simulation 
	sim.ranksL = rank(sims[rranksinv[1],], ties.method="max")
	sim.ranksU = rank(sims[rranksinv[1],], ties.method="min")
	# ids of those simulations that give the kth lowest and highest value
	rankL = which( sim.ranksL == sort(sim.ranksL)[RANKLO] )  
	rankU = which( sim.ranksU == sort(sim.ranksU)[RANKUP] )
	# first handle the lower envelope
	if(length(rankL)!=1) { # if there are ties 
		apu3 = sample(1:length(rankL), size=1) # sample one of the ties
		sims.type[rankL[apu3]] = "contributing" # characterize the chosen simulation as "contributing"
		sims.NofStrangeR[rankL[apu3]] = sims.NofStrangeR[rankL[apu3]] + 1 # increase the number "strange", i.e. contributing or outside, for this simulation
		sims.type.false[rankL[-apu3]] = 1  # other simulations (for which the tie occurred) except apu3 classified as "false strange"
		rankL = rankL[apu3] # the chosen
	}
	else { # no ties: characterize the simulation as "contributing" and increase the number "strange", i.e. contributing or outside, for this simulation
		sims.type[rankL] = "contributing"
		sims.NofStrangeR[rankL] = sims.NofStrangeR[rankL] + 1
	}
	# the same for the upper envelope 
	if(length(rankU)!=1) {
		apu3 = sample(1:length(rankU), size=1) 
		sims.type[rankU[apu3]] = "contributing"
		sims.NofStrangeR[rankU[apu3]] = sims.NofStrangeR[rankU[apu3]] + 1
		sims.type.false[rankU[-apu3]] = 1 
		rankU = rankU[apu3]
	}	
	else {
		sims.type[rankU] = "contributing"
		sims.NofStrangeR[rankU] = sims.NofStrangeR[rankU] + 1
	}
	# values for the envelopes
	envelopeL = envelopeU = vector(length=length(x$r)) # extract at the end those values that are chosen by [rmin, rmax]
	envelopeL[rranksinv[1]] = sims[rranksinv[1],rankL[1]] # [1] since all envelopeL.sims give the same value (if length(rankL) > 1)
	envelopeU[rranksinv[1]] = sims[rranksinv[1],rankU[1]] # [1] since all envelopeL.sims give the same value (if length(rankU) > 1)
	# ids of those simulations that are strictly outside envelopes
    outside = which( sim.ranksL < RANKLO | sim.ranksU > RANKUP )
	if(length(outside)>0) {
		# characterize these simulations as "outside"
		sims.type[outside] = "outside"
		# and increase the number "strange", i.e. contributing or outside, for this simulation
		sims.NofStrangeR[outside] = sims.NofStrangeR[outside] + 1
	}
	# then go through other values of r
	for(i in rranksinv[-1]) { 
		sim.ranksL = rank(sims[i,], ties.method="max")
		sim.ranksU = rank(sims[i,], ties.method="min")
		rankL = which( sim.ranksL == sort(sim.ranksL)[RANKLO] )
		rankU = which( sim.ranksU == sort(sim.ranksU)[RANKUP] )
		if(length(rankL)!=1) { # if there are ties
			# check first which of the simulations is the most "strange", i.e. has most small or large *strange* (*outside* or *contributing*) ranks for all r[k]
			apu = which( sims.NofStrangeR[rankL] == max(sims.NofStrangeR[rankL]) )
			if(length(apu)==1) {
				if(sims.type[rankL[apu]]!="outside") sims.type[rankL[apu]] = "contributing" # if the simulation has not already been classified as "outside" (or contributing), give it the status "contributing", i.e. do not make change "outside" -> "contributing"
				sims.NofStrangeR[rankL[apu]] = sims.NofStrangeR[rankL[apu]] + 1             # increase the number "strange", i.e. contributing or outside, for this simulation
				# if the other simulations are "normal", i.e. inside envelopes for other r[k], classify them as "false strange"
				for(j in rankL[-apu]) if(sims.type[j]=="normal") sims.type.false[j] = 1 				
				rankL = rankL[apu] # the chosen
			}
			# if none of the simulatios is the most "very strange" nor "strange", take one at random
			else {
				apu2 = sample(1:length(apu), size=1)
				if(sims.type[rankL[apu[apu2]]]!="outside") sims.type[rankL[apu[apu2]]] = "contributing" # if the simulation has not already been classified as "outside" (or contributing), give it the status "contributing", i.e. do not make change "outside" -> "contributing"
				sims.NofStrangeR[rankL[apu[apu2]]] = sims.NofStrangeR[rankL[apu[apu2]]] + 1       # increase the number "strange", i.e. contributing or outside, for this simulation
				# if the other simulations are "normal", i.e. inside envelopes for other r[k], classify them as "false strange"
				for(j in rankL[apu[-apu2]]) if(sims.type[j]=="normal") sims.type.false[j] = 1		
				rankL = rankL[apu[apu2]] # the chosen
			}
		}
		else {
			if(sims.type[rankL]!="outside") sims.type[rankL] = "contributing"  # if the simulation has not already been classified as "outside" (or contributing), give it the status "contributing", i.e. do not make change "outside" -> "contributing"
			sims.NofStrangeR[rankL] = sims.NofStrangeR[rankL] + 1                # increase the number "strange", i.e. contributing or outside, for this simulation
		}
		# the same for rankU
		if(length(rankU)!=1) {
			apu = which( sims.NofStrangeR[rankU] == max(sims.NofStrangeR[rankU]) )
			if(length(apu)==1) {
				if(sims.type[rankU[apu]]!="outside") sims.type[rankU[apu]] = "contributing"
				sims.NofStrangeR[rankU[apu]] = sims.NofStrangeR[rankU[apu]] + 1
				for(j in rankU[-apu]) if(sims.type[j]=="normal") sims.type.false[j] = 1 				
				rankU = rankU[apu] # the chosen
			}
			# if none of the simulatios is the most "very strange" nor "strange", take one at random
			else {
				apu2 = sample(1:length(apu), size=1)
				if(sims.type[rankU[apu[apu2]]]!="outside") sims.type[rankU[apu[apu2]]] = "contributing"
				sims.NofStrangeR[rankU[apu[apu2]]] = sims.NofStrangeR[rankU[apu[apu2]]] + 1
				for(j in rankU[apu[-apu2]]) if(sims.type[j]=="normal") sims.type.false[j] = 1		
				rankU = rankU[apu[apu2]] # the chosen
			}
		}
		else {
			if(sims.type[rankU]!="outside") sims.type[rankU] = "contributing"
			sims.NofStrangeR[rankU] = sims.NofStrangeR[rankU] + 1
		}
		envelopeL[i] = sims[i,rankL]
		envelopeU[i] = sims[i,rankU]
		outside = which( sim.ranksL < RANKLO | sim.ranksU > RANKUP )
		sims.type[outside] = "outside" 
		sims.NofStrangeR[outside] = sims.NofStrangeR[outside] + 1 
	}
	# ids of thosen simulations that are outside the envelopes
	ids = which(sims.type == "contributing" | sims.type=="outside")
	# Type I error approximation: typeIerror 
	
	# false sims.type.false
	# If some of the "false" are "outside"/"contributing" for smaller r[k], change "false" back to 0 ("not false")
	sims.type.false[sims.type!="normal"] = 0
	
	res = list(r=x$r[rranks], lower=envelopeL[rranks], upper=envelopeU[rranks], k=k, 
			typeIerror=length(ids)/NofSimulations,
			sims.type = sims.type, 
			sims.NofStrangeR = sims.NofStrangeR, sims.type.false=sims.type.false,
			call=match.call())
	oldClass(res) = "envs.typeIerror"
	res
}

print.envs.typeIerror = function(x, ...) {
	with(x, cat("Type I error: ", typeIerror, "\n"))	
	#print(table(x$sims.type))
}
