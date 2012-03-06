require(INLA)

#Load up the data - in the INLA package
data(Seeds)

# Define your formula. response ~ Latent field
#Random effects are formulated through the f( ) function

formula = r~ x1 + x2 +x1*x2 + f(plate,model="iid")

#Run INLA.  
#`family' defines the link function
#`Ntrials' - the number of trials (just for binomial)
#data - the data ;)
mod.seeds = inla(formula,family="binomial",Ntrials=n,data=Seeds)

#View the results 
summary(mod.seeds)

#We can make the inference better
hyp.seeds = inla.hyperpar(mod.seeds)
summary(hyp.seeds)
plot(hyp.seeds)
