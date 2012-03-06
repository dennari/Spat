require(INLA)
data(Leuk)
g = system.file("demodata/Leuk.graph", package="INLA")

formula = inla.surv(Leuk$time, Leuk$cens) ~ sex + age +
    f(inla.group(wbc), model="rw1")+
    f(inla.group(tpi), model="rw2")+
    f(district, model="besag", graph.file = g)

result = inla(formula, family="coxph", data=Leuk)

source(system.file("demodata/Leuk-map.R", package="INLA"))
Leuk.map(result$summary.random$district$mean)
plot(result)