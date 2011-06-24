# run this code in ../../R/
source("sf1.R")   # change here
#source("sf2.R")   # change here
source("simBCF.R")

sim.out <- sim.bcf(y,x,gama=2,nu=2,sample.size=100)
#margin.out <- margin(y,sample.size=100)
#out <- para.sim(y,x,margin.out,gama=2,nu=2)

par(mfrow=c(1,2))
x <- seq(0,1,length=200)
margin.out <- margin(y,sample.size=100)
out <- para.sim(y,x,margin.out,gama=2,nu=2)
y.hat <- out$y
y.exa <- apply(y.hat,2,mean)
y.up <- apply(y.hat,2,quantile,prob=0.95)
y.low <- apply(y.hat,2,quantile,prob=0.05)

plot(x,y,pch=".",xlab="Time",ylab="Simulation",main="True & fitted curves")
lines(x,y.exa,col=2)
lines(x,y.true,col=4)
