source("simSIRmultiPopSingleIwithICUcomp.R")
intMat=readRDS("interactMatsymmetric.RDS")
popTable=readRDS("popParams.RDS")


#crushing
out=simulateSIRmultiPopICUcompartment( R0default = 2.7,R0init=0.9,  popParams=popTable, restriction.effect = 1, interactionMat = intMat, free.bracket.fraction = 1, free.bracket=17, forceNormalAt = round(365*1.25) , step.up.days=round(365*.5),nstep.max = 3*365, nstep.min = 3*365)

#plot mortality
plot(apply(out$data[, "M",],1,sum), type="l", ylim=c(0, 2.5e6))

#flatten
out=simulateSIRmultiPopICUcompartment( R0default = 2.7,R0init=1.6,  popParams=popTable, restriction.effect = 1, interactionMat = intMat, free.bracket.fraction = 1, free.bracket=17, forceNormalAt = round(365*1.25) , step.up.days=round(365*.5),nstep.max = 3*365, nstep.min = 3*365)
lines(apply(out$data[, "M",],1,sum), col=2)

#no interaction matrix
out=simulateSIRmultiPopICUcompartment( R0default = 2.7,R0init=1.6,  popParams=popTable, restriction.effect = 1, interactionMat = NULL, free.bracket.fraction = 1, free.bracket=17, forceNormalAt = round(365*1.25) , step.up.days=round(365*.5),nstep.max = 3*365, nstep.min = 3*365)
lines(apply(out$data[, "M",],1,sum), col=4) #will be off the scale



#heterogeneous
#relax through 50 (free.bracket=10), row 10 in popTable
#relax only 70% of the <50 population (free.bracket.faction)
#rest is at 10% interactions (restriction.effect=0.1)
#R0default and R0init are the same
out=simulateSIRmultiPopICUcompartment( R0default = 2.7,R0init=2.7,  popParams=popTable, restriction.effect = 0.3, interactionMat = intMat, free.bracket.fraction = 0.7, free.bracket=10, forceNormalAt = round(365*1.25) , step.up.days=round(365*.5),nstep.max = 3*365, nstep.min = 3*365)
lines(apply(out$data[, "M",],1,sum), col=3)
