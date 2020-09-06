#' Simulates the COVID-19 dynamics of interacting populations, with a single infected population (exponential infectious period) and a seperate ICU compartment.
#' 
#' @param R0default The default R0 value wih no mititgations
#' @param R0init The initial R0 value for the unrestricted population. Useful for modeling homogeneous mitigations. 
#' @param popParamsIn A data.frame of population bracket parameters.  Must contain "age" (character name for the age bracket),"hosp" (fraction of cases needing hospitalization), "crit" (fraction of hospitalizations needing ICU), "mort" case mortality, "fraction" (fraction of the population in that bracket). Example in populationParams.txt. 
#' @param restriction.effect=.3 The mulitplier for interaction frequency applied to the contact matrix when mitigations are present
#' @param free.bracket The bracket ID for the unrestricted population. All brackets up to free.bracket are unrestricted at time 0.
#' @param free.bracket.fraction=1 the fraction of the population in brackets up.to free.bracket that is unrestricted. Some inidividuals are still restricted because of risk factors other than age.
#'  @param interactionMat The unitless, symmetric interaction frequency matrix for population brackets defined in popParamsIn. Example in interactionMatrix.txt.
#'  @param step.up.days Length in days of the ramp up period to normal transmision
#'  @param forceNormalAt Time in days when normal transmission resumes
#'  @param  step.size Step size in days
#'  @param nstep.max Maximum number of steps
#'  @param nstep.min Minimum number of steps
#'  @param init.vec The vector of initial infected individuals
#'  @param beta0 Hard coded beta0 (otherwise it is computed to match R0default)
#'
#'



simulateSIRmultiPopICUcompartment=function(R0default,  R0init, popParamsIn, restriction.effect=.3, free.bracket=4, free.bracket.fraction=1, interactionMat=NULL,  step.up.days=90, forceNormalAt=NULL,  step.size=1, nstep.max=365*3, nstep.min=NULL,  init.vec=NULL, beta0=NULL){
  
  if(is.null(interactionMat)){ #default interaction matrix is just all 1s, all interactions are equally likely
    message("Constant interaction matrix")
    interactionMat=matrix(1, nrow=nrow(popParamsIn), ncol=nrow(popParamsIn)) 
  }
  #number of populations we are modeling
  npop=nrow(popParamsIn)
  
  #hard coded parameters"
  I=100000 #total number of infected individuals
  year=365 #year
  recoverytime=14 #recovery time in days. Taken to be 14 day
  N=totalPop=324356000 #US population, last census
  
  
  iif=1:free.bracket #index of the free age bracket
  iir=(free.bracket+1):npop #index of the restricted bracket
  
  #if free.bracket fraction is less than 1 we have to split those age groups in two
  if(free.bracket.fraction<1){ 
    #set up expanded populations
    popParamsInFreeR=popParamsInFree=popParamsIn[iif,]
    popParamsInFreeR$age=paste(popParamsInFreeR$age, "R")
    tmpfrac=popParamsInFree$fraction
    #redistribute the population fraction  across free and free-R (technically free but restricted due to factors other than age)
    popParamsInFree$fraction=tmpfrac*free.bracket.fraction 
    popParamsInFreeR$fraction=tmpfrac*(1-free.bracket.fraction) 
    #this is the restricted brackets
    popParamsInR=popParamsIn[iir,]
    #the final popParamsIn is a combination of free, free-R, and R
    popParamsInThis=rbind(popParamsInFree, popParamsInFreeR,popParamsInR)
    
    #expand the interaction matrix in the same way
    #the proportions are facored out of the interaction matrix so there is need to do any adjustments, we can just copy
    A=interactionMat[iif,iif]
    B=interactionMat[iif, iir]
    C=interactionMat[iir, iir]
    
    interactionMat=rbind(cbind(A,A,B), cbind(A,A,B), cbind(t(B),t(B),C))
    rownames(interactionMat)=colnames(interactionMat)=popParamsInThis$age
    
  }
  else{ #otherwise popParamsIn and interactionMat can be left as is
    popParamsInThis=popParamsIn
  }
  
  #adjust the number of population brackets
  npop=nrow(popParamsInThis)
  
  #total number of individuals in each bracket
  N0Vec=N*popParamsInThis$fraction
  
  
  
  
  
  #constants for leaving various states
  alpha=1/recoverytime #total rate of leaving the infected state
  leave.critical=1/(10/14*recoverytime) #hard coded to be 10 when recovery is 14
  leave.critical.hold=1/(4/14*recoverytime) # individuals that end up in C spends 4 days as infectius
  
  
  #compute the default beta corresponding to R=Rdefault
  Qmat=interactionMat%*%diag(popParamsInThis$fraction) # age distribution adjusted interaction matrix
  #the initial "steady state" distribution and the R0 multiplier are eigenvectors and eigenvalue of Q matrix
#left eigenvalues because I is multipled on the left
    ee=eigen(t(Qmat))
  ss=Re(ee$vectors[,1])
  ss=ss*sign(ss[1]) #make positive
  ss=ss/sum(ss) #normalize to sum=1
  if(is.null(beta0)){
  beta0=R0default*alpha/(Re(ee$values[1]))
  }
  else{
    message("beta0 is given")
  }
  
  #in order to adjust from R0 default we compute multiplicative adjustments
  #this is for global mitigations in the initial phase
  multInit=R0init/R0default
  
  
  
  #distribute infectoins according to steady state
  if(is.null(init.vec)){
    Ivec=I0Vec=I*ss
  }
  else{
    Ivec=I0Vec=init.vec
  }

  #setting up initial conditions
  Svec=N0Vec-I0Vec #number susceptible
  Rvec=double(length(Svec))
  Mvec=double(length(Svec))
  
  
  #adjust the interaction matrix to account for mitigations
  interactionMatFree=interactionMat #interactionMatFree*beta0 will give R0default transmission
  interactionMat=interactionMatFree*restriction.effect*multInit #initial transmission
  interactionMat[iif, iif]=interactionMatFree[iif,iif]*multInit #release the free bracked [indexed by iff]
  interactionMatR=interactionMat #save the restricted interaction matrix for output, the interactionMat itself will change
  
  
  #so the output will be in order
  # S, R, I, M, Chold, C for each age group
  
  output=array(data = 0,dim=c(nstep.max+1, 7, nrow(popParamsInThis)), dimnames = list(1:(nstep.max+1), c("S", "R", "I", "M", "Chold", "C", "Itotal"), popParamsInThis$age))
  
  #initialize infected and susceptible, everything else is 0
  output[1,"S",]=Svec
  output[1,"I",]=Ivec
  

  output[1, "Chold",]=rep(0, nrow(popParamsInThis))
  output[1, "C",]=rep(0, nrow(popParamsInThis))
  peak.reached=F
  peak.val=-1
  restrictions=T
  s=0
  rlifted=Inf
  smin=0
  
  if(is.null(nstep.min)){
    smin=forceNormalAt*1.25
  }
  else{
    smin=nstep.min
  }
  
  #not keeping track of hospitalizations, translate input to overal critical care rate.
  popParamsInThis$crit=popParamsInThis$crit*popParamsInThis$hosp
  
  
  #stopping criteria is number of infections < 10000 or nstep.max is reached
  while((sum(output[dim(output)[1],"I",])>10000 && s<nstep.max)||s<smin){
    s=s+1 
    day=s*step.size
    
    
    
    
    curTime=output[s,,] #get the current slice
    if(any(is.na(curTime))){
      stop()
    }
    
    #compute the next step according to the SIR model
    
    
    #  newI=curTime["I",]+((tmpIplus<-rowSums(outer(curTime["I",], curTime["S",])*beta0/N))-alpha*curTime["I",])*step.size
#people in Chold are still infectious
  #  newI=curTime["I",]+((
    
    #compute the temporary newly infected people
    tmpIplusAll=((curTime["I",]+curTime["Chold",])%*%(interactionMat)%*%diag(curTime["S",])*beta0/N)
    #split them up into I and Chold
    tmpCholdPlus=tmpIplusAll*popParamsInThis$crit
    tmpIplus=tmpIplusAll-tmpCholdPlus
    newI=curTime["I",]+(tmpIplus-alpha*curTime["I",])*step.size
    
    #update the other variables
    newS=curTime["S",]-tmpIplusAll*step.size
newChold=curTime["Chold",]+(tmpCholdPlus-curTime["Chold",]*leave.critical.hold)*step.size
newC=curTime["C",]+(curTime["Chold",]*leave.critical.hold-curTime["C",]*leave.critical)*step.size
    newR=curTime["R",]+(curTime["I",]*alpha+curTime["C",]*leave.critical)*step.size
    newM=newR*popParamsInThis$mort
    
    
    
      
    #handle the forced normal point
    
    #begin lifting restrictions
    if(abs(forceNormalAt-day-step.up.days)<1e-6){
      rlifted=s
      message(paste("Starting at", s))
    }
  
    #ramp up perriod
    if(forceNormalAt-day<step.up.days && restrictions==T){
      
    #  message(paste("Force normal step up", day))
     
      daydiff=day-rlifted
      interactionMat=interactionMatFree*(multInit*restriction.effect*(1-daydiff/step.up.days)+(daydiff/step.up.days)) #ramp up from init*restriction to release
      
      interactionMat[iif,iif]=(interactionMatFree*(multInit*(1-daydiff/step.up.days)+(daydiff/step.up.days)))[iif,iif] #ramp up from init to 1
      
      #show(interactionMat[8,8])
      restrictions=T
    }
    
    if(abs(day-forceNormalAt)<1e-6){
      message("Force normal final")
      interactionMat=interactionMatFree
      interactionMat[iif,iif]=interactionMatFree[iif,iif] #set0
      restrictions=F
    }
    
  
  #create the next time step in order
  #c("S", "R", "I", "M", "Chold", "C", "Itotal")
  # variables Chold and C are not part of the model, they are computed from newI
  nextTime=rbind(newS, newR, newI, newM, newChold, newC, newI+newC+newChold)
  rownames(nextTime)=c("S", "R", "I", "M", "Chold", "C", "Itotal")
#  output=abind(output, nextTime, along=1)    #append to output
output[s+1,,]=nextTime
  }
  
  #change the names so plotting is consistent
  dimnames(output)[[2]]=c("S", "R", "InotC", "M", "Chold", "C", "I")
  
#show(c(max(apply(output[, "C", ],1,sum)), totalICU))
#save the relaxation time
stopifnot(s>=smin)
return(list(data=output, beta0=beta0, eigenvec=ss,rlifted=rlifted, interactionMatR=interactionMatR, interactionMatFree=interactionMatFree, free.bracket=free.bracket, free.bracket.fraction=free.bracket.fraction, restriction.effect=restriction.effect, R0default=R0default, R0init=R0init,  step.up.days=step.up.days))

}
