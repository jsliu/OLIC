library(MASS)                                          # to use mvnorm()
dyn.load("../src/BCF.so")


# Simulating the changepoints and model choices
sim.bcf <- function(y,                                 # observations
                    x,                                 # variates
                    sample.size=100,                   # sample size of simulation
                    order=3,                           # model orders
                    delta=c(1e+2,1e+3,1e+4)^2,         # delta matrix
                    gama=2,nu=2,                  # parameters for prior dist
                    tranprob=.004,                     # transition probability
                    threshold=1e-6,                    # threshold for resampling
                    prior.model=c(0.5,0.5),            # prior probability of model
                    continuity=T,                      # if the model continuous or not
                    resample=T)                        # doing resampling
{
    if(length(delta)!=order)
      cat("The length of deltas should be consistant with order!\n")
    else
    {
      size <- length(y)
      .C("simBCF",
         as.double(y),as.double(x),
         as.integer(size),
         as.integer(sample.size),
         as.integer(order),
         as.double(delta),
         as.double(gama),
         as.double(nu),
         as.double(tranprob),
         as.double(threshold),
         as.double(prior.model),
         as.integer(continuity),
         as.integer(resample))
    }
  }

# plot the marginal distribution and return the simulation result
margin <- function(y,sample.size=100)
{
  size <- length(y)
  cpt <- scan("../output/ChangePoints")
  mod <- scan("../output/ModelAtChpts")
  
  cpts <- cpt[cpt!=0]
  mods <- mod[cpt!=0]

  ucpt <- sort(unique(cpts))
  umod <- sort(unique(mods))
  num.cpt <- length(ucpt)
  num.mod <- length(umod)

  cpt.prob <- tapply(cpts,cpts,length)/sample.size
  mod.prob <- numeric(length(ucpt))
  joint.prob <- matrix(nrow=num.mod, ncol=num.cpt)
  
  for(i in 1:length(umod))
    {
     for(j in 1:length(ucpt))
       {
         tmp <- mod[cpt==ucpt[j]]
         mod.prob[j] <- sum(tmp==umod[i])/length(tmp)
       }
     joint.prob[i,] <- cpt.prob*mod.prob
   }

  if(length(umod)==1)
    {
      plot(ucpt[-1]/size,joint.prob[,-1],type="h",xlim=range(0,1),ylim=range(0,1),
           xlab="Time", ylab="Joint Posterior",col=umod+1,lty=rep(1,length(umod)),
           main="Posterior of positions and types of change points")
    }
  else
    {
      matplot(ucpt[-1]/size,t(joint.prob[,-1]),type="h",xlim=range(0,1),ylim=range(0,1),
              xlab="Time", ylab="Joint Posterior",col=umod+1,lty=rep(1,length(umod)),
              main="Posterior of positions and types of change points")
      #  for(i in 1:length(umod))
      #      plot(ucpt[-1]/size,t(joint.prob[i,-1]),type="h",xlim=range(0,1),ylim=range(0,1),
      #          xlab="Time", ylab="Joint Posterior",col=i+1,lty=rep(1,length(umod)),
      #          main="Posterior of positions and types of change points")
    }

    return(list(changepoint=cpt,modelchoice=mod))
}


# With this function you can forget the code below
# which is used to estimate regression parameters
fit.plot <- function(x,y)
{
    data <- read.table("../output/FitData")
    n <- length(data)
    fit.data <- apply(data,2,mean)
    up.data <- apply(data,2,quantile,0.95)
    low.data <- apply(data,2,quantile,0.05)

    plot(x,y,pch=".")
    lines(x,fit.data,col=2)
    lines(x,up.data,col=3)
    lines(x,low.data,col=3)
}


########################################################################
#
#                      For continuous model
#
########################################################################

# Parameter estimation within exact inference
mix.para.est <- function(y,x,                   # response and variates
                         cptSample,             # simulation of changepoints
                         modSample,             # result of model choice
                         ord,                   # largest order
                         gama,nu,               # hyper prior
                         delta=rep(0,ord))      # variance of regression parameter
{
  n <- length(y)

  cptSample <- rev(c(n,cptSample))
  modSample <- rev(modSample)

  # looking for discontinuous change point
  # and count how many continuous change points are in between two discontinuous change points
  discpt <- 0
  num <- 0               # number of continuous change points
  num.concpt <- 0
  len <- length(cptSample)
  for(i in 2:len){
      num.concpt <- num.concpt+1
      if(modSample[i]==1 || cptSample[i]==n){
          discpt <- c(discpt,cptSample[i])
          num <- c(num,num.concpt)
          num.concpt <- 0
      }
  }

  beta.final <- NULL
  sigma2.final <- NULL
  len.discpt <- 0
  tmp <- 0
  for(i in 2:length(discpt)){
      # basis matrix, the number of rows is equal
      len.discpt <- c(len.discpt,discpt[i]-discpt[i-1])
      GS <- matrix(0,nrow=len.discpt[i],ncol=1+(ord-1)*num[i])
      GS[,1] <- 1
      num.cpt <- sum(num[1:(i-1)])
      pre.discpt <- sum(len.discpt[1:(i-1)])

      if(ord>1)
      {
          for(k in 1:num[i])
              for(m in 2:ord)
              {
                  GS[(cptSample[k+num.cpt]+1):cptSample[k+num.cpt+1]-pre.discpt,m+(k-1)*(ord-1)] <- (x[(cptSample[k+num.cpt]+1):cptSample[k+num.cpt+1]]-x[cptSample[k+num.cpt]+1])^(m-1)
                  if(cptSample[k+num.cpt+1]<discpt[i])
                      GS[(cptSample[k+num.cpt+1]+1):discpt[i]-pre.discpt,m+(k-1)*(ord-1)] <- (x[cptSample[k+num.cpt+1]]-x[cptSample[k+num.cpt]+1])^(m-1)
              }
      }

      # variances of beta
      if(ord==1)
      {
          delta.dash <- delta
          invDS <- diag(1/delta.dash,1,1)
      }
      else
      {
          delta.dash <- c(delta,rep(delta[-1],num[i]-1))
          invDS <- diag(1/delta.dash)
      }

      # observed data
      yS <- as.matrix(y[(discpt[i-1]+1):discpt[i]])

      MS <- solve(t(GS)%*%GS+invDS)
      PS <- diag(len.discpt[i])-GS%*%MS%*%t(GS)
      y.norm2 <- t(yS)%*%PS%*%yS

      # simulate sigma^2 from posterior probability first
      inv.sigma2 <- rgamma(1,(len.discpt[i]+nu)/2,1)/((y.norm2+gama)/2)
      sigma2.hat <- as.vector(1/inv.sigma2)

      # simulate beta given sigma^2 from posterior distribution
      ME <- MS%*%t(GS)%*%yS
      vbeta <- sigma2.hat*MS
      beta.hat <- mvrnorm(1,ME,vbeta)
      beta.dash <- matrix(0,nrow=num[i],ord)
      beta.dash[1,] <- beta.hat[1:ord]

      # calculate all the parameters
      if(num[i]>1){
          for(k in 2:num[i])
          {
              if(ord>1)
              {
                  beta.dash[k,-1] <- beta.hat[2:ord+(k-1)*(ord-1)]
                  beta.dash[k,1] <- GS[cptSample[k+num.cpt]+1-pre.discpt,2:(1+(k-1)*(ord-1))]%*%as.vector(t(beta.dash[1:(k-1),-1]))+beta.dash[1,1]
              }else{
                  beta.dash[k,1] <- beta.dash[1,]
              }
          }
      }

      beta.final <- rbind(beta.final,beta.dash)
      sigma2.final <- c(sigma2.final,sigma2.hat)
  }

  return(list(betas=beta.final,sigma2s=sigma2.final))
}


########################################################################
#
#                      For discontinuous model
#
########################################################################

# Parameter estimation
dist.para.est <- function(y,x,                   # response and variates
                          cptSample,             # simulation of change points
                          modSample,             # simulation of model choices
                          ord,                   # largest model order
                          gama,nu,               # hyper prior
                          delta=rep(0,ord))      # variance of regression parameter
{
  n <- length(y)
  cptSample <- rev(c(n,cptSample))
  modSample <- rev(modSample)

  beta.hat <- matrix(0,nrow=length(modSample),ncol=ord)
  sigma2.hat <- NULL

  for(i in 2:length(cptSample))
  {
      inter <- cptSample[i]-cptSample[i-1]
      mod <- modSample[i-1]

      # observations in a segment
      yS <- y[(cptSample[i-1]+1):cptSample[i]]

      # basic matrix in a segment
      GS <- matrix(0,nrow=inter,ncol=mod)
      for(m in 1:mod)
          GS[,m] <- (x[(cptSample[i-1]+1):cptSample[i]]-x[cptSample[i-1]+1])^(m-1)

      # variance matrix in a segment
      invDS <- diag(1/delta[1:mod],nrow=mod,ncol=mod)

      MS <- solve(t(GS)%*%GS+invDS)
      PS <- diag(inter)-GS%*%MS%*%t(GS)
      y.norm2 <- t(yS)%*%PS%*%yS

      # simulate sigma^2 from posterior probability first
      inv.sigma2 <- rgamma(1,(inter+nu)/2,(y.norm2+gama)/2)
      sigma2.dash <- 1/inv.sigma2
      sigma2.hat <- c(sigma2.hat,sigma2.dash)

      # simulate beta given sigma^2 from posterior distribution
      ME <- MS%*%t(GS)%*%yS
      vbeta <- sigma2.dash*MS
      beta.dash <- mvrnorm(1,ME,sigma2.dash*MS)
      beta.hat[i-1,1:mod] <- beta.dash
  }

  return(list(betas=beta.hat,sigma2s=sigma2.hat))
}


###################################################################################
#
#    Parameters simulation for both models
#
###################################################################################

# parameter simulation base on the sample of change points within exact inference
para.sim <- function(y,x,                           # response and variates
                     sim.out,                       # sample result
                     ord=3,                         # largest order of models
                     gama=2,nu=2,              # parameters for prior dist
                     delta=c(1e+2,1e+3,1e+4)^2,     # variance of regression parameter
                     continuity=T)
{
  plot(x,y,pch=".",xlab="Time",ylab="Simulation",main="Fitted curve")

  start <- 1                                        # starting point for change point
  end <- 1
  k <- 1

  n <- length(y)
  cpts <- sim.out$changepoint
  model <- sim.out$modelchoice
  y.hat <- matrix(0,nrow=sum(cpts==0),ncol=n)

  num.cpt <- NULL
  while(end<=length(cpts)){
    # to separate one simulation from the results which is cpts
    while(cpts[end]!=0)
      end <- end+1
    cptSample <- cpts[start:end]
    modSample <- model[start:end]
    SampLen <- length(cptSample)

    # record the number of change points
    num.cpt <- c(num.cpt,length(cptSample))
    
    if(continuity==T)
        out <- mix.para.est(y,x,cptSample,modSample,ord=ord,gama=gama,nu=nu,delta=delta)
    else
        out <- dist.para.est(y,x,cptSample,modSample,ord=ord,gama=gama,nu=nu,delta=delta)

    beta.hat <- out$betas
    sigma2.hat <- out$sigma2s
    #cat(sigma2.hat,"\n")
    # obtain the fitted value and draw the graph
    cptSample <- rev(c(n,cptSample))

    j <- 1
    for(i in 1:n){
        if(i>cptSample[j+1])
            j <- j+1

        xx <- (x[i]-x[cptSample[j]+1])^(1:length(beta.hat[j,])-1)
        y.hat[k,i] <- xx%*%beta.hat[j,]
    }

    #lines(x,y.hat[k,],col=k)
    #locator(1)
    start <- end+1
    end <- end+1
    k <- k+1
  }

  y.exa <- apply(y.hat,2,mean)
  y.up <- apply(y.hat,2,quantile,prob=0.95)
  y.low <- apply(y.hat,2,quantile,prob=0.05)

  lines(x,y.exa,col=2)
  #lines(x,y.up,col=3,lty=2)
  #lines(x,y.low,col=3,lty=2)
  #lines(x,y.true,col=4)

  return(list(y=y.hat,beta=beta.hat,sigma2=sigma2.hat,num.cpt=num.cpt))
}


# calculate the emprical prior of sigma^2, i.e. gamma and nu
# by moving average

empirical <- function(y,neighbour=2,init=2)
{
    size <- length(y)
    y.guess <- rep(0,size)

    for(i in 1:size)
        if(i<=neighbour)
            y.guess[i] <- median(y[1:(i+neighbour)])
        else if(i>size-neighbour)
            y.guess[i] <- median(y[(i-neighbour):size])
        else
            y.guess[i] <- median(y[(i-neighbour):(i+neighbour)])

    #plot(y,pch=".")
    #lines(c(1:length(y)),y.guess,col=2)
    sigma.guess <- mean((y-y.guess)^2)

    alpha <- init
    beta <- alpha*sigma.guess

    return(list(nu=alpha, gama=beta))
}
