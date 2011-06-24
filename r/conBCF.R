#######################################################################
#                                                                     #
#  Almost identical function. Replace linear regression by piecewise  #
#  polynomials which is identical to Denison, et al 1998              #
#                                                                     #
#######################################################################

library(MASS)                                          # to use mvnorm()
dyn.load("BCF.so")


# Simulating the changepoints and model choices
sim.bcf <- function(y,                                 # observations
                    x,                                 # variates
                    sample.size=1000,                  # sample size of simulation
                    order=3,                           # model orders
                    delta=c(1e+2,1e+3,1e+4)^2,         # delta matrix
                    gama=1000,nu=1000,                 # parameters for prior dist
                    tranprob=.004,                     # transition probability
                    threshold=1e-6,                    # threshold for resampling
                    prior.model=c(.5,.5),              # prior probability of model
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

        cpt <- scan("../output/ChangePoints")
        mod <- scan("../output/ModelAtChpts")

        cpts <- cpt[cpt!=0]
        mods <- mod[cpt!=0]

        ucpt <- sort(unique(cpts))
        umod <- sort(unique(mods))
        cpt.prob <- tapply(cpts,cpts,length)/sample.size
        xcpt <- as.numeric(names(cpt.prob))
        mod.prob <- numeric(length(ucpt))
        joint.prob <- matrix(nrow=length(umod), ncol=length(ucpt))

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
        }
        return(list(changepoint=cpt,modelchoice=mod))
    }  # if-else
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
    y.fit <- y
    n <- length(y)

    cptSample <- rev(c(n,cptSample))
    modSample <- rev(modSample)

    # looking for discontinuous change point
    # and count how many continuous change points are between two discontinuous change points
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
                        # this is the case of cubic spline
                        GS[(cptSample[k+num.cpt+1]+1):discpt[i]-pre.discpt,m+(k-1)*(ord-1)] <- (x[(cptSample[k+num.cpt+1]+1):discpt[i]]-x[cptSample[k+num.cpt]+1])^(m-1)
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

        y.fit[(discpt[i-1]+1):discpt[i]] <- GS%*%beta.hat
    }  # for i

    return(y.fit)
}



###################################################################################
#
#    Parameters simulation for continuous models
#
###################################################################################

# parameter simulation base on the sample of change points within exact inference
para.sim <- function(y,x,                               # response and variates
                     sim.out,                           # sample result
                     ord=3,                             # largest order of models
                     gama=1000,nu=1000,                 # parameters for prior dist
                     delta=c(1e+2,1e+3,1e+4)^2)         # variance of regression parameter
{
    plot(x,y,pch=".",xlab="Time",ylab="Simulation",main="Fitted curve")

    start <- 1                                        # starting point for change point
    end <- 1
    k <- 1

    n <- length(y)
    cpts <- sim.out$changepoint
    model <- sim.out$modelchoice
    y.hat <- matrix(0,nrow=sum(cpts==0),ncol=n)

    while(end<=length(cpts))
    {
        # to separate one simulation from the results which is cpts
        while(cpts[end]!=0)
            end <- end+1

        cptSample <- cpts[start:end]
        modSample <- model[start:end]

        y.hat[k,] <- mix.para.est(y,x,cptSample,modSample,ord=ord,gama=gama,nu=nu,delta=delta)

        start <- end+1
        end <- end+1
        k <- k+1
    }

    y.exa <- apply(y.hat,2,mean)
    y.up <- apply(y.hat,2,quantile,prob=0.95)
    y.low <- apply(y.hat,2,quantile,prob=0.05)

    #cin <- sum(y.true<=y.up & y.true>=y.low)
    #cat("\nProportion of times within Credible Interval =",cin/n,"\n\n")

    lines(x,y.exa,col=2)
    lines(x,y.up,col=3)
    lines(x,y.low,col=3)

    return(y.hat)
}







