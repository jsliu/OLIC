################################################################
#
# The importance weights for each simulation, i.e
#                  P(cpt,q|data)/P.hat(cpt,q|data)
# where the P and P.hat are joint posterior of change points
#
#    P(cpt,q|data) \propto P(data|cpt,q)P(cpt)P(q)
#
################################################################

# Fisrtly, calculate the exact posterior P(cpt|data)
exact.post <- function(y,x,                   # response and variates
                       cptSample,             # simulation of changepoints
                       modSample,             # result of model choice
                       tran,                  # transition probability
                       ord,                   # largest order
                       gama,nu,               # hyper prior
                       delta=rep(0,ord))      # variance of regression parameter
{
    # calculate log likelihood
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
    sigma2.final <- NULL
    len.discpt <- 0
    tmp <- 0
    loglik <- 0

    for(i in 2:length(discpt))
    {
        # To calculate logP(s,t,q)
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
        }else{
            delta.dash <- c(delta,rep(delta[-1],num[i]-1))
            invDS <- diag(1/delta.dash)
        }

        # observed data
        yS <- as.matrix(y[(discpt[i-1]+1):discpt[i]])

        MS <- solve(t(GS)%*%GS+invDS)
        PS <- diag(len.discpt[i])-GS%*%MS%*%t(GS)
        y.norm2 <- t(yS)%*%PS%*%yS

        part1 <- log(det(MS))-sum(log(delta.dash))
        part2 <- nu*log(gama)-(len.discpt[i]+nu)*log(y.norm2+gama)
        part3 <- lgamma((len.discpt[i]+nu)/2)-lgamma(nu/2)

        logP <- (-len.discpt[i]*log(pi)+part1+part2)/2+part3
        loglik <- loglik+logP
    } # for

    logpri <- (len-2)*log(tran)+(n-(len-2)-1)*log(1-tran)+(len-1)*log(0.5)
    logpost <- loglik+logpri
    #cat(loglik,"\t",logpri,"\n")
    return(logpost)
}

# Secondly the joint posterior by approximated model
# for the independence assumption
apprx.post <- function(y,x,cptSample,modSample)
{
    size <- length(y)
    particle <- read.table(file="../output/CptPosAppr",header=F,sep="")

    len <- nrow(particle)
    new.table <- cbind(c(1:len),particle)
    index <- new.table[is.na(new.table[,2,]),1]

    posterior <- 1
    ptr <- 1
    i <- 1

    while(ptr<=size)
    {
        if(ptr<size)
            subpart <- new.table[(index[ptr]+1):(index[ptr+1]-1),-1]
        else
            subpart <- new.table[(index[ptr]+1):nrow(particle),-1]

        tmp<- subpart[subpart[,2]==cptSample[i],]
        posterior <- tmp[tmp[,3]==modSample[i],1]*posterior
        ptr <- size-cptSample[i]+1
        #cat(cptSample[i],"\t",modSample[i],"\t")
        #cat(tmp[tmp[,3]==modSample[i],1],"\n")
        i <- i+1
    }

    return(log(posterior))
}


# Calculate the imporatance sampling weights
ISweight <- function(y,x,cpts,mods,
                     tran=0.004,
                     ord=3,
                     gama=1000,nu=1000,
                     delta=rep(0,ord))
{
    size <- length(y)
    lg.epost <- NULL
    lg.apost <- NULL

    start <- 1
    end <- 1
    k <- 1

    weis <- scan("../output/CptPosAppr")
    
    while(end<=length(cpts))
    {
        #cat(k,"\n")
        # to separate one simulation from the results which is cpts
        while(cpts[end]!=0)
            end <- end+1

        cptSample <- cpts[start:end]
        modSample <- mods[start:end]
        weiSample <- weis[start:end]

        lg.epost <- c(lg.epost,exact.post(y,x,cptSample,modSample,tran,ord,gama,nu,delta))
        lg.apost <- c(lg.apost,log(prod(weiSample))) 
        #lg.apost <- c(lg.apost,apprx.post(y,x,cptSample,modSample))
        if(k%%10==0)
            cat(k,"\n",cptSample,"\nlog(epost)=",lg.epost[k],"\nlog(apost)=",lg.apost[k],"\n")

        start <- end+1
        end <- end+1
        k <- k+1
    }

    plot(c(1:(k-1)),lg.epost-lg.apost,type="l",col=2)
    return(list(elgpost=lg.epost,algpost=lg.apost))
}


# normalising the weights
normalise <- function(logprobs)
{
    block1 <- max(logprobs[[1]])
    lg.epost <- logprobs[[1]]-block1
    epost <- exp(lg.epost)
    epost <- epost/sum(epost)

    apost <- exp(logprobs[[2]])
    apost <- apost/sum(apost)
    wei <- epost/apost
    #wei <- wei/sum(wei)

    #plot(wei,type="l",col=2)
    hist(wei,breaks=100,xlab="Importance Weights",main="")
    return(wei)
}


# effective sample size
curve.res <- function(y.sample,wei)
{
    N <- nrow(y.sample)
    wei <- wei/sum(wei)
    index <- sample(1:N,prob=wei,replace=T)

    y.res <- y.sample
    for(i in 1:N)
        y.res[i,] <- y.sample[index[i],]

    y.fit <- apply(y.res,2,mean)
    plot(y.fit,type="l",col=2)
    return(y.fit)
}
