# block function
blocks <- function(x,noise=1,std=1)
{
    sig <- 2.839887
    cpts <- c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
    hs <- c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)

    size <- length(x)
    ncpt <- length(cpts)

    y <- numeric(size)
    for(i in 1:size)
    {
        for(j in 1:ncpt)
            y[i] <- y[i]+sig*0.5*hs[j]*(1+sign(x[i]-cpts[j]))

        y[i] <- y[i]+noise*rnorm(1,sd=std)
    }

    return(y)
}
#sim.out <- sim.bcf(y,x,order=1,delta=100^2,sample.size=100)
#margin.out <- margin(y,sample.size=100)
#out <- para.sim(y,x,y.true,margin.out,ord=1,delta=100^2)

#sim.out <- sim.bcf(y,x,sample.size=100,continuity=F,prior.model=rep(1/3,3))
#margin.out <- margin(y,sample.size=100)
#out <- para.sim(y,x,y.true,margin.out,continuity=F)

# bump function
bumps <- function(x,noise=1,std=1)
{
    sig <- 9.697620
    cpts <- c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
    hs <- c(4,5,3,4,5,4.2,2.1,4.3,3.1,5.1,4.2)
    ws <- c(0.005,0.005,0.006,0.01,0.01,0.03,0.01,0.01,0.005,0.008,0.005)

    size <- length(x)
    ncpt <- length(cpts)

    y <- numeric(size)
    for(i in 1:size)
    {
        for(j in 1:ncpt)
            y[i] <- y[i]+sig*hs[j]*(1+abs((x[i]-cpts[j])/ws[j]))^(-4)

        y[i] <- y[i]+noise*rnorm(1,sd=std)
    }

    return(y)
}


# heavisine function
heavi.sine <- function(x,noise=1,std=1)
{
    size <- length(x)

    sig <- 2.268244

    f <- sig*(4*sin(4*pi*x)-sign(x-.3)-sign(.72-x))+noise*rnorm(size,sd=std)
    #f <- sig*(4*sin(4*pi*x)-sign(x-.3))+noise*rnorm(size,sd=std)
    return(f)
}


#doppler function
doppler <- function(x,noise=1,std=1)
{
    size <- length(x)

    sig <- 23.889505;

    y = sig*(x*(1-x))^(0.5)*sin(2.0*pi*1.05/(x+0.05)) + noise*rnorm(size,sd=std);
}

# smooth function
sm <- function(x,noise=1)
{
    beta <- matrix(nrow=5,ncol=3)
    beta[1,] <- c(0,10,-40)
    beta[2,] <- c(1.985,-10,10)
    beta[3,] <- c(0.83,-5,5)
    beta[4,] <- c(-1.27,2,-0.8)
    beta[5,] <- c(-2.104,4,-2)

    size <- length(x)
    cpts <- c(90,180,320,410,size)/size
    y <- numeric(size)

    j <- 1
    ord <- 3
    sigma <- 0.09

    for(i in 1:size)
    {
        if(x[i]>cpts[j])
           j <- j+1

        y[i] <- (x[i]^(1:ord-1))%*%beta[j,]+noise*rnorm(1,0,sigma)
    }

    return(y)
}


# other smooth functions in Denison/Mallick/Smith, JRSS B, 60, 1998
# x <- seq(-2,2,length=200)
sf1 <- function(x,noise=1)
{
    size <- length(x)
    f <- x+2*exp(-16*x^2)+noise*rnorm(size,0,.4)
    return(f)
}

#sim.out <- sim.bcf(y,x,gama=100,nu=1000)
#out <- para.sim(y,x,sim.out,gama=100,nu=1000)

# x <- seq(-2,2,length=200)
sf2 <- function(x,noise=1)
{
    size <- length(x)
    f <- sin(2*x)+2*exp(-16*x^2)+noise*rnorm(size,0,.3)
    return(f)
}

#sim.out <- sim.bcf(y,x,gama=50,nu=1000)
#out <- para.sim(y,x,sim.out,gama=50,nu=1000)

# function 1,2,3 in Mao/Zhao JRSS B, 65, 2003
g2 <- function(x,noise=1)
{
    n <- length(x)
    f <- 1.5*dnorm((x-0.35)/0.15)-dnorm((x-0.8)/0.04)+noise*rnorm(n,0,0.054)
    return(f)
}

g3 <- function(x,noise=1)
{
    n <- length(x)
    f <- numeric(n)
    for(t in 1:n)
    {
        if(x[t]<0.4079)
            f[t] <- 3*(3*(x[t]-0.2)^2+0.5)
        else if(x[t]<0.8)
            f[t] <- 3*(-1.2*(x[t]-0.65)^2+0.7)
        else
            f[t] <- 3*(-1.2*(x[t]-0.65)^2+0.7-0.07)

        f[t] <- f[t]+noise*rnorm(1,0,0.054)
    }
    return(f)
}

#sim.out <- sim.bcf(y,x,gama=3,nu=1000)
#out <- para.sim(y,x,sim.out,gama=3,nu=1000)
