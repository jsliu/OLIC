/***********************************************
 *
 *  The code is originated from Paul Fearnhead
 *  I wrote some function in C++
 *
 ***********************************************/


#include <cstdio>
#include <cmath>
using namespace std;

#include "rand.hpp"

/*********************************
 *                               *
 * Routine to read seed for file *
 *                               *
 *********************************/

/* global variables */
FILE *seed_file;

/* seed */
unsigned long seed[2];

void seed_set(FILE *file)
{
  long i,temp;

  unsigned long max=4294967295ul;
  char ch;
 
  for(i=0;i<2;i++){
    /*ignore white space*/
    
    while(1){
      ch = fgetc(file);
      if(ch != '\n' && ch != '\n' && ch != '\t') break;
    }
    
    temp= ch-'0';
    
    while(1){
      ch = fgetc(file);
      if(ch == '\n' || ch == '\n' || ch == '\t' || ch ==EOF) break;
      temp= (10*temp+ch-'0') % max;
    }
    *(seed+i)=temp;
  }
}

/***********************
 *                     *
 * put seed in file    *
 * file_r and file_w   *
 * point to same file  *
 * but for read (r)    *
 * and write (w)       *
 *                     *
 ***********************/

void seed_put(FILE *file)
{
  (void)fprintf(file,"%lu\n%lu", *seed,  *(seed+1));
}   

/* set the random generation enviorment */
void getRNG(void)
{
  /* obtain seed */
  seed_file = fopen(".seed","r");
  if(seed_file==NULL){
    (void) fprintf(stderr,"Can't open .seed; Using default values.\n");
    *seed = 1349761;
    *(seed+1) = 2498265;
  }else{
    seed_set(seed_file);
    (void)fclose(seed_file);
  }
}

void putRNG(void)
{
   /* change .seed file */
  seed_file = fopen(".seed","w");
  seed_put(seed_file);
  (void)fclose(seed_file);
}


/* generate a binomially distributed r.v. with parameters n and p */
int rbinom(int n, double p)
{
#define KK 10
  int i,k,X;
  double theta;
  double V;

  k=n;
  X=0;
  theta=p;
  while(k > KK){
    i=(int)(1+k*theta);
    V=rbeta((double) i, (double)(k+1-i) );
    if(theta<V){
      theta=theta/V;
      k=i-1;
    }else{
      X=X+i;
      theta=(theta-V)/(1-V);
      k=k-i;
    }
  }
  for(i=0;i<k;i++)
      if(runif()<theta) X++;
  
  return X;
}

/*generate a beta(a,b) r.v. a beta is gamma(a)/(gamma(a)+gamma(b))*/
double rbeta(double a, double b)
{
  double X,Y;
  X=rgamma(a);
  Y=rgamma(b);
  return X/(X+Y);
}

/*S+ runif routine*/
double runif(void)
{
        unsigned long n, lambda=69069;
        double temp;
        *seed = *seed * lambda;
        *(seed+1) ^= *(seed+1) >> 15;
        *(seed+1) ^= *(seed+1) << 17;
        n = *seed ^ *(seed+1);
        temp=(((n>>1) & 017777777777)) / 2147483648. ;
        temp+= 1/4294967296.0; /*never return 0 or 1*/
        return(temp);
}

/*standard normal random variable*/
/*Box-Muller*/
double stnorm(void)
{
  double E,R,TH;

  TH=6.2831853*runif();
  E=-log(runif());
  R=sqrt(2.0*E);
  return(R*cos(TH));
 
}


/*gamma random variable - shape a */

double rgamma(double a)
{
  double U1,U2;
  double e=2.718281828;
  double X;
  
  if(a<=1){
    double c=e/(a+e);
    while(1){
      U1=runif();
      U2=runif();
      if(U1< c ){
	X=log( U1/c )/a;
	if(log(U2)< -exp(X) ) return(exp(X));
      }
      else{
	X= 1-log( (1-U1)/(1-c) );
	if(log(U2)< (a-1)*log(X)) return(X);
      }
    }
  }
  else{
    double c1=a-1.0;
    double c2=(a-1.0/(6.0*a))/c1;
    double c3=2.0/c1;
    double c4=2.0/(a-1)+2.0;
    double c5=1/sqrt(a);
    double W;

    while(1){
      U1=-1;
      while(U1<0 || U1>1){
	U1=runif();
	U2=runif();
	if(a>2.5) 
	  U1=U2+c5*(1.0-1.86*U1);
      }
      W=c2*U2/U1;
      if(c3*U1+W+1/W <= c4) 
	return(c1*W);
      if(c3*log(U1)-log(W)+W <1) 
	return(c1*W);
    }
  }
}

/* sample from a discrete density - unnormalized weights w */
/* returns index of sampled value, 0,1,..M-1 */
/* sample of size N, use uniforms U WHICH ARE ORDERED */
void rdunif(double *w,double sumw,int *index,double *U,int N,int M)
{
  int i=-1;
  int j= 0;
  double sum = 0;
  while(j<N){
    i++;
    if(i==M){
      (void)fprintf(stderr,"Error in rdunif; run out of weights\n");
      abort();
    }
    sum+=w[i];
    while(U[j]*sumw<sum && j<N){
      index[j]=i;
      j++;
    }
  }
}


/* calculates the complementary error function (error 10^-7) -
  see Press et al. p221. 
  cerf(x)=2/(root(pi)) * int_{x to infinity} exp (-t^2) dt */
double erfcc(double x)
{
  double t,z,ans;

  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 + t*(-0.18628806 + t*(0.27886807)+t*(-1.13520398 + t*(1.48851587 + t*(-0.82215223 + t*(0.17087277)))))))));
  return x>=0 ? ans : 2.0 - ans;
}

/* calculate error function using series method - see Press et al. p 218ff */

double erffc(double x)
{
 return x <0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);

}

/* calculate incomplete gamma function */

double gammp(double a, double x)
{
  double gamser,gammcf,gln;

  if(x<0.0 || a<=0.0){
    fprintf(stderr,"Invalid argument in routine gammp \n");
    abort();
  }
  if(x< (a+1.0)){
    gser(&gamser,a,x,&gln);
    return gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf;
  }
}
 
/* returns incomplete gamma function evaluate by series representation */
void gser(double *gamser,double a, double x, double *gln)
{
  int n;
  double sum,del,ap;

  if(a!=0.5)
    *gln=gammln_p(a);
  else
    *gln=0.572364942; /*log of root pi*/
  if(x<=0.0){
    if(x<0.0){
      fprintf(stderr,"x less than 0 in routine gser\n");
      abort();
    }
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for(n=1;n<=ITMAX;n++){
      ++ap;
      del *= x/ap;
      sum += del;
      if( fabs(del) < fabs(sum)*EPS){
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    fprintf(stderr,"a too large, ITMAX too small in gser\n");
    *gamser=sum*exp(-x+a*log(x)-(*gln));
    return;
  }

}

/* returns the incomplete gamma function evaluated by continued fraction */
void gcf(double *gammcf,double a, double x, double *gln)
{
  int i;
  double an,b,c,d,del,h;

  if(a!=0.5)
    *gln=gammln_p(a);
  else
    *gln=0.572364942;
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for(i=1;i<ITMAX;i++){
    an= -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if ( fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if( fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) 
    fprintf(stderr,"a too large, ITMAX too small in gcf\n");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;

}

/* returns log(gamma_fn(xx) for xx>0 */
double gammln_p(double xx)
{
  double x,y,tmp,ser,ans;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155, 0.1208650973866179e-02,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for(j=0;j<=5;j++) ser += cof[j]/++y;
  ans=-tmp+log(2.5066282746310005*ser/x);
  return( ans );
}

/*********************************************
 *
 *    The C++ code thereafter
 *
 *********************************************/

/* generate a negative binomial distribution with parameter size and pr 
   where size is the target for number of successful trials */
int rnbinom(int size, double pr)
{
  double uu;
  int *u = new int[size];

  int nb = 0;
  for(int i=0;i<size;i++){
    uu = runif();
    u[i] = (int) floor(log(uu)/log(1-pr))+1;
    nb += u[i];
  }

  delete []u;
  return nb;
}


/* generate a multivariate normal distribution, with a mean mu
   and a variance matrix sigma */
Matrix<double> mult_norm(const Matrix<double> &mean, const Matrix<double> &var)
{
		int dim=mean.num_rows();
		Matrix<double> Mnorm(dim,1);
		for(int i=0;i<dim;i++)
				Mnorm[i][0] = stnorm();

		Cholesky<double> cho_var(var);

		Matrix<double> L;
		L = cho_var.getL();

		Matrix<double>tmp = L*Mnorm;
		Matrix<double> Umat = mean+tmp;

		return Umat;
} 
