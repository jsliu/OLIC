#ifndef _RAND_H
#define _RAND_H

#include <iostream>
#include "tnt/tnt_matrix.h"
#include "tnt/tnt_linalg.h"
using namespace TNT;
using namespace Linear_Algebra;


#define ITMAX 100 /*constants for calculating error function*/
#define EPS 3.0e-7
#define FPMIN 1.0e-30


/* file for random seed */
void seed_set(FILE *a);
void seed_put(FILE *a);
void getRNG(void);
void putRNG(void);

/* uniform random number generator */
double runif(void);
/* standard normal */
double stnorm(void);
/* gamma shape a */
double rgamma(double a);
/* discrete density - unnormalized weights w */
void rdunif(double *w,double sumw,int *index, double *U, int N,int M);
/* calculates the complementary error function (error 10^-7) */
double erfcc(double x);
/* calculate erf using series expansion */
double gammln_p(double xx);
void gcf(double *gammcf,double a, double x, double *gln);
void gser(double *gamser,double a, double x, double *gln);
double gammp(double a, double x);
double erffc(double x);
double rbeta(double a,double b);
int rbinom(int n, double p);
/* negative binomial distribution */
int rnbinom(int size, double pr);
/* multivariate normal distribution */
Matrix<double> mult_norm(const Matrix<double>&, const Matrix<double>&);

#endif // _RAND_H
