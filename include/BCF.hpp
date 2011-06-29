#ifndef _PARTICLE_LIST
#define _PARTICLE_LIST

#include <list>
#include<vector>
#include<iostream>
#include<fstream>

#include "mat_utsil.hpp"
#include "tnt_addon.hpp"

// define particle first
struct Particle
{
  int pos;                                          // position for change point
  int mod;                                          // type of change point
  double wei;                                       // posterior probability Pr(C_t|y_{1:t})
  double mean_beta[3], var_beta[3][3];              // posterior mean and variance of beta   
  double gamma_sigma,nu_sigma;                      // posterior parameters of sigma2
  double meanbeta02;                                // predict value at time t 
  double delta02;                                   // delta[0] in the continuous model
  //double extra_gamma,extra_nu;                    // extra parameters for continuous model
  double inv_sigma;                                 // inverse sigma square
};

class BCF
{
public:
  // defining functions here
  BCF(double*, double*, int, int, double*, double, double, double);
  virtual ~BCF();
  void add_to_tail(int);
  void copy_list(int);
  void remove_from_list();
  bool search();
  bool is_bigger_minsum(double);               // criterion for src
  void src(double);                            // stratified rejection control
  void sor(int, int);                          // stratified optimal resampling
  void normalize();                            // normalize the weights 
  int list_length();
  virtual void simu_bcf(int);                  // doing simulation of change points
  virtual void expect_mean_var(int) = 0;

protected:
  // calculate Pr(C_t|y_{1:t}) and block
  virtual void calc_weight(int, int) = 0;
  list<Particle> LParticle;                    // create a list of paritcles 
  list<Particle>::iterator LP_iter;            // the iterator
  list<Particle> *storage;                     // store all diffrents lists
  Particle tmp_part;                           // new node
  
  double *y,*x;                                // observations and variates
  int sz, ode;                                 // data size and model order
  double *del;                                 // deltas
  double ga,vu;                                // parameters for prior distribution
  double tp;                                   // transition probability
  double *expect_mean;                         // expected mean of beta02
  double *delta_for_beta0;                     // delta_0 for beta_0 in case of continuous model
  double **extra_sigma;                        // extra parameters for sigma, depending on data
  int qmod;                                    // number of model choices
  double *pmod;                                // prior probability of model choices
  fstream File_sigma_time;                     // file to sigma at each time
};


class DistBCF:public BCF
{
public:
  DistBCF(double*, double*, int, int, double*, double, double, double, double*);
  virtual ~DistBCF();
  
  // calculate Pr(C_t|y_{1:t}) and block, also return the logliklihood
  void calc_weight(int, int);
  
  // parameter update 
  void para_update(Matrix<double>&, Matrix<double>&, double&, double&, double&, int, int, int);

 // calculate the expected mean and variance
  void expect_mean_var(int);
};


class MixBCF:public BCF
{
public:
  MixBCF(double*, double*, int, int, double*, double, double, double, double*);
  virtual ~MixBCF();
  
  // calculate Pr(C_t|y_{1:t}) and block also return the logliklihood
  void calc_weight(int, int);                         

  // calculate the expected mean and variance
  void expect_mean_var(int);                         

  // doing simulation of change points from smoothing distribution
  void simu_bcf(int); 
};


#endif    // _PARTICLE_LIST
