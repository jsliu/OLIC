#include "rand.hpp"
#include "BCF.hpp"
   
BCF::BCF(double *observ, double *variates, int size, int order, 
	 double *del_var, double gamma, double nu, double tran_prob)
    :y(observ),x(variates),sz(size),ode(order),del(del_var),ga(gamma),vu(nu),tp(tran_prob)
{
  // store differnt list at differnt time
  storage = new list<Particle>[sz];
  
  // a matrix to store the expected mean and variance
  // and extra parameters for sigma
  extra_sigma = new double*[sz]; 
  expect_mean = new double[sz];
  delta_for_beta0 = new double[sz];
  
  for(int i=0;i<sz;i++)
    {
      extra_sigma[i] = new double[2];
      for(int j=0;j<2;j++)
	extra_sigma[i][j] = 0.0;
      
      expect_mean[i] = 0.0;
      delta_for_beta0[i] = 0.0;
    }

  // file to store the sigma at each time for MixBCF
  File_sigma_time.open("../output/sigma_time",ios::out);

  // obtain the seed
  getRNG();
}


// don't miss the virtual destructor
BCF::~BCF()
{	
  //free
  for(int i=0;i<sz;i++)
    delete []extra_sigma[i];
  delete []extra_sigma;
  
  delete []expect_mean;
  delete []delta_for_beta0;
  
  delete []storage;

  // close the file
  File_sigma_time.close();

  // save the seed
  putRNG();
}


void BCF::add_to_tail(int time)
{
  // upgrade the existed particle
  for(LP_iter = LParticle.begin();LP_iter!=LParticle.end();LP_iter++)
    calc_weight(LP_iter->pos,time);

  // add new nodes
  for(int m=1;m<=qmod;m++)
    {
      tmp_part.mod = m;
      calc_weight(time-1,time);
      LParticle.push_back(tmp_part);
    }          
  
  // normalising the weights
  normalize();
}


void BCF::copy_list(int time)
{
  storage[time-1] = LParticle;
  
  if(time % 100==0)
    cout<<time<<" data has been analysed with Particles: "<<storage[time-1].size()<<"\n";
}


void BCF::remove_from_list()
{
  LP_iter = LParticle.begin();
  while(LP_iter->wei!=0.0)
    LP_iter++;
  
  if(LP_iter->wei==0.0)
    LParticle.erase(LP_iter);
  
  else  
    cout<<"No node can be deleted."<<endl;		
}


bool BCF::search()
{
    LP_iter = LParticle.begin();
    while(LP_iter!=LParticle.end())
        if(LP_iter->wei==0.0)
            return true;
        else 
            LP_iter++;
		return false;
}


bool BCF::is_bigger_minsum(double hold)
{
    double sum_min = 0.0;
    for(LP_iter=LParticle.begin();LP_iter!=LParticle.end();LP_iter++)
        sum_min += min(1.0,LP_iter->wei/hold);	
    sum_min = floor(sum_min+1);
    
    return LParticle.size() > sum_min ? true:false;
}


void BCF::src(double hold)
{ 
    // a kind of stratified sampling
    double U = runif();
    for(LP_iter=LParticle.begin();LP_iter!=LParticle.end();LP_iter++)
        if(LP_iter->wei<hold)
        {
            U -= LP_iter->wei/hold;
            if(U<0)
            {
                LP_iter->wei = hold;
                U++;
            }
            else 
                LP_iter->wei = 0.0;
        }
}


void BCF::sor(int m, int n)
/**********************************************
 *
 *  produce n particles from m particles (n<m)
 *
 **********************************************/
{
    // calcualte the threshold so that "c satisfy sum_{i=1}^{m}min(1,cq_i)=n"
    double c = n;
    double sum_prob, p;
    long l, pre_l;            // number of particles whose probabilities are greater than 1/c
    do{
        sum_prob = 0.0;
        pre_l = l;
        l = 0;
        p=1/c;
        for(LP_iter=LParticle.begin();LP_iter!=LParticle.end();LP_iter++)
            if(LP_iter->wei<p)
            {
                sum_prob += LP_iter->wei;
                l++;
            }
        c = (n+l-m)/sum_prob;
    }
    while(l!=pre_l && l!=n-m);
    
    double thresh = 1/c;
    src(thresh);
}


void BCF::normalize()
{
    double sum_prob = 0.0;
    for(LP_iter = LParticle.begin();LP_iter!=LParticle.end();LP_iter++)
        sum_prob += LP_iter->wei;
 
    for(LP_iter = LParticle.begin();LP_iter!=LParticle.end();LP_iter++)
        LP_iter->wei /= sum_prob;
}


int BCF::list_length()
{
    return LParticle.size();
}


// simulation of change points and model choices at those points
// this is a backward simulation
void BCF::simu_bcf(int sample_size)
{        
  // output
  fstream File_cpt, File_mod, File_fit, File_appr;
  File_cpt.open("../output/ChangePoints",ios::out);
  File_mod.open("../output/ModelAtChpts",ios::out);
  File_appr.open("../output/CptPosAppr",ios::out);

  if(File_cpt.fail() || File_mod.fail() || File_appr.fail())
    {
      cout<<"Open one of files failed"<<endl;
      exit(1);
    }

  // no type of changepoints
  int *chpt_pos = new int[sample_size];
  int *model_cho = new int[sample_size];
  int *num_chpt = new int[sample_size];
  
  // dynamic 2-d array
  // store change points and model choices
  vector<int> *chpt = new vector<int>[sample_size];   
  vector<int> *modl = new vector<int>[sample_size];
  
  // store the fitted data
  //double **fit_data = new double*[sample_size];
  //for(int i=0;i<sample_size;i++)
  //  fit_data[i]= new double[sz];
  
  // one row of basis matrix 
  //double *basis = new double[ode];
  
  int last_chpt = sz;
  for(int i=0;i<sample_size;i++)
    {
      chpt_pos[i] = sz;
      num_chpt[i] = 0;
    }
  
  while(last_chpt>0)
    {
      int nt = 0;
      for(int i=0;i<sample_size;i++)
	if(chpt_pos[i]==last_chpt)
	  nt++;
      
      int len = storage[last_chpt-1].size();
      double *weights = new double[len];
      vector<Particle> vParticle(len);
      int index_cpt, index_mod;
      
      int count = 0;
      for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
	{
	  weights[count] = LP_iter->wei;
	  vParticle[count] = *LP_iter;
	  count++;
	}
      
      // only simulate when there are different change points
      if(nt>0)
	{
	  // sampling change point and model choices nt times
	  double U[1];
	  for(int i=0;i<sample_size;i++)
	    if(chpt_pos[i]==last_chpt)
	      {
		U[0] = runif();
		
		// sampling change points and model choices
		rdunif(weights,1.0,&index_cpt,U,1,len);
		chpt_pos[i] = vParticle[index_cpt].pos;
		model_cho[i] = vParticle[index_cpt].mod;
		
		// recorde them
		if(chpt_pos[i]<sz)
		  {
		    chpt[i].push_back(chpt_pos[i]);
		    modl[i].push_back(model_cho[i]);
		  }
	      }
	} // nt
      
      File_appr<<"NA"<<"\t"<<last_chpt<<"\t"<<last_chpt<<"\n";
      for(int i=0;i<len;i++)
	File_appr<<vParticle[i].wei<<"\t"<<vParticle[i].pos<<"\t"<<vParticle[i].mod<<"\n";
      
      // free
      delete []weights;
      last_chpt--;
    } //while
  
      // write down the parameters and change points
  for(int i=0;i<sample_size;i++)
    {
      for(int j=0;j<chpt[i].size();j++)
	{
	  File_cpt<<chpt[i][j]<<"\t";
	  File_mod<<modl[i][j]<<"\t";
	}
      File_cpt<<endl;
      File_mod<<endl;
    }
  
  // free
  delete []chpt_pos;
  delete []model_cho;
  delete []num_chpt;
  
  delete []chpt;
  delete []modl;	


  File_cpt.close();
  File_mod.close(); 
  File_appr.close();
}


/*********************************************************************************
 *
 *                    Derived Class DistBCF
 *
 *********************************************************************************/


DistBCF::DistBCF(double *observ, double *variates, int size, int order, 
		 double *del_var, double gamma, double nu, double tran_prob, double *prior)
  :BCF(observ, variates, size, order, del_var, gamma, nu, tran_prob)
{  
  assert(ode>0 && ode<=5);
  
  // initialize the prior probability of models
  qmod = ode;		
  pmod = new double[qmod];
  for(int i=0;i<qmod;i++)
    pmod[i] = prior[i];
}


DistBCF::~DistBCF()
{	
  delete []pmod;
}


void DistBCF::calc_weight(int cpt_time, int cur_time)
{		
  double log_lik;
  double alpha_sigma,beta_sigma;
  int model;
  int interval = cur_time-cpt_time;

  if(interval==1)
    {
      model = tmp_part.mod;
      Matrix<double> mu_beta(model,1),Sigma_beta(model,model);
      tmp_part.pos = cpt_time;

      for(int i=0;i<model;i++)
	Sigma_beta[i][i] = del[i];

      alpha_sigma = vu;
      beta_sigma = ga;
 
      para_update(mu_beta,Sigma_beta,alpha_sigma,beta_sigma,log_lik,model,cpt_time,cur_time);

      for(int i=0;i<model;i++)
	{
	  tmp_part.mean_beta[i]=mu_beta[i][0];
	  for(int j=0;j<model;j++)
	    tmp_part.var_beta[i][j]=Sigma_beta[i][j];
	}

      tmp_part.nu_sigma = alpha_sigma;
      tmp_part.gamma_sigma = beta_sigma;
      tmp_part.wei = exp(log_lik)*tp*pmod[model-1];	

    }else{
      model = LP_iter->mod;
      Matrix<double> mu_beta(model,1),Sigma_beta(model,model);    
  
      for(int i=0;i<model;i++)
	{
	  mu_beta[i][0]= LP_iter->mean_beta[i];
	  for(int j=0;j<model;j++)
	    Sigma_beta[i][j] = LP_iter->var_beta[i][j];
	}

      alpha_sigma = LP_iter->nu_sigma;
      beta_sigma = LP_iter->gamma_sigma;

      para_update(mu_beta,Sigma_beta,alpha_sigma,beta_sigma,log_lik,model,cpt_time,cur_time);

      for(int i=0;i<model;i++)
	{
	  LP_iter->mean_beta[i]=mu_beta[i][0];
	  for(int j=0;j<model;j++)
	    LP_iter->var_beta[i][j]=Sigma_beta[i][j];
	}

      LP_iter->nu_sigma = alpha_sigma;
      LP_iter->gamma_sigma = beta_sigma; 
      LP_iter->wei = exp(log_lik)*(1-tp)*LP_iter->wei;
    }
}

void DistBCF::para_update(Matrix<double>& mu_beta, Matrix<double>& Sigma_beta, 
			  double& alpha_sigma, double& beta_sigma, double& log_lik,
			  int model, int cpt_time, int cur_time)
{
  const double pi = 3.141592653589793238462643;
  double y_current = y[cur_time-1];

  // set identity matrix
  Matrix<double> identity(model,model);
  for(int i=0;i<model;i++)
    identity[i][i] = 1.0;   
  
  // basis matrix, G
  Matrix<double> basis_seg(1,model);
  for(int j=0;j<model;j++)  
    basis_seg[0][j] = pow(x[cur_time-1]-x[cpt_time],j);
  
  // Firstly, inverse of the Sigma
  Cholesky<double> cho_Sigma(Sigma_beta);
  Matrix<double> inv_Sigma = cho_Sigma.solve(identity);
  double det_Sigma = cho_Sigma.det();
  
  // Secondly, G^t*G
  Matrix<double> GG = transpose_mult<double>(basis_seg,basis_seg);
  
  // Thirdly, calculate G^t*y
  Matrix<double> Gy(model,1);
  for(int j=0;j<model;j++)
    Gy[j][0] = basis_seg[0][j]*y_current;
  
  // M = (G^t*G+invSigma)^(-1)
  Matrix<double> Mtmp = GG+inv_Sigma;
  Cholesky<double> cho_Mtmp(Mtmp);
  Matrix<double> M = cho_Mtmp.solve(identity);
  Cholesky<double> cho_M(M);
  double det_M = cho_M.det();
  
  // N = G^t*y+invSigma*mu
  Matrix<double> N = Gy+inv_Sigma*mu_beta;
  
  // mu^t*invSigma*mu and N^t*M*N
  Matrix<double> tmp1 = transpose_mult<double>(mu_beta,inv_Sigma);
  Matrix<double> part1 = tmp1*mu_beta;
  Matrix<double> tmp2 = transpose_mult<double>(N,M);
  Matrix<double> part2 = tmp2*N;
  double y_norm2 = y_current*y_current+part1[0][0]-part2[0][0];		
  
  // calculation of the log likelihood here
  double lgMD = log(det_M)-log(det_Sigma);
  double lgexp1 = (alpha_sigma+1)*log(beta_sigma+y_norm2);
  double lgexp2 = alpha_sigma*log(beta_sigma); 
  double lgamma1 = gammln_p((alpha_sigma+1)/2.0);
  double lgamma2 = gammln_p(alpha_sigma/2.0);
  log_lik = (-log(pi)+lgMD-lgexp1+lgexp2)/2.0+lgamma1-lgamma2;
  
  // update parameters of beta		
  mu_beta = M*N;
  Sigma_beta = M;
  
  // update parameters of sigma2
  alpha_sigma = alpha_sigma+1;       
  beta_sigma = beta_sigma+y_norm2;
}

void DistBCF::expect_mean_var(int current){}

/*********************************************************************************
 *
 *                     Derived Classs MixBCF
 *
 *********************************************************************************/


MixBCF::MixBCF(double *observ, double *variates, int size, int order, 
	       double *del_var, double gamma, double nu, double tran_prob, double *prior)
  :BCF(observ,variates,size,order,del_var,gamma,nu,tran_prob)
{
  // assign prior for modelc choices
  qmod = 2;
  pmod = new double[qmod];
  for(int i=0;i<qmod;i++)
    pmod[i] = prior[i];          
}


MixBCF::~MixBCF()
{
  delete []pmod;
}

// The likelihood is the weight
void MixBCF::calc_weight(int cpt_time, int cur_time)
{
  const double pi = 3.141592653589793238462643;

  // set the identity matrix
  Matrix<double> identity(ode,ode);
  for(int i=0;i<ode;i++)
    identity[i][i] = 1.0;   
  
  Matrix<double> mu_beta(ode,1),Sigma_beta(ode,ode);
  double alpha_sigma,beta_sigma;
  double log_lik;                         
  double len_sigma = 0.0;
  double norm2_sigma= 0.0;
  int model;
  double y_current = y[cur_time-1];
  
  int interval = cur_time-cpt_time;
  
  if(interval==1)
    {
      
      model = tmp_part.mod;
      
      // assign spaces for posterior mean and variance of beta
      for(int i=0;i<ode;i++)
	{
	  tmp_part.mean_beta[i] = 0.0;
	  for(int j=0;j<ode;j++)
	    tmp_part.var_beta[i][j] = 0.0;
	}
      tmp_part.pos = cpt_time;
      
      for(int i=0;i<ode;i++)
	//(*tmp_part.var_beta)[i][i]= del[i];
	tmp_part.var_beta[i][i]= del[i];    

      for(int i=0;i<ode;i++)
	{
	  mu_beta[i][0] = tmp_part.mean_beta[i];
	  for(int j=0;j<ode;j++)
	    Sigma_beta[i][j] = tmp_part.var_beta[i][j];
	}

      // extra parameters for beta_0 and sigma2
      if(model==2 && cpt_time!=0)
	{
	  mu_beta[0][0] = expect_mean[cpt_time-1];
	  Sigma_beta[0][0] = delta_for_beta0[cpt_time-1];
	}
      
      if(cpt_time!=0)
	{
	  double m1=extra_sigma[cpt_time-1][0];
	  double m2=extra_sigma[cpt_time-1][1];
	 
	  len_sigma= 2.0*m1*m1/(m2-m1*m1);
	  norm2_sigma= 2.0*m1/(m2-m1*m1); //IS FACTOR OF 2 CORRECT i.e. gamma(vu/2,ga/2)?
	  
	  //CHANGE -- OLD CODE
	  /*
	  len_sigma = extra_sigma[cpt_time-1][0];
	  norm2_sigma = extra_sigma[cpt_time-1][1];
	  */
	}
      
      alpha_sigma = vu+len_sigma;
      beta_sigma = ga+norm2_sigma;
    }else{

      for(int i=0;i<ode;i++)
	{
	  mu_beta[i][0]= (LP_iter->mean_beta)[i];
	  for(int j=0;j<ode;j++)
	    Sigma_beta[i][j] = (LP_iter->var_beta)[i][j];
	}

      alpha_sigma = LP_iter->nu_sigma;
      beta_sigma = LP_iter->gamma_sigma;
    }
  
  // basis matrix, G
  Matrix<double> basis_seg(1,ode);
  for(int j=0;j<ode;j++)  
    basis_seg[0][j] = pow(x[cur_time-1]-x[cpt_time],j);
  
  // Firstly, inverse of the Sigma
  Cholesky<double> cho_Sigma(Sigma_beta);
  Matrix<double> inv_Sigma = cho_Sigma.solve(identity);
  double det_Sigma = cho_Sigma.det();
  
  // Secondly, G^t*G
  Matrix<double> GG = transpose_mult<double>(basis_seg,basis_seg);
  
  // Thirdly, calculate G^t*y
  Matrix<double> Gy(ode,1);
  for(int j=0;j<ode;j++)
    Gy[j][0] = basis_seg[0][j]*y_current;
  
  // M = (G^t*G+invSigma)^(-1)
  Matrix<double> Mtmp = GG+inv_Sigma;
  Cholesky<double> cho_Mtmp(Mtmp);
  Matrix<double> M = cho_Mtmp.solve(identity);
  Cholesky<double> cho_M(M);
  double det_M = cho_M.det();
  
  // N = G^t*y+invSigma*mu
  Matrix<double> N = Gy+inv_Sigma*mu_beta;
  
  // mu^t*invSigma*mu and N^t*M*N
  Matrix<double> tmp1 = transpose_mult<double>(mu_beta,inv_Sigma);
  Matrix<double> part1 = tmp1*mu_beta;
  Matrix<double> tmp2 = transpose_mult<double>(N,M);
  Matrix<double> part2 = tmp2*N;
  double y_norm2 = y_current*y_current+part1[0][0]-part2[0][0];		
  
  // calculation of the log likelihood here
  double lgMD = log(det_M)-log(det_Sigma);
  double lgexp1 = (alpha_sigma+1)*log(beta_sigma+y_norm2);
  double lgexp2 = alpha_sigma*log(beta_sigma); 
  double lgamma1 = gammln_p((alpha_sigma+1)/2.0);
  double lgamma2 = gammln_p(alpha_sigma/2.0);
  log_lik = (-log(pi)+lgMD-lgexp1+lgexp2)/2.0+lgamma1-lgamma2;
		
  // update parameters of beta		
  mu_beta = M*N;
  Sigma_beta = M;
  
  // update parameters of sigma2
  alpha_sigma = alpha_sigma+1;       
  beta_sigma = beta_sigma+y_norm2;
  
  // calculate E(z_t|y_1:t-1), where z_t is the signal
  // borrow information from previous segment, so there is no need for the discontinuous model

  // new basis matrix, as time goes to t from t-1
  if(cur_time<sz)
    for(int j=0;j<ode;j++)  
      basis_seg[0][j] = pow(x[cur_time]-x[cpt_time],j);

  // intercept = (beta0, beta1, beta2)^T*(1, x, x^2)  	
  Matrix<double> expect_mu = basis_seg*mu_beta;
  Matrix<double> expect_var = basis_seg*transpose_mult2(Sigma_beta,basis_seg); 

  
  // update and record
  if(interval==1)
    {
      tmp_part.meanbeta02 = expect_mu[0][0];
      tmp_part.delta02 = expect_var[0][0];
      //tmp_part.extra_nu = alpha_sigma-vu;
      //tmp_part.extra_gamma = beta_sigma-ga;
      
      for(int i=0;i<ode;i++)
	{
	  tmp_part.mean_beta[i]=mu_beta[i][0];
	  for(int j=0;j<ode;j++)
	    tmp_part.var_beta[i][j]=Sigma_beta[i][j];
	}

      tmp_part.nu_sigma = alpha_sigma;
      tmp_part.gamma_sigma = beta_sigma;
      //tmp_part.inv_sigma = rgamma(alpha_sigma/2.0)/(beta_sigma/2.0);
      tmp_part.wei = exp(log_lik)*tp*pmod[model-1];    
    }else{
      LP_iter->meanbeta02 = expect_mu[0][0];
      LP_iter->delta02 = expect_var[0][0];
      //LP_iter->extra_nu = alpha_sigma-vu;
      //LP_iter->extra_gamma = beta_sigma-ga;
      
      for(int i=0;i<ode;i++)
	{
	  LP_iter->mean_beta[i]=mu_beta[i][0];
	  for(int j=0;j<ode;j++)
	    LP_iter->var_beta[i][j]=Sigma_beta[i][j];
	}

      LP_iter->nu_sigma = alpha_sigma;
      LP_iter->gamma_sigma = beta_sigma; 
      //LP_iter->inv_sigma = rgamma(alpha_sigma/2.0)/(beta_sigma/2.0);
      LP_iter->wei = exp(log_lik)*(1-tp)*LP_iter->wei;
    }	
}


void MixBCF::expect_mean_var(int current)
{
  double nu2 = 0.0, gamma2=0.0;
  double inv_sigma=0.0;  

    for(LP_iter=LParticle.begin();LP_iter!=LParticle.end();LP_iter++)
    {
      // extra parameters for sigma2
      //extra_sigma[current-1][0] += LP_iter->extra_nu*LP_iter->wei/LP_iter->extra_gamma; //mean

      //extra_sigma[current-1][1] += (2.0+LP_iter->extra_nu)*LP_iter->extra_nu*LP_iter->wei/(LP_iter->extra_gamma*LP_iter->extra_gamma);                                                      //mean square -- is factor of 2 correct (see comment in MixBCF)
      //mean
      extra_sigma[current-1][0] += (LP_iter->nu_sigma-vu)*LP_iter->wei/(LP_iter->gamma_sigma-ga);
      
      //mean square -- is factor of 2 correct (see comment in MixBCF)
      extra_sigma[current-1][1] += (2.0+LP_iter->nu_sigma-vu)*(LP_iter->nu_sigma-vu)*LP_iter->wei/((LP_iter->gamma_sigma-ga)*(LP_iter->gamma_sigma-ga));                                          
      

      //OLD CODE 
      /*
      extra_sigma[current-1][0] += LP_iter->nu_sigma*LP_iter->wei;
      extra_sigma[current-1][1] += LP_iter->gamma_sigma*LP_iter->wei;
      */

      // expect mean and delta02 for the continuous model
      expect_mean[current-1] += LP_iter->meanbeta02*LP_iter->wei;
      delta_for_beta0[current-1] += LP_iter->delta02*LP_iter->wei; 

    }	
    double m1 = extra_sigma[current-1][0];
    double m2 = extra_sigma[current-1][1];
    double len_s= 2.0*m1*m1/(m2-m1*m1);
    double norm_s = 2.0*m1/(m2-m1*m1);

    /*
    for(int i=0;i<100;i++)
      {
	  inv_sigma = rgamma(len_s/2)/(norm_s/2);
	  File_sigma_time<<1/inv_sigma<<"\t";
      }
    File_sigma_time<<"\n";
    */

    File_sigma_time<<len_s<<"\t"<<norm_s<<"\n";
}


void MixBCF::simu_bcf(int sample_size)
{
  const double pi = 3.141592653589793238462643;

  // output
  fstream File_cpt, File_mod, File_fit, File_sigma, File_appr;
  File_cpt.open("../output/ChangePoints",ios::out);
  File_mod.open("../output/ModelAtChpts",ios::out);
  File_sigma.open("../output/Variance",ios::out);
  File_appr.open("../output/CptPosAppr",ios::out);

  if(File_cpt.fail() || File_mod.fail() || File_sigma.fail() || File_appr.fail())
    {
      cout<<"Open one of files failed"<<endl;
      exit(1);
    }
  // considering the types of change points  
  // sampling change point and model choices nt times
  for(int sp=0;sp<sample_size;sp++)
    {
      int last_chpt = sz;
      int last_mod = 1;
      int index_cpt;
      double inv_sigma2, beta0,deltax;
      Matrix<double> bs(1,ode); 
      
      while(last_chpt>0)
	{
	  int len = storage[last_chpt-1].size();
	  vector<Particle> vParticle(len);
	  double *weights = new double[len];
	  double *logwei = new double[len];
	  
	  double U[1];
	  U[0] = runif();
	  double sum_wei = 0.0;
	  int count = 0;
	  
	  /******************* previous code ***************************

	           if(last_chpt==sz){
		   for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		   {
		   weights[count] = LP_iter->wei;
		   vParticle[count] = *LP_iter;
		   
		   sum_wei += weights[count];
		   count++;
		   }
		   rdunif(weights,sum_wei,&index_cpt,U,1,len);
		   
		   // simulate 1/sigma2 and beta_0
		   inv_sigma2 = rgamma(vParticle[index_cpt].nu_sigma)/vParticle[index_cpt].gamma_sigma;
		   Matrix<double> mean_tmp = *vParticle[index_cpt].mean_beta;
		   Matrix<double> var_tmp = mult(*vParticle[index_cpt].var_beta,1/inv_sigma2);
		   Matrix<double> beta = mult_norm(mean_tmp,var_tmp);
		   
		   beta0 = beta(1,1);
		   }else{
		   for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		   {
		   double invsigma_loglik = LP_iter->nu_sigma*log(LP_iter->gamma_sigma)-gammln_p(LP_iter->nu_sigma)+(LP_iter->nu_sigma-1)*log(inv_sigma2)-LP_iter->gamma_sigma*inv_sigma2;
		   
		   // discontinuous case
		   if(last_mod==1)
		   {
		   Matrix<double> mean_tmp = *LP_iter->mean_beta;
		   Matrix<double> var_tmp = mult(*LP_iter->var_beta,1/inv_sigma2);
		   Matrix<double> beta = mult_norm(mean_tmp,var_tmp);                
		   beta0 = beta(1,1);
		   
		   logwei[count] = log(LP_iter->wei)+invsigma_loglik;
		   }else{ // continuous case
		   Matrix<double> mean_tmp = *LP_iter->mean_beta;
		   Matrix<double> var_tmp = mult(*LP_iter->var_beta,1/inv_sigma2);
		   double mean0 = mean_tmp(1,1);
		   double var0 = var_tmp(1,1);  
		   double beta0_loglik = -log(sqrt(2*pi*var0))-pow(beta0-mean0,2.0)/(2*var0);
		   logwei[count] = log(LP_iter->wei)+invsigma_loglik+beta0_loglik;
		   }
		   vParticle[count] = *LP_iter;
		   count++;
		   }
		   
		   // normalise
		   count = 0;
		   double max_logwei = logwei[0];
		   for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		   {
		   if(logwei[count]>max_logwei)
		   max_logwei = logwei[count];
		   count++;
		   }
		   
		   count = 0;
		   for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		   {
		   weights[count] = exp(logwei[count]-max_logwei);
		   sum_wei += weights[count];
		   count++;
		   }
		   rdunif(weights,sum_wei,&index_cpt,U,1,len);
		   } // if-else

	  *******************************************************************/
	  //**
	  // sampling change points and model choices with prob P(C_T|y_{1:T})
	  if(last_chpt==sz)  
	    {            
	      for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		{
		  weights[count] = LP_iter->wei;
		  vParticle[count] = *LP_iter;
		  
		  sum_wei += weights[count];
		  count++;
		}
	      rdunif(weights,sum_wei,&index_cpt,U,1,len);
	      
	      // simulate 1/sigma2 and beta_0
	      double nu_tmp = vParticle[index_cpt].nu_sigma/2.0;
	      double gamma_tmp = vParticle[index_cpt].gamma_sigma/2.0;
	      inv_sigma2 = rgamma(nu_tmp)/gamma_tmp;
	      File_sigma<<1/inv_sigma2<<"\n";
	      
	      if(vParticle[index_cpt].mod==2)
		{ 
		  //only for continuous changepoint
		  Matrix<double> mean_tmp(ode,1);
		  Matrix<double> var_tmp(ode,ode);
		  for(int i=0;i<ode;i++)
		    {
		      mean_tmp[i][0]= vParticle[index_cpt].mean_beta[i];
		      for(int j=0;j<ode;j++)
			var_tmp[i][j]= vParticle[index_cpt].var_beta[i][j]/inv_sigma2;
		    }
		  
		  Matrix<double> beta = mult_norm(mean_tmp,var_tmp); 
		  beta0 = beta(1,1);
		}
	    }else{
	      for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		{
		  double invsigma_loglik = LP_iter->nu_sigma*log(LP_iter->gamma_sigma)-gammln_p(LP_iter->nu_sigma)+(LP_iter->nu_sigma-1)*log(inv_sigma2)-LP_iter->gamma_sigma*inv_sigma2;
		  
		  //contribution for both cases
		  logwei[count] = log(LP_iter->wei)+invsigma_loglik;
		  
		  //continuous case
		  if(last_mod==2){
		    //Matrix<double> mean_tmp = *LP_iter->mean_beta;
		    //Matrix<double> var_tmp = mult(*LP_iter->var_beta,1/inv_sigma2);
		    
		    // change in x across segment
		    //deltax = x[last_chpt]-x[LP_iter->pos];
		    
		    // current basis matrix of 1 X 3
		    //for(int k=0;k<ode;k++)
		    // bs[0][k] = pow(deltax,k);
		    
		    //double mean0 = (bs*mean_tmp)(1,1);
		    //double var0 = (bs*transpose_mult2(var_tmp,bs))(1,1);
		    double mean0 = LP_iter->meanbeta02;
		    double var0 = LP_iter->delta02/inv_sigma2;
		    
		    double beta0_loglik = -log(sqrt(2.0*pi*var0))-pow(beta0-mean0,2.0)/(2.0*var0);
		    double tmp = log(LP_iter->wei);
		    logwei[count] = log(LP_iter->wei)+invsigma_loglik+beta0_loglik;
		  }
		  vParticle[count] = *LP_iter;
		  count++;
		}
	      
	      // normalise
	      count = 0;
	      double max_logwei = logwei[0];
	      for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		{
		  if(logwei[count]>max_logwei)
		    max_logwei = logwei[count];
		  count++;
		}
	      
	      count = 0;
	      for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		{
		  weights[count] = exp(logwei[count]-max_logwei);
		  sum_wei += weights[count];
		  count++;
		}
	      
	      count = 0;
	      for(LP_iter=storage[last_chpt-1].begin();LP_iter!=storage[last_chpt-1].end();LP_iter++)
		{
		  weights[count] /= sum_wei;
		  count++;
		}
     
	      rdunif(weights,1.0,&index_cpt,U,1,len);
	      
	      //now simulate new beta0
	      if(vParticle[index_cpt].mod==2){
		if(last_mod==1){
		  //discontinuous case -- simulate from distribution stored for the 
		  // particle, given value of inv_sigma
		  Matrix<double> mean_tmp(ode,1);
		  Matrix<double> var_tmp(ode,ode);		
		  for(int i=0;i<ode;i++)
		    {
		      mean_tmp[i][0]= vParticle[index_cpt].mean_beta[i];
		      for(int j=0;j<ode;j++)
			var_tmp[i][j]= vParticle[index_cpt].var_beta[i][j]/inv_sigma2;
		    }
		  
		  Matrix<double> beta = mult_norm(mean_tmp,var_tmp);
		  beta0 = beta(1,1);		 
		}
		
		if(last_mod==2){
		  //As above but need to condition on beta0 for the previous segment.
		  //observation is that beta0=b0+dx*b1+dx^2*b2
		  //if x sim N(mu,Q) and Ax=b
		  Matrix<double> mu(ode,1);
		  Matrix<double> var(ode,ode);
		  for(int i=0;i<ode;i++)
		    {
		      mu[i][0]= vParticle[index_cpt].mean_beta[i];
		      for(int j=0;j<ode;j++)
			var[i][j]= vParticle[index_cpt].var_beta[i][j]/inv_sigma2;
		    }
		  
		  // conditional mean and variance
		  deltax = x[last_chpt]-x[vParticle[index_cpt].pos];
		  for(int k=0;k<ode;k++)
		    bs[0][k] = pow(deltax,k);
		  Matrix<double> tmp = transpose_mult2(var,bs); // QA^T
		  double cent = (bs*tmp)(1,1); // AQA^T, but here is a 1x1 matrix
		  Matrix<double> tmp2 = mult(tmp,1.0/cent); // QA^T(AQA^T)^{-1}
		  
		  //mu* = mu-QA^T(AQA^T)^{-1}(Amu-b)
		  Matrix<double> tmp3 = mult(tmp2,(bs*mu)(1,1)-beta0);
		  Matrix<double> mean_tmp = mu-tmp3;
		  
		  //Q* = Q-QA^T(AQA^T)^{-1}AQ
		  Matrix<double> var_tmp = var-transpose_mult2(tmp2,tmp);
		  
		  //now simulate
		  Matrix<double> beta = mult_norm(mean_tmp,var_tmp);
		  
		  beta0 = beta(1,1);
		}
	      } // if mod==2
	    } // if-else
	  
	  // store the values
	  last_chpt = vParticle[index_cpt].pos;
	  last_mod = vParticle[index_cpt].mod;
	  
	  File_cpt<<last_chpt<<"\t";
	  File_mod<<last_mod<<"\t";
	  File_appr<<weights[index_cpt]<<"\t";

	  // free
	  delete []weights;
	  delete []logwei;
	} //while
      File_cpt<<endl;
      File_mod<<endl;
      File_appr<<endl;
    } // for sp

  // close files
  File_cpt.close();
  File_mod.close(); 
  File_sigma.close();
  File_appr.close();
}
