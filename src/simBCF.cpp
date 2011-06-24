#include "rand.hpp"
#include "BCF.hpp"

// forward declaration
void recur_calc(double*, double*, int, int, int, double, double, double, double, double, bool, bool);

void recur_calc(double *obs,double *vat, int ssiz,int siz,int ord,double *dva,
		double gam, double nus,double trp,double the,double *pm,bool con, bool resamp)
/**************************************************************************************
 *
 *  obs=observations, vat=variates, ssize=sampling size, siz=size, ord=order, 
 *  trp=transition probability, gam=gamma, nus=nu, dva=delta matrix, 
 *  pmo=prior prob of model choices
 *
 **************************************************************************************/
{
  BCF *bcf;
 
  if(con==true)
    bcf = new MixBCF(obs,vat,siz,ord,dva,gam,nus,trp,pm);   
  else
    bcf = new DistBCF(obs,vat,siz,ord,dva,gam,nus,trp,pm);
 
  for(int ctime=1;ctime<=siz;ctime++)
    {
      bcf->add_to_tail(ctime);
      
      // doing resampling 
      if(resamp == true)
	{
	  if(bcf->is_bigger_minsum(the))
	    {
	      bcf->src(the);
	      while(bcf->search())
		bcf->remove_from_list();
	      bcf->normalize();
	    }
	}
      
      // calculating the expected mean and variance for beta02
      // this is only used in the mixture model
      bcf->expect_mean_var(ctime);				
      
      // store the list 
      bcf->copy_list(ctime);
    }
  
  bcf->simu_bcf(ssiz);
  
  delete bcf;
}


/*********************************************************
int main(int argc, char* argv[])
{

  // read the data
  if(argc != 2)
    {
      cerr<<"usage: %s [i/p file]"<<argv[0]<<"\n";
      exit(1);
    }
  
  int sample_size,size;
  double gam,nu,th;
  char cont,resamp;

  // get user input information
  cout<<"\nOnline algorithm for Bayesian curve fitting\n";
  cout<<"============================================\n\n";
  cout<<"Size of the data:";
  cin>>size;
  cout<<"Sample size of simulations:";
  cin>>sample_size;
  cout<<"Hyper paramters for sigma2 (gamma):";
  cin>>gam;
  cout<<"Hyper paramters for sigma2 (nu):";
  cin>>nu;
  cout<<"Threshold for resampling:";
  cin>>th;
  cout<<"Allowing for continuous change point? (y/n):";
  cin>>cont;
  cout <<"Doing resampling? (y/n):";
  cin>>resamp;
  cout<<"\n\n";
  
  int order = 3;
  double dv[3] = {1e+4,1e+6,1e+8};          // remember to change these for different dataset
  double tp = 0.004;
  bool con,resam;
  double *pm;
  
  switch(cont)
    {
    case 'y' : 
      con = true;
      pm = new double[2];
      for(int i=0;i<2;i++)
	pm[i] = 0.5;
      break;
    case 'n' : 
      con = false;
      pm = new double[order];
      for(int i=0;i<order;i++)
	pm[i] = (double) 1.0/order;
      break;
    default:
      cout<<"please choose yes or no.";
      break;
    }
  
  switch(resamp)
    {
    case 'y' : 
      resam = true;
      break;
    case 'n' : 
      resam = false;
      break;
    default:
      cout<<"please choose yes or no.";
      break;
    }
  
  fstream File;
  File.open(argv[1],ios::in);
  
  if(File.fail()){
    cout<<"File can't be open\n";
    exit(1);
  }
  
  double *data = new double[size];
  double *coor = new double[size];
  
  int i = 0;
  while(!File.eof()){
    File >> coor[i] >> data[i];
    i++;
  }
  
  recur_calc(data,coor,sample_size,size,order,dv,gam,nu,tp,th,pm,con,resam);
  
  // free
  delete []data;
  delete []coor;
  
  File.close();
  
  return 0;
  
}
*********************************************************************/

//**
extern "C" {

void simBCF(double *data,double *coor,int *size,int *sample_size,int *order,double *dv,
	    double *gam, double *nu, double *tp,double *th,double *pm,bool *cont,bool *res)
{
  recur_calc(data,coor,*sample_size,*size,*order,dv,*gam,*nu,*tp,*th,pm,*cont,*res);
}
  
}
//**/

  
