#include "integrator.hpp"


using namespace std;


integrator::integrator(vector<double>* Y0, void (*dY_func)(vector<double>*,vector<double>*),
		       array<double,2> t_span, double h) {
  this-> Y0 = Y0;
  this-> dY_func = dY_func;
  this-> h = h;
  this-> t_span = t_span;
  numOfYs = Y0->size();
  for (int i=0;i<numOfYs;++i){k1.push_back(0);}
  dt = t_span[1]-t_span[0];
  steps = (int)(ceil(dt/h));
}

integrator::~integrator(){
  cout << "integrator destructor" << endl;
}

void integrator::solve(){
  for (int i=0;i<getSteps();++i){
    dY_func(Y0,&k1);
    for (int j=0;j<getNumOfYs();++j){
      Y0->data()[j] += h* (k1[j]);
    }
  }
}

rk2::rk2(vector<double>* Y0, void (*dY_func)(vector<double>*,vector<double>*),
	   array<double,2> t_span, double h): integrator(Y0,dY_func,t_span,h){
  for (int i=0;i<getNumOfYs();++i){k2.push_back(0);}
  for (int i=0;i<getNumOfYs();++i){y_half_h.push_back(0);}
}

rk2::~rk2(){
  cout << "rk2 destructor" << endl;
}




void rk2::solve(){
  for (int i=0;i<getSteps();++i){
    dY_func(Y0,getK1()); // Calculates k1
    for (int j=0;j<getNumOfYs();++j){
      y_half_h[j] = Y0->data()[j] + (h*0.5)* (getK1()->data()[j]);
    }
    dY_func(&y_half_h,&k2); // Calculates k2 from half step
    for (int j=0;j<getNumOfYs();++j){
      Y0->data()[j] += h* (k2[j]);
    }    
  }
}

double Y0_exact(double t) {
  double y = 3.0 * exp(-2.0 * t);
  return y;
}


void dY0_func(vector<double>* y, vector<double>* dy) {
  for (int i =0;i<(int)y->size();++i){
    dy->data()[i] = -2.0 * y->data()[i];
  }
}

