// header guard // 
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <array>
#include <string>
#include <cmath>

using namespace std;

class integrator
{
public:
  vector<double>* Y0;
  double h;
  array<double,2> t_span;
  void (*dY_func)(vector<double>* y,vector<double>* dy);
  void solve();
  string getName(){return name;}
  int getSteps(){return steps;}
  int getNumOfYs(){return numOfYs;}
  vector<double>* getK1(){return &k1;}
  integrator(vector<double>*, void (*)(vector<double>*,vector<double>*),array<double,2>, double);
  ~integrator();
private:
  double dt;
  int steps;
  string name = "Euler";
  int numOfYs;
  vector<double> k1;
};



class rk2 : public integrator {
public:
  rk2(vector<double>*, void (*)(vector<double>*,vector<double>*),array<double,2>, double);
  ~rk2();
  void solve();
private:
  string name= "RK2";
  vector<double> y_half_h;
  vector<double> k2;
};



#endif
