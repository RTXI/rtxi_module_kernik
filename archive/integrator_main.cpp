#include <iostream>
#include "integrator.cpp"
#include <iomanip> // cout precision 
#include <chrono> // clocking



int main() {
  cout << setprecision(10);
  vector<double> y = {3};
  array<double,2> t = {0,2};
  double h = 0.0001;
  rk2 myintegrator(&y, (*dY0_func), t, h);

  // Solve and clock //
  auto t1 = chrono::high_resolution_clock::now();
  myintegrator.solve();
  auto t2 = chrono::high_resolution_clock::now();
  chrono::duration<double, std::milli> fp_ms = t2 - t1;


  // Calculate Error //
  double exact =Y0_exact(2.0);
  double err = Y0_exact(2.0)-y[0];
  double rel_err = abs(exact-y[0])/exact;

  // Outputs // 
  cout << "y*(t=2) = "<< y[0] << " y(t=2) = "<< exact <<endl;
  cout << "err= "<< err<< " rel_err= " << rel_err <<endl;
  cout << myintegrator.getName()<<" solve():" << fp_ms.count() << " ms " << endl;
  cout << "t_span= {" << t[0] << "," << t[1] <<"} h= " << h << endl;
  return 0;
}
