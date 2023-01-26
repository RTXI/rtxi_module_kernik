// header guard // 
#ifndef KERNIK_H
#define KERNIK_H

#include <vector>
#include <array>
#include <cmath>
#include "./include/RealTimeMath.h" // RealTimeMath library 

using namespace std;


class kernik
{
public:
  kernik();
  ~kernik();

  void initialize();

  void scale_IK1(double);
  void scale_Ito(double);
  void scale_IKr(double);
  void scale_IKs(double);
  void scale_IK1_ishi(double);
  void scale_ICaL(double);

  double get_Iout(){return Iout;}
  void reset_Iout(){this->Iout=0.0;}

  double getV(){return V;}
  void setV(double V){this->V=V;}
  
  void set_h(double h){
    this->h=h;
    steps = (int)(ceil(dt/h));
  }
  void set_dt(double dt){
    this->dt=dt;
    steps = (int)(ceil(dt/h));
  }


  void set_integrator(int id){integrator_id=id;}

private:
  RealTimeMath *RTmath;

  // Reversal Potentials
  double V; // model units in mV 
  double T,R,F;
  double Ki,Ko;
  double Nai,Nao;
  double Cao;
  double E_K,E_Ca,E_Na;

  // K+ channels 
  
  array<double,6> K1_params;
  array<double,11> Kr_params;
  array<double,6> Ks_params;
  array<double,11> to_params;

  // IK1
  double alpha_xK1,beta_xK1,XK1_inf,i_K1;

  // IKr
  array<double,2> Xrx;
  void dXrx_func(array<double,2>* y, array<double,2>* dy);
  void Xrx_solve();
  double i_Kr;
  double Xr1_3,Xr2_3,Xr1_4,Xr2_4;
  double alpha_Xr1,beta_Xr1,Xr1_inf,tau_Xr1;
  double alpha_Xr2,beta_Xr2,Xr2_inf,tau_Xr2;

  // IKs
  array<double,1> Xs;
  void dXs_func(array<double,1>* y, array<double,1>* dy);
  void Xs_solve();
  double i_Ks;
  double ks3,ks4;
  double alpha_Xs,beta_Xs,Xs_inf,tau_Xs;
  
  // Ito
  array<double,2> SR;
  double i_to;
  double r3,r4,s3,s4;
  double alpha_s, beta_s, s_inf, tau_s;
  double alpha_r, beta_r, r_inf, tau_r;
  void dSR_func(array<double,2>* y, array<double,2>* dy);
  void SR_solve();

  // IK1 Ishihara
  array<double,1> y1;
  double i_K1_ishi;
  double ishi_alpha,ishi_beta,Kd_spm_l,Kd_mg,fo,y2;
  double gK1_ishi, SPM_in, Mg_in, phi;
  void dy1_func(array<double,1>* y, array<double,1>* dy);
  void y1_solve();

  // Ca2+ System
  array<double,10> y_Ca;
  void dy_Ca_func(array<double,10>* y, array<double,10>* dy);
  void Ca_solve();
  double Vc, V_SR;
  double Buf_SR, Kbuf_SR, Ca_SR_bufSR;  // CaSR (millimolar) parameters
  double Buf_C, Kbuf_C, Cai_bufc;   // Cai (millimolar) parameters

  // L-type Ca2+
  array<double,11> CaL_params;
  double i_CaL;
  double d3, d4, f3, f4;
  double alpha_d, beta_d, d_inf, tau_d;
  double alpha_f, beta_f, f_inf, tau_f;
  double alpha_fCa, beta_fCa, gamma_fCa, fCa_inf, tau_fCa, k_fca;
  double p_CaL_shannonTot, p_CaL_shannonCap, p_CaL_shannonNap, p_CaL_shannonKp;
  double p_CaL_Ca, p_CaL_Na, p_CaL_K, ibarca, i_CaL_Ca, ibarna, i_CaL_Na, ibark, i_CaL_K;

  // T-type Ca2+
  double dcat_inf,tau_dcat,fcat_inf,tau_fcat;
  double g_CaT, i_CaT;

  // Na+/Ca2+ Exchanger
  double KmCa, KmNai, Ksat, gamma, alpha, kNaCa;
  double i_NaCa;
  
  // SR Uptake/SERCA (J_up)
  double Kup, VmaxUp;
  double j_up;

  // SR Leak (J_leak)
  double V_leak;
  double j_leak;

  // SR Release/RYR (J_rel)
  double ks, koCa, kiCa, kom, kim, ec50SR, MaxSR, MinSR;
  double kCaSR, koSRCa, kiSRCa, RI;
  double j_rel;

  // Background Calcium (I_bCa)
  double g_b_Ca;
  double i_b_Ca;

  // Calcium SL Pump (I_pCa)
  double g_PCa, KPCa, i_PCa;
  
  // Integrator variables 
  int integrator_id;
  double h;
  double dt;
  int steps;
  double Iout; // Variable to store total Current to inject
  
};


#endif
