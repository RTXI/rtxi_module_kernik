#include <iostream>
#include "kernik19.hpp"
#include <cmath>

#define fastEXP RTmath->fastEXP // shortcut for fastEXP function
#define fastPOW RTmath->fastPOW // shortcut for fastPOW function

using namespace std;

kernik::kernik(){
  initialize();
  RTmath = new RealTimeMath();
}

kernik::~kernik(){
  delete RTmath; 
}

void kernik::initialize(){
  // Integrator Variables 
  integrator_id = 1; // Default Euler
  dt = 0.1; // integration interval 
  h = 0.01; // integrator step
  steps = (int)(ceil(dt/h));
  
  // Physical Parametes
  T = 310.0;   // kelvin (in model_parameters)
  R = 8.314472;   // joule_per_mole_kelvin (in model_parameters)
  F = 96.4853415;   // coulomb_per_mmole (in model_parameters)

  V = -75.596602; // millivolt
  
  Ko = 5.4; // mM (in model_parameters)
  Ki = 130; // mM from Quach 2018
  E_K = R*T/F*log(Ko/Ki); // millivolt

  Nao = 137.0; // mM from Quach 2018
  Nai = 10.0; // mM from Quach 2018

  
  //***** Ca2+ Channel parameters *****//
  
  // Ca2+ Handling System ICs
  // y_Ca[0]=Cai; y_Ca[1]=CaSR; y_Ca[2-4]=d,f,fCa; y_Ca[5-6]=dCaT,fCaT; y_Ca[7-9] = R,O,I
  y_Ca = {3.894379e-04, 3.109153e-01,3.609301e-01, 3.144104e-01,6.186005e-01,
	  3.701136e-01, 4.903681e-01,1.423304e-02, 2.594838e-04,  5.738016e-02};

  // Ca2+ Specific Physical Parametes
  Cao = 2.0; // mM from Quach 2018; Note: Kernik 2019 Cao = 1.8 mM
  E_Ca = 0.5*R*T/F*log(Cao/y_Ca[0]); // millivolt
  Vc = 3960. * ( 16404. / 17498. );
  V_SR = 3960. * ( 1094. / 17498. );

  // L-type Ca2+ parameters
  i_CaL = 0.0;
  tau_fCa = 2.0; // ms
  CaL_params = {3.0802769e-01, 1.2966294e+01, 7.0791460e+00, 4.4909416e-02,-6.9098804e+00,
		5.1258983e-04,-4.9505712e+01, 1.9312112e+03, 5.7300275e+00, 1.6582469e+00,
		1.0046256e+02};
  
  p_CaL_shannonTot = 5.4e-4 + 1.5e-8 + 2.7e-7;
  p_CaL_shannonCap = 5.4e-4 / p_CaL_shannonTot;
  p_CaL_shannonNap = 1.5e-8 / p_CaL_shannonTot;
  p_CaL_shannonKp = 2.7e-7 / p_CaL_shannonTot;
  
  p_CaL_Ca = p_CaL_shannonCap * CaL_params[0];
  p_CaL_Na = p_CaL_shannonNap * CaL_params[0];
  p_CaL_K = p_CaL_shannonKp * CaL_params[0];

  // parameter-dependent values: ICaL
  d3 = CaL_params[3] * CaL_params[1];
  d4 = 1. / ( ( 1. / CaL_params[2] ) + ( 1. / CaL_params[4] ) );
  f3 = CaL_params[7] * CaL_params[5];
  f4 = 1. / ( ( 1. / CaL_params[6] ) + ( 1. / CaL_params[8] ) );

  // T-type Ca2+ parameters
  i_CaT =0.0;
  g_CaT = 0.185; // nS_per_pF (in i_CaT)

  // Na+/Ca2+ Exchanger parameters
  KmCa = 1.38;
  KmNai = 87.5;
  Ksat = 0.1;
  gamma = 0.35*2;
  alpha = 2.5*1.1;
  kNaCa = 1000. * 1.1;
  
  // SR Uptake/SERCA (J_up) parameters
  Kup = 0.00025 * 0.702;  
  VmaxUp = 0.000425 * 0.26;
  
  // SR Leak (J_leak) parameters
  V_leak = 0.00008*0.02;

  // SR Release/RYR (J_rel) parameters
  ks = 12.5;
  koCa = 56320. * 11.43025;
  kiCa = 54. * 0.3425;
  kom = 1.5 * 0.1429;
  kim = 0.001 * 0.5571;
  ec50SR = 0.45;
  MaxSR = 15.;
  MinSR = 1.;

  // Background Calcium (I_bCa) parameters
  g_b_Ca = 0.000592 * 0.62;
  
  // Calcium SL Pump (I_pCa) parameters
  g_PCa = 0.025 * 10.5;
  KPCa = 0.0005;

  // CaSR (millimolar) parameters
  Buf_SR = 10.0 * 1.2;
  Kbuf_SR = 0.3;

  // Cai (millimolar) parameters
  Buf_C = 0.06;
  Kbuf_C = 0.0006;
  
  // K+ Channel parameters
  K1_params = {1.3378578e-01,4.7799497e-01,2.7242756e+01,4.9250233e+00,8.7222376e+00,56.636197};
  Kr_params = {2.1802500e-01,5.7488524e-03,1.3623493e+01,4.7630571e-02,-7.0680874e+00,
	       1.2456641e-02,-2.5994458e+01,3.7342633e+01,2.2091964e+01,5.0000000e+01,
	       0.0000000e+00};
  Ks_params = {0.0077,0.0011655845,66726.839,0.28045891,-18.866972,4.74115e-06};

  to_params = {0.1178333333,0.0553614182,11.684202343,3.9891810804,-11.0471393012,0.0003442309,
	       -17.6344722898,186.7605369097,8.1809338733,0.6967584212,11.2244577239};

  // Ishihara 2009 params
  gK1_ishi = 2.5*pow((Ko/5.4),0.4); // nS/pF
  SPM_in = 0.005; // in mM
  Mg_in = 1.; // in mM
  phi = 0.9; // fraction of high affinity channels
  fo = 0.; // initial condition
  y2=0.; // initial condition
  y1 = {0.5}; // initial condition
  ishi_alpha = 0.;
  ishi_beta = 0.;
  i_K1_ishi = 0.;
    
  // K+ Channel Gates ICs
  Xs = {1.5378828e-01};
  Xrx = {3.0976749e-01,4.5057719e-01}; // Xr1 = Xrx[0]  Xr2 = Xrx[1]
  SR = {0.74680281,0.00026759783};
  
  // parameter-dependent values: IKs
  ks3 = Ks_params[3] * Ks_params[1];
  ks4 = 1. / ( ( 1. / Ks_params[2] ) + ( 1. / Ks_params[4] ) );
  
  // parameter-dependent values: IKr
  Xr1_3 = Kr_params[3] * Kr_params[1];
  Xr2_3 = Kr_params[7] * Kr_params[5];
  Xr1_4 = 1. / ( ( 1. / Kr_params[2] ) + ( 1. / Kr_params[4] ) );
  Xr2_4 = 1. / ( ( 1. / Kr_params[6] ) + ( 1. / Kr_params[8] ) );

  // parameter-dependent values: Ito
  r3 = to_params[3] * to_params[1];
  r4 = 1. / ( ( 1. / to_params[2] ) + ( 1. / to_params[4] ) );
  s3 = to_params[7] * to_params[5];
  s4 = 1. / ( ( 1. / to_params[6] ) + ( 1. / to_params[8] ) );


}	

void kernik::dXrx_func(array<double,2>* y, array<double,2>* dy){
  // Rapid Delayed Rectifier Current (Ikr):
  // Xr1 (dimensionless) (activation in i_Kr_Xr1)
  alpha_Xr1 = Kr_params[1] * fastEXP( ( V ) / Kr_params[2] );
  beta_Xr1  = Xr1_3 * fastEXP( ( V ) / Xr1_4 );
  Xr1_inf   = alpha_Xr1 / ( alpha_Xr1 + beta_Xr1 );
  tau_Xr1   = ( ( 1. / ( alpha_Xr1 + beta_Xr1 ) ) + Kr_params[9] );
  dy->data()[0] = ( Xr1_inf - y->data()[0] ) / tau_Xr1;
  
  // Xr2 (dimensionless) (inactivation in i_Kr_Xr2)
  alpha_Xr2 = Kr_params[5] * fastEXP( ( V ) / Kr_params[6] );
  beta_Xr2  = Xr2_3 * fastEXP( ( V ) / Xr2_4 );
  Xr2_inf   = alpha_Xr2 / ( alpha_Xr2 + beta_Xr2 );
  tau_Xr2   = ( ( 1. / ( alpha_Xr2 + beta_Xr2 ) ) + Kr_params[10] ) ;
  dy->data()[1] = ( Xr2_inf - y->data()[1] ) / tau_Xr2;
  
}


void kernik::dXs_func(array<double,1>* y, array<double,1>* dy) {
  // Xs (dimensionless) (activation in i_Ks)
  alpha_Xs = Ks_params[1] * fastEXP( ( V ) / Ks_params[2] );
  beta_Xs  = ks3 * fastEXP( ( V ) / ks4 );
  Xs_inf   = alpha_Xs / ( alpha_Xs + beta_Xs );
  tau_Xs   = ( 1. / ( alpha_Xs + beta_Xs ) ) + Ks_params[5];
  dy->data()[0] = ( Xs_inf - y->data()[0] ) / tau_Xs;

}

void kernik::dSR_func(array<double,2>* y, array<double,2>* dy){
  // 17: s (dimensionless) (inactivation in i_to)
  alpha_s = to_params[5] * fastEXP( ( V ) / to_params[6] );
  beta_s = s3 * fastEXP( ( V ) / s4);
  s_inf = alpha_s / ( alpha_s + beta_s );
  tau_s = ( ( 1. / ( alpha_s + beta_s ) ) + to_params[10] );
  dy->data()[0]  = ( s_inf - y->data()[0] ) / tau_s;

  // 18: r (dimensionless) (activation in i_to)
  alpha_r = to_params[1] * fastEXP( ( V ) / to_params[2] );
  beta_r = r3 * fastEXP( ( V ) / r4 );
  r_inf = alpha_r / ( alpha_r + beta_r );
  tau_r = ( 1. / ( alpha_r + beta_r ) ) + to_params[9];
  dy->data()[1] = ( r_inf - y->data()[1] ) / tau_r;
}

void kernik::dy1_func(array<double,1>* y, array<double,1>* dy){
  // y1 from Ishihara 2009 IK1 model
  ishi_alpha = (0.17*fastEXP(-0.07*((V-E_K) + 8*Mg_in)))/(1+0.01*fastEXP(0.12*(V-E_K)+8*Mg_in));
  ishi_beta = (SPM_in*280*fastEXP(0.15*(V-E_K)+8*Mg_in))/(1+0.01*fastEXP(0.13*(V-E_K)+8*Mg_in));
  Kd_spm_l = 0.04*fastEXP(-(V-E_K)/9.1);
  Kd_mg = 0.45*fastEXP(-(V-E_K)/20);
  fo = 1/(1 + (Mg_in/Kd_mg));
  y2 = 1/(1 + SPM_in/Kd_spm_l);

  dy->data()[0]=(ishi_alpha*(1-y->data()[0]) - ishi_beta*fastPOW(fo,3)*y->data()[0]);
  
}

// Ca2+ Handling System

void kernik::dy_Ca_func(array<double,10>* y, array<double,10>* dy) {
  // CaL gates
  // y_Ca[2]: d (dimensionless) (activation in i_CaL)
  alpha_d = CaL_params[1] * fastEXP(V / CaL_params[2]);
  beta_d = d3 * fastEXP( V / d4 );
  d_inf = alpha_d / ( alpha_d + beta_d );
  tau_d = ( ( 1. / ( alpha_d + beta_d ) ) + CaL_params[9] );
  dy->data()[2] = ( d_inf - y->data()[2] ) / tau_d;

  // y_Ca[3]: f (dimensionless) (activation in i_CaL)
  alpha_f = CaL_params[5] * fastEXP( ( ( V ) ) / CaL_params[6] );
  beta_f = f3 * fastEXP( ( ( V ) ) / f4 );
  f_inf = alpha_f / ( alpha_f + beta_f );
  tau_f = ( ( 1. / ( alpha_f + beta_f ) ) + CaL_params[10] );
  dy->data()[3] = ( f_inf - y->data()[3] ) / tau_f;

  // y_Ca[4]:  (dimensionless) (activation in i_CaL)
  alpha_fCa = 1.0 / ( 1.0 + fastPOW( ( ( 1.2 * y->data()[0] ) / 0.000325 ), 8. ) ) ; //^8.0);
  beta_fCa = 0.1 / ( 1.0 + fastEXP( ( 1.2 * y->data()[0] - 0.0005 ) / 0.0001 ) );
  gamma_fCa = 0.2 / ( 1.0 + fastEXP( ( 1.2 * y->data()[0] - 0.00075 ) / 0.0008 ) );
  fCa_inf = ( ( alpha_fCa + beta_fCa + gamma_fCa + 0.23 ) / ( 1.46 ) );
  if( fCa_inf > y->data()[4] && V > -60. ) {
    k_fca = 0;
  } else {
    k_fca = 1;
  } //end
  dy->data()[4] = k_fca * ( fCa_inf - y->data()[4] ) / tau_fCa;

  // Calculate i_CaL
  ibarca= (p_CaL_Ca*4.0*V*F*F/(R*T) * (0.341*y->data()[0]*fastEXP(2.0*V*F/(R*T))-0.341*Cao)
	   / (fastEXP(2.0*V*F/(R*T))-1.0));
  i_CaL_Ca =  ibarca*y->data()[2]*y->data()[3]*y->data()[4];
    
  ibarna = (p_CaL_Na*V*F*F/(R*T) * (0.75*Nai*fastEXP(V*F/(R*T))-0.75*Nao)
	   / (fastEXP(V*F/(R*T))-1.0));
  i_CaL_Na = ibarna*y->data()[2]*y->data()[3]*y->data()[4];
    
  ibark = (p_CaL_K*V*F*F/(R*T) * (0.75*Ki*fastEXP(V*F/(R*T)) -0.75*Ko)
	   / ( fastEXP( V * F / ( R * T ) )-1.0));
  i_CaL_K = ibark*y->data()[2]*y->data()[3]*y->data()[4] ;

  i_CaL = i_CaL_Ca + i_CaL_Na + i_CaL_K;

  // CaT gates
  // 5: dCaT (activation in i_CaT)
  dcat_inf = 1. / ( 1. + fastEXP( -( ( V ) + 26.3 ) / 6. ) );
  tau_dcat = 1. / ( 1.068 * fastEXP( ( ( V ) + 26.3 ) / 30. )
			   + 1.068 * fastEXP( -( ( V ) + 26.3 ) / 30. ) );
  dy->data()[5] = ( dcat_inf - y->data()[5] ) / tau_dcat;
  
  // 6: fCaT (inactivation in i_CaT)
  fcat_inf = 1. / ( 1. + fastEXP( ( ( V ) + 61.7 ) / 5.6 ) ) ;
  tau_fcat = 1. / ( 0.0153 * fastEXP( -( ( V ) + 61.7 ) / 83.3 )
			   + 0.015 * fastEXP( ( ( V ) + 61.7 ) / 15.38 ) );
  dy->data()[6] = ( fcat_inf - y->data()[6] ) / tau_fcat;
  
  // Calculate i_CaT
  i_CaT = g_CaT * ( V - E_Ca ) * y->data()[5] * y->data()[6]; // nS_per_pF (in i_CaT)

  // Na+/Ca2+ Exchanger current (INaCa) pA_per_pF (in i_NaCa)
  i_NaCa=kNaCa*((fastEXP(gamma*V*F/(R*T))*(Nai*Nai*Nai)*Cao)
		 -(fastEXP((gamma-1.0)*V*F/(R*T))*(Nao*Nao*Nao) * y->data()[0] * alpha))
	  /(((KmNai*KmNai*KmNai)+(Nao*Nao*Nao))*(KmCa+Cao)*(1.0+Ksat*fastEXP((gamma-1.0)*V*F/(R*T))));

  // SR Uptake/SERCA (J_up)
  j_up = VmaxUp / (1.0 + Kup*Kup / (y->data()[0]*y->data()[0])); // millimolar_per_milisecond

  // SR Leak (J_leak)
  j_leak = ( y->data()[1] - y->data()[0] ) * V_leak; // millimolar_per_milisecond

  // SR Release/RYR (J_rel)
  kCaSR = MaxSR-(MaxSR-MinSR)/(1.+fastPOW((ec50SR/y->data()[1]),2.5));
  koSRCa = koCa / kCaSR;
  kiSRCa = kiCa * kCaSR;
  RI = 1. - y->data()[7] - y->data()[8] - y->data()[9] ;

  dy->data()[7] = ( ( kim * RI - kiSRCa * y->data()[0] * y->data()[7] )
	     - ( koSRCa * y->data()[0] * y->data()[0] * y->data()[7]  - kom * y->data()[8] ) );   // R
  dy->data()[8] = ( ( koSRCa * y->data()[0] * y->data()[0] * y->data()[7] - kom * y->data()[8] )
	     - ( kiSRCa * y->data()[0] * y->data()[8] - kim * y->data()[9] ) ); // O
  dy->data()[9] = ( ( kiSRCa * y->data()[0] * y->data()[8] - kim * y->data()[9] )
	     - ( kom * y->data()[9] - koSRCa * y->data()[0] * y->data()[0] * RI ) );   // I

  j_rel = ks*y->data()[8]*( y->data()[1]-y->data()[0])*(V_SR/Vc); // millimolar_per_milisecond

  // Background Calcium (I_bCa)
  i_b_Ca = g_b_Ca * ( V - E_Ca ); // nS_per_pF (in i_b_Ca)

  // Calcium SL Pump (I_pCa)
  i_PCa = g_PCa * y->data()[0] / ( y->data()[0] + KPCa ); // nS_per_pF (in i_b_Ca)
  
  // 1: CaSR (millimolar)
  Ca_SR_bufSR = (1./(1.0+Buf_SR * Kbuf_SR/((y->data()[1] + Kbuf_SR)*(y->data()[1]+Kbuf_SR))));

  dy->data()[1] = Ca_SR_bufSR * Vc/V_SR * (j_up-(j_rel+j_leak));

  // 0: Cai (millimolar)
  Cai_bufc = 1. / (1.0+Buf_C*Kbuf_C / ((y->data()[0]+Kbuf_C)*(y->data()[0]+Kbuf_C)));
  
  dy->data()[0] = Cai_bufc*(j_leak-j_up+j_rel-(i_CaL_Ca+i_CaT+i_b_Ca+i_PCa-2.*i_NaCa)*60.0/(2.0*Vc*F));
}


void kernik::scale_IK1(double sc){
  // No solve on IK1
  alpha_xK1 = K1_params[1] * fastEXP( ( V + K1_params[3] ) / K1_params[2] );
  beta_xK1  = fastEXP( ( V + K1_params[5] ) / K1_params[4] );
  XK1_inf   = alpha_xK1 / ( alpha_xK1 + beta_xK1 );
  i_K1 = K1_params[0] * XK1_inf * (V-E_K) * sqrt(Ko/5.4); // pA_per_pF
  Iout += sc * i_K1;
}

void kernik::scale_IKs(double sc){
  i_Ks = Ks_params[0] * (V-E_K) * Xs[0] * Xs[0]; // pA_per_pF
  Iout += sc * i_Ks;
  Xs_solve(); // Update Xs for next step
}

void kernik::scale_IKr(double sc){
  i_Kr = Kr_params[0] * (V-E_K) * Xrx[0]* Xrx[1] * sqrt(Ko/5.4); // pA_per_pF
  Iout += sc * i_Kr;
  Xrx_solve(); // Update Xrx for next step
}

void kernik::scale_Ito(double sc){
  i_to = to_params[0] * (V-E_K) * SR[0] * SR[1]; // pA_per_pF
  Iout += sc * i_to;
  SR_solve();
}

void kernik::scale_IK1_ishi(double sc){
  i_K1_ishi = gK1_ishi*(V-E_K)*(phi*fo*y1[0] + (1-phi)*y2); // in pA_per_pF
  Iout += sc * i_K1_ishi;
  y1_solve();
}

void kernik::scale_ICaL(double sc){
  Iout += sc * i_CaL;
  Ca_solve();
}

void kernik::Xrx_solve(){
  array<double,2> k1,k2,y_half_h;
  if (integrator_id == 1) { // Euler Method
    for (int i=0;i<steps;++i){
      dXrx_func(&Xrx,&k1); // Calculate K1 from Y0
      for (int j=0;j<(int)Xrx.size();++j){
	Xrx[j] += h* (k1[j]); // Update Y
      }
    }
  } else if (integrator_id == 2) { // RK2
    for (int i=0;i<steps;++i){
      dXrx_func(&Xrx,&k1); // Calculates k1
      for (int j=0;j<(int)Xrx.size();++j){
	y_half_h[j] = Xrx[j] + (h*0.5)* (k1[j]);
      }
      dXrx_func(&y_half_h,&k2); // Calculates k2 from half step
      for (int j=0;j<(int)Xrx.size();++j){
	Xrx[j] += h* (k2[j]);
      }    
    }
  } 
}


void kernik::Xs_solve(){
  array<double,1> k1,k2,y_half_h;
  if (integrator_id == 1) { // Euler Method
    for (int i=0;i<steps;++i){
      dXs_func(&Xs,&k1); // Calculate K1 from Y0
      for (int j=0;j<(int)Xs.size();++j){
	Xs[j] += h* (k1[j]); // Update Y
      }
    }
  } else if (integrator_id == 2) { // RK2
    for (int i=0;i<steps;++i){
      dXs_func(&Xs,&k1); // Calculates k1
      for (int j=0;j<(int)Xs.size();++j){
	y_half_h[j] = Xs[j] + (h*0.5)* (k1[j]);
      }
      dXs_func(&y_half_h,&k2); // Calculates k2 from half step
      for (int j=0;j<(int)Xs.size();++j){
	Xs[j] += h* (k2[j]);
      }    
    }
  } 
}

void kernik::SR_solve(){
  array<double,2> k1, k2, y_half_h;
  if (integrator_id == 1) { // Euler Method
    for (int i=0;i<steps;++i){
      dSR_func(&SR,&k1); // Calculate K1 from Y0
      for (int j=0;j<(int)SR.size();++j){
	SR[j] += h* (k1[j]); // Update Y
      }
    }
  } else if (integrator_id == 2) { // RK2
    for (int i=0;i<steps;++i){
      dSR_func(&SR,&k1); // Calculates k1
      for (int j=0;j<(int)SR.size();++j){
	y_half_h[j] = SR[j] + (h*0.5)* (k1[j]);
      }
      dSR_func(&y_half_h,&k2); // Calculates k2 from half step
      for (int j=0;j<(int)SR.size();++j){
	SR[j] += h* (k2[j]);
      }    
    }
  } 
}

void kernik::y1_solve(){
  array<double,1> k1, k2, y_half_h;
  if (integrator_id == 1) { // Euler Method
    for (int i=0;i<steps;++i){
      dy1_func(&y1,&k1); // Calculate K1 from Y0
      for (int j=0;j<(int)y1.size();++j){
	y1[j] += h* (k1[j]); // Update Y
      }
    }
  } else if (integrator_id == 2) { // RK2
    for (int i=0;i<steps;++i){
      dy1_func(&y1,&k1); // Calculates k1
      for (int j=0;j<(int)y1.size();++j){
	y_half_h[j] = y1[j] + (h*0.5)* (k1[j]);
      }
      dy1_func(&y_half_h,&k2); // Calculates k2 from half step
      for (int j=0;j<(int)y1.size();++j){
	y1[j] += h* (k2[j]);
      }    
    }
  } 
}

void kernik::Ca_solve(){
  array<double,10> k1, k2, y_half_h;
  if (integrator_id == 1) { // Euler Method
    for (int i=0;i<steps;++i){
      dy_Ca_func(&y_Ca,&k1); // Calculate K1 from Y0
      for (int j=0;j<(int)y_Ca.size();++j){
	y_Ca[j] += h* (k1[j]); // Update Y
      }
    }
  } else if (integrator_id == 2) { // RK2
    for (int i=0;i<steps;++i){
      dy_Ca_func(&y_Ca,&k1); // Calculates k1
      for (int j=0;j<(int)y_Ca.size();++j){
	y_half_h[j] = y_Ca[j] + (h*0.5)* (k1[j]);
      }
      dy_Ca_func(&y_half_h,&k2); // Calculates k2 from half step
      for (int j=0;j<(int)y_Ca.size();++j){
	y_Ca[j] += h* (k2[j]);
      }    
    }
  }
  E_Ca = 0.5*R*T/F*log(Cao/y_Ca[0]); // millivolt
}
