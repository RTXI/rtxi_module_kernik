#include <iostream>
#include "rtxi_module_kernik.hpp"

using namespace std;

// Required by RTXI
extern "C" Plugin::Object *createRTXIPlugin(void) {
    return new rtxi_module_kernik();
}

// These variables will automatically populate the Default GUI
// { "GUI LABEL", "GUI TOOLTIP", "VARIABLE TYPE",},
static DefaultGUIModel::variable_t vars[] = {    
    { "Voltage Input (V)", "Voltage Input (V)", DefaultGUIModel::INPUT, }, // input(0)
    { "Current Output (A)", "Current output (A)", DefaultGUIModel::OUTPUT, }, // output(0)
    { "Vm", "Membrane Potential (mV)", DefaultGUIModel::STATE, },
    { "I", "Injected Current (pA)", DefaultGUIModel::STATE, },
    { "BCL", "Basic cycle length for pacing (ms)", DefaultGUIModel::PARAMETER | DefaultGUIModel::INTEGER, },
    { "Stim Mag", "Stimulus magnitude (nA)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Stim Length", "Stimulus length (ms)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Pace", "(0: Pace off) (1: Pace on)", DefaultGUIModel::PARAMETER | DefaultGUIModel::INTEGER, },
    { "Cm", "Membrane capacitance (pF)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "IK1", "scaling coefficient (dimensionless)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Ito", "scaling coefficient (dimensionless)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "IKr", "scaling coefficient (dimensionless)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "IKs", "scaling coefficient (dimensionless)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "IK1_ishi", "scaling coefficient (dimensionless)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "ICaL", "scaling coefficient (dimensionless)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "ljp", "Liquid Junction Potential (mV)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
};

// Necessary RTXI variable *DO NOT EDIT*
static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);


/*** Model Constructor ***/
rtxi_module_kernik::rtxi_module_kernik(void) :
  DefaultGUIModel("rtxi_module_kernik", ::vars, ::num_vars) {
    
  createGUI(vars, num_vars); // Function call to create graphical interface *can be overloaded for custom GUI*
  show(); // Show GUI
  model.initialize(); // Calls the model initialize
  model.set_integrator(1);
  cell_initialize(); // Initializes Parameters for dynamic clamp
  refresh(); // Refresh GUI
}
// End constructor

/*** Model Destructor ***/
rtxi_module_kernik::~rtxi_module_kernik(void) {
} // End destructor

/*** Non-Realtime Update Function ***/
void rtxi_module_kernik::update(DefaultGUIModel::update_flags_t flag) {
  switch (flag) {

  case UNPAUSE: // Called when pause button is untoggled
    dt = RT::System::getInstance()->getPeriod()*1e-6; // RTXI thread period; converts to ms (from ns)
    h = dt/100.0;
    model.set_dt(dt);
    model.set_h(h);
    
    //scale = getParameter("Scale Factor").toDouble();
    break;
    
  case PAUSE: // Called when pause button is toggled
    
    out = 0;
    
    break;
    
  case MODIFY: // Called when modify button is hit                
    // Update Paramters
    bcl = getParameter("BCL").toInt();
    stimMag = getParameter("Stim Mag").toDouble();
    stimLength = getParameter("Stim Length").toDouble();
    pace = getParameter ("Pace").toInt();
    Cm = getParameter("Cm").toDouble();
    ljp = getParameter("ljp").toDouble();
    sc_k1 = getParameter("IK1").toDouble();
    sc_to = getParameter("Ito").toDouble();
    sc_kr = getParameter("IKr").toDouble();
    sc_ks = getParameter("IKs").toDouble();
    sc_k1_ishi = getParameter("IK1_ishi").toDouble();
    sc_CaL = getParameter("ICaL").toDouble();

    dt = RT::System::getInstance()->getPeriod()*1e-6; // RTXI thread period (ms; from ns)
    h = dt/100.0;
    model.set_dt(dt);
    model.set_h(h);
    bclConverted = bcl/dt; 
    stimLengthConverted = stimLength / dt;
    
    break;
  default :
    abort();
  }

}


/*** Parameter Initialization ***/
void rtxi_module_kernik::cell_initialize() {

  // default initialization
  flag = 0;
  ljp = -2.8; // in mV based off of Quach 2018 solutions
  pace = 0; // 0 = no pacing 
  bcl = 1000; // in ms
  stimMag = 0.5; // in nA
  stimLength = 2; // in ms
  Cm = 60;
  V = model.getV();
  I =0.0;
  sc_k1 = 0;
  sc_to = 0;
  sc_kr = 0;
  sc_ks = 0;
  sc_k1_ishi = 0;
  sc_CaL = 0;

  // Set States and Parameters, connects GUI with actual variables
  setState("Vm", V);
  setState("I", I);
  setParameter("ljp", ljp);
  setParameter("BCL", bcl);
  setParameter("Stim Mag", stimMag);
  setParameter("Stim Length", stimLength);
  setParameter("Pace", pace);
  setParameter("Cm", Cm);
  setParameter("IK1", sc_k1);
  setParameter("Ito", sc_to);
  setParameter("IKr", sc_kr);
  setParameter("IKs", sc_ks);
  setParameter("IK1_ishi", sc_k1_ishi);
  setParameter("ICaL", sc_CaL);

  update(MODIFY); // Update user input parameters and rate dependent variables
}

void rtxi_module_kernik::calc_I(){ // function that calculates the current to inject

  if ( pace == 1 ) { // Pacing: directly add current to out
    if ( count < stimLengthConverted ) {
      out += (stimMag * 1.0e-9); // Convert units (nA -> A)
    }
    if ( count >= bclConverted ) {// Reset count at end of beat
      count = 0;
    } else {
      count++;
    }
  }

  // Current Scaling
  if (sc_k1!=0.0) {
    model.scale_IK1(sc_k1);
  }
  if (sc_to!=0.0) {
    model.scale_Ito(sc_to);
  }
  if (sc_kr!=0.0) {
    model.scale_IKr(sc_kr);
  }
  if (sc_ks!=0.0) {
    model.scale_IKs(sc_ks);
  }
  if (sc_k1_ishi!=0.0) {
    model.scale_IK1_ishi(sc_k1_ishi);
  }
  if (sc_CaL!=0.0) {
    model.scale_ICaL(sc_CaL);
  }
  
  I=model.get_Iout(); // model units in pA_per_pF
  model.reset_Iout();
  
}

/*** Realtime Execute Function ***/
void rtxi_module_kernik::execute() {
  out = 0; 
  V = input(0) * 1000 + ljp;
  model.setV(V);
  calc_I();
  out += (-1 * Cm * I * 1e-12); // -1 for sign change for amplifier polarity; Convert to A
  output(0) = out; 

}
