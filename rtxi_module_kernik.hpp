

/*** Header Guard ***/
#ifndef rtxi_module_kernik_H
#define rtxi_module_kernik_H

#include <default_gui_model.h> // Standard RTXI Model
#include "kernik19.hpp" 

/*** rtxi_module_kernik class ***/
class rtxi_module_kernik : public DefaultGUIModel {

 public:
  rtxi_module_kernik(void); // Constructor
  ~rtxi_module_kernik(void); // Destructor
  
  void execute();
  void update(update_flags_t flag); // Non-realtime update function

  
 private:
  kernik model;
  void cell_initialize();
  void calc_I(); // function that calculates the current to inject
  // Dynamic Clamp Experimental Parameters
  int flag;
  double dt,h;
  double sc_k1, sc_to, sc_kr, sc_ks, sc_k1_ishi, sc_CaL;
  double count;
  double bcl;
  double Cm;
  int pace;
  double stimMag;
  double stimLength;
  double bclConverted;
  double stimLengthConverted;
  double out;
  double V,I;
  double ljp;
};

#endif
