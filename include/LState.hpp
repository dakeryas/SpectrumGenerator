#ifndef LORENTZ_STATE
#define LORENTZ_STATE

#include "RState.hpp"

class LState: public RState { //class for representing real states with a width

  double w; //full half width of the state
  
public:
  LState(const FourMomentum& P, int PDG, const std::vector<State*>& children, double w);
  LState(const FourMomentum& P, int PDG, const std::vector<State*>& children, const std::vector<double>& r, double w);
  LState(double m0, int PDG, const std::vector<State*>& children, double w);
  LState(double m0, int PDG, const std::vector<State*>& children, const std::vector<double>& r, double w);
  ~LState(){}
  void AdaptMass(TRandom3&);//sets P to a random rest value according to w /!\ the resulting may be negative ; then acts on the children (VStates), to propagate the random Mass ; 
  
};

#endif