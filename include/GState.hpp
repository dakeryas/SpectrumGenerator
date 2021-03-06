#ifndef GAUSS_STATE
#define GAUSS_STATE

#include "RState.hpp"

class GState: public RState { //class for representing real states with a width

  double sigma; //sigma of the gaussian used to pick the mass = 1/(2 ln2 * fwhm)
  
public:
  GState(const FourMomentum& P, int PDG, const std::vector<State*>& children, double w);
  GState(const FourMomentum& P, int PDG, const std::vector<State*>& children, const std::vector<double>& r, double w);
  GState(double m0, int PDG, const std::vector<State*>& children, double w);
  GState(double m0, int PDG, const std::vector<State*>& children, const std::vector<double>& r, double w);
  ~GState(){}
  void AdaptMass(TRandom3&);//sets P to a random rest value according to w /!\ the resulting may be negative ; then acts on the children (VStates), to propagate the random Mass ; 
  
};

#endif