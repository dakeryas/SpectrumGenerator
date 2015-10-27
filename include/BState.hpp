#ifndef BETA_STATE
#define BETA_STATE

#include <Hist2d.hpp>
#include "VState.hpp"

class BState: public VState { //variables used in the decay from a particle/state of mass m1 and width w1 to particles of masses m2 and m3 and corresponding widths
  
  Hist2d kinetic_map;//kinetic energy distribution of the electron as a function of the effective qBeta (which changes for resonnant states)
  double GetQbeta();
  double GetRndmElectronKinetic(TRandom3& ran);//get a random kinetic energy for the electron; slice a TH2 at Yvalue and returns a random number distributed according to the resulting slice, ie the TH1 named "chosent dist"
  
public:
  BState(const FourMomentum& P, const std::vector<State*>& children, const Hist2d& kinetic_map);
  BState(double m0, const std::vector<State*>& children, const Hist2d& kinetic_map);
  ~BState(){}
  void Process(TRandom3& ran);
 
};

#endif