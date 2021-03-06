#include "LState.hpp"

LState::LState(const FourMomentum& P, int PDG, const std::vector<State*>& children, double w):RState(P, PDG, children),w(w){
  
  IncorporateRatios(std::vector<double>(children.size(), 1/children.size()));
  
}

LState::LState(const FourMomentum& P, int PDG, const std::vector<State*>& children, const std::vector<double>& r, double w):RState(P, PDG, children, r),w(w){
  
}

LState::LState(double m0, int PDG, const std::vector<State*>& children, double w):LState(FourMomentum(m0), PDG, children, w){
 
  IncorporateRatios(std::vector<double>(children.size(), 1/children.size()));
  
}

LState::LState(double m0, int PDG, const std::vector<State*>& children, const std::vector<double>& r, double w):LState(FourMomentum(m0), PDG, children, r, w){
  
}

void LState::AdaptMass(TRandom3& ran){

  P.SetRest(ran.BreitWigner(m0,w));
  while(GetMass()<m0-10*w || GetMass()>m0+10*w) P.SetRest(ran.BreitWigner(m0,w));//the TH2's are limited to 10w
  SetChildrenP();
  
}