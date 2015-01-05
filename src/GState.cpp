#include "GState.hpp"

GState::GState(const FourMomentum& P, const int PDG, const vector<State*>& children, const double w):RState(P, PDG, children),sigma(w/2/sqrt(2*log(2))){
  
  IncorporateRatios(vector<double>(children.size(), 1/children.size()));
  
}

GState::GState(const FourMomentum& P, const int PDG, const vector<State*>& children, const vector<double>& r, const double w):RState(P, PDG, children, r),sigma(w/2/sqrt(2*log(2))){
  
}

GState::GState(const double m0, const int PDG, const vector<State*>& children, const double w):GState(FourMomentum(m0), PDG, children, w){
 
  IncorporateRatios(vector<double>(children.size(), 1/children.size()));
  
}

GState::GState(const double m0, const int PDG, const vector<State*>& children, const vector<double>& r, const double w):GState(FourMomentum(m0), PDG, children, r, w){
  
}

void GState::AdaptMass(TRandom3& ran){

  P.SetRest(ran.Gaus(m0,sigma));
  while(GetMass()<m0-10*sigma || GetMass()>m0+10*sigma) P.SetRest(ran.Gaus(m0,sigma));//the TH2's are limited to 10w
  SetChildrenP();
  
}