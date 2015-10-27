#include "RState.hpp"

RState::RState(const FourMomentum& P, int PDG):State(P),PDG(PDG){
  
}

RState::RState(const FourMomentum& P, int PDG, const std::vector<State*>& children):State(P, children), PDG(PDG){
  
  IncorporateRatios(std::vector<double>(children.size(), 1/children.size()));
  
}

RState::RState(const FourMomentum& P, int PDG, const std::vector<State*>& children, const std::vector<double>& r):State(P, children, r),PDG(PDG){
  
}

RState::RState(double m0, int PDG):RState(FourMomentum(m0), PDG){

}

RState::RState(double m0, int PDG, const std::vector<State*>& children):RState(FourMomentum(m0), PDG, children){
 
  IncorporateRatios(std::vector<double>(children.size(), 1/children.size()));
  
}

RState::RState(double m0, int PDG, const std::vector<State*>& children, const std::vector<double>& r):RState(FourMomentum(m0), PDG, children, r){
  
}

void RState::SetChildrenP(){
  
  for(State* child : children) child->SetP(P);
  
}

bool RState::IsStable(){
  
  return rsum.empty();
  
}

int RState::GetPDG() const {
  
  return PDG;
  
}

void RState::Process(TRandom3&){
  
  SetChildrenP();
  
}