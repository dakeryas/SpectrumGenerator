#include "BState.hpp"

BState::BState(const FourMomentum& P, const std::vector<State*>& children, const Hist2d& kinetic_map):VState(P, children), kinetic_map(kinetic_map){
  
}

BState::BState(const double m0, const std::vector<State*>& children, const Hist2d& kinetic_map):VState(m0, children), kinetic_map(kinetic_map){
  
}

double BState::GetQbeta(){
  
  return GetMass()-GetHeaviestChild()->GetMass();//works for atomic masses
  
}

double BState::GetRndmElectronKinetic(TRandom3& ran){

  return kinetic_map.getRandomFromInterpolatedSliceAtY(ran, GetQbeta());

}

void BState::Process(TRandom3& ran){

  State* electron = GetChildOfPDG(11);
  double E = GetRndmElectronKinetic(ran)+electron->GetMass();
  electron->SetP(FourMomentum(E, sqrt(pow(E,2)-pow(electron->GetMass(),2)), ran)); //pick a random direction on the surface of the sphere of radius p for(px, py, pz)

}