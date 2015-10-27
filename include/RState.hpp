#ifndef REAL_STATE
#define REAL_STATE

#include "State.hpp"

class RState: public State { //variables used in the decay from a particle/state of mass m1 and width w1 to particles of masses m2 and m3 and corresponding widths

protected:
  int PDG;
  void SetChildrenP();//propagate the kinetic energy (through P) to boost indeed the virtual children that will effectively decay;propagate the random mass
  
public:
  RState(const FourMomentum& P, int PDG);
  RState(const FourMomentum& P, int PDG, const std::vector<State*>& children);
  RState(const FourMomentum& P, int PDG, const std::vector<State*>& children, const std::vector<double>& r);
  RState(double m0, int PDG);
  RState(double m0, int PDG, const std::vector<State*>& children);
  RState(double m0, int PDG, const std::vector<State*>& children, const std::vector<double>& r);
  ~RState(){}
  bool IsStable();
  int GetPDG() const;
  void Process(TRandom3&);
  
};

#endif