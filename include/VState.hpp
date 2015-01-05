#ifndef VIRTUAL_STATE
#define VIRTUAL_STATE

#include "State.hpp"
#include "InvMassGenerator.hpp"

class VState: public State { //variables used in the decay from a particle/state of mass m1 and width w1 to particles of masses m2 and m3 and corresponding widths

protected:
  void InitialiseStdRatios();//sets r to 1 for unstable states and 0 for the others
  State* GetChildOfPDG(const int pdg);
  State* GetHeaviestChild();
  void TwobodyDecay(TRandom3& ran);
  
public:
  VState(const FourMomentum& P, const vector<State*>& children);
  VState(const FourMomentum& P, const vector<State*>& children, const vector<double>& r);
  VState(const double m0, const vector<State*>& children);
  VState(const double m0, const vector<State*>& children, const vector<double>& r);
  virtual ~VState(){}
  static double TrFct(const double x, const double y, const double z);//triangle function used in particle physics
  unsigned GetNbody() const;
  vector<State*> GetStableChildren();
  bool CanDecay();//checks if the state can kinematically decay
  void Process(TRandom3& ran);
  
};

#endif