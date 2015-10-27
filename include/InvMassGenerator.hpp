#ifndef INV
#define INV

#include "VState.hpp"

class VState;//declaration because of the mutual inclusion

class InvMassGenerator {
  
  unsigned n;//number of bodies involved in the decay
  std::vector<double> m;//to store the masses of the n particles to decay plus the mass of the State s to decay
  TRandom3 ran;
  std::vector<double> Mmin;//lower bound for the invariant mass set
  std::vector<double> M;//invariant mass set to be returned
  double y;//the current y to be compared to GetDistribution()
  double dmax;//value set by ComputeMaxDistribution() when the phase space is largest, used to renormalize GetDistribution()
  void InitialiseMasses(const VState* s);
  double GetDistribution();//returns the "probability" of the current mass set M
  void ComputeMaxDistribution();//get the distribtion of the mass set for which the phase space is largest
  void Generate();

public:
  InvMassGenerator();
  InvMassGenerator(const VState* s, TRandom3& ran);
  std::vector<double> GetMassSet(); //returns a copy of M without the first (mother) and the last elements (last child)
 
};

#endif