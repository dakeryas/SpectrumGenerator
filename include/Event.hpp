#ifndef EVENT_H
#define EVENT_H

#include <fstream>
#include "Point.hpp"
#include "LState.hpp"
#include "GState.hpp"
#include "BState.hpp"

class Event {

  State* start;
  TRandom3 ran;
  Point position; //to store the whole vertex position(only one per event at the moment)
  std::vector<State*> road;
  std::vector<State*> stable;
  void ChooseDecayRoad();//fills 'road' with the states on the decay path of 'start'
  void AdaptRoad();//set a random mass for the real states and check that the configuraiton is kinematically possible
  void DecayRoadStates();//decay all the virtual states and fills the stable states on the way
  
public:
  Event();
  Event(State* start, const Point& position);
  Event(State* start, double R, double h);
  Event(State* start, double R1, double R2, double h1, double h2);
  void SetPosition(const Point& position);
  void SetInCylinder(double R, double h);
  void SetInHollowCylinder(double R1, double R2, double h1, double h2);
  void Process();
  void Fill(TH1F& h);
  void Fill(TH1F& h, int PDG);//fill only particles with pdg = PDG
  const Point& GetPosition() const;
  const std::vector<State*>& GetStableStates() const;
  
  
};

std::ostream& operator<<(std::ostream& output, const Event& Ev);


#endif