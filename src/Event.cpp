#include "Event.hpp"

ostream& operator<<(ostream& output, const Event& Ev){//for input masses in MeV sets the file into GeV

  output<<Ev.GetStableStates().size()<<"\n";
  for(State* state : Ev.GetStableStates())
    output<<"1"<<" "
    <<state->GetPDG()<<" 0"<<" 0"<<" "<<state->GetP().Getpx()*1.e-3<<" "<<state->GetP().Getpy()*1.e-3<<" "<<state->GetP().Getpz()*1.e-3<<" "<<state->GetMass()*1.e-3<<" "<<"0"<<" "
    <<Ev.GetPosition().x<<" "<<Ev.GetPosition().y<<" "<<Ev.GetPosition().z<<"\n";
  
  return output;
  
}

Event::Event(){
    
}

Event::Event(State* start, const Point& position):start(start), ran(TRandom3 (0)), position(position){
 
}

Event::Event(State* start, const double R, const double h):start(start),ran(TRandom3 (0)), position(Point(ran, R, h)){
  
}

Event::Event(State* start, double R1, double R2, double h1, double h2):start(start),ran(TRandom3 (0)), position(Point(ran, R1, R2, h1, h2)){
  
}

void Event::SetPosition(const Point& position){
  
  this->position = position;

}

void Event::SetInCylinder(double R, double h){//fills the position member with a random point in a cylinder of radius R and height h
  
  position = Point(ran, R, h);
  
}

void Event::SetInHollowCylinder(double R1, double R2, double h1, double h2){//fills the position member with a random point between two cylinders of respective radii R1 and R2 and of respective height h1 and h2
  
  position = Point(ran, R1, R2, h1, h2);
  
}

void Event::ChooseDecayRoad(){
 
  road.clear();
  stable.clear();
  
  State* old = start;
  State* current = old->Next(ran);
  road.push_back(old);
  
  while(!current->IsStable()){
     
    old = current;
    current = old->Next(ran);
    road.push_back(old);
    for(State* child : old->GetStableChildren()) stable.push_back(child);
    
  }
  
}

void Event::AdaptRoad(){

  for(State* state : road) state->AdaptMass(ran);//acts on RStates and propagates the mass to the children (VStates)

  while(!all_of(road.begin(), road.end(), [](State* s){return s->CanDecay();})) //VStates must be kinematically possible
    AdaptRoad();//otherwise reshuffle the road

}

void Event::DecayRoadStates(){
  
  for(State* state : road) state->Process(ran); //decay virtual states

}

void Event::Process(){
  
  ChooseDecayRoad();
  AdaptRoad();
  DecayRoadStates();
  
}

void Event::Fill(TH1F& h){

  for(State* s : stable) h.Fill(s->GetKinetic());
  
}

void Event::Fill(TH1F& h, int PDG){

  for(State* s : stable) if(s->GetPDG() == PDG) h.Fill(s->GetKinetic());
  
}

const Point& Event::GetPosition() const{
  
  return position;

}

const vector<State*>& Event::GetStableStates() const{
  
  return stable;

}

