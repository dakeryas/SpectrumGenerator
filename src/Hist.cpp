# include "Hist.hpp"

ostream& operator<<(ostream& output, const Hist& h){
  
  output<<h.getName()<<"\n";
  for(unsigned k = 0; k<h.getNumberOfBins(); ++k) output<<"["<<setw(6)<<internal<<h.getBinLowEdge(k)<<", "<<setw(6)<<internal<<h.getBinUpEdge(k)<<"]"<<setw(8)<<left<<" "<<"-->"<<setw(8)<<left<<" "<<setw(12)<<left<<h.getBinContent(k)<<"\n";
  return output;
  
}

Hist::Hist(const string& name, const vector<double>& edge, const vector<double>& value):name(name),edge(edge),value(value){
  
}

Hist::Hist(const string& name, const vector<double>& edge):Hist(name, edge, vector<double> (edge.size()-1, 0)){

}

Hist::Hist(const string& name, const unsigned numberOfBins, const double lowerEdge, const double upperEdge):Hist(name, vector<double>(numberOfBins+1), vector<double>(numberOfBins, 0)){

  buildConstantBining(numberOfBins, lowerEdge, upperEdge);
  
}

void Hist::buildConstantBining(const unsigned int numberOfBins, const double lowerEdge, const double upperEdge){
  
  if(numberOfBins == getNumberOfBins() && upperEdge - lowerEdge >0) edge = vector<double>(numberOfBins+1, (upperEdge - lowerEdge)/numberOfBins);

}


const string& Hist::getName() const{
  
  return name;
  
}

const double& Hist::getBinLowEdge(const unsigned i) const{
  
  return edge.at(i);

}

const double& Hist::getBinUpEdge(const unsigned i) const{
  
  return edge.at(i+1);

}

double Hist::getBinWidth(const unsigned int i) const{
  
  return edge.at(i+1) -edge.at(i);

}

unsigned Hist::getNumberOfBins() const{
  
  return value.size();

}

const double& Hist::getBinContent(const unsigned i) const{

  return value.at(i); 
  
}

unsigned Hist::getIndexOfBinContaining(const double x) const{

  return find_if(edge.begin()+1, edge.end(), [&x](const double upEdge){return x<upEdge;}) - 1 - edge.begin();
  
}

double Hist::getLowEdgeOfBinContaining(const double x) const{
  
  return *(find_if(edge.begin()+1, edge.end(), [&x](const double upEdge){return x<upEdge;}) - 1);

}

double Hist::getUpEdgeOfBinContaining(const double x) const{

  return *find_if(edge.begin()+1, edge.end(), [&x](const double upEdge){return x<upEdge;});
  
}

double Hist::getValueAt(const double x) const{
  
  unsigned lowIndex = getIndexOfBinContaining(x);
  return value[lowIndex];

}

double Hist::getIntegral() const{
  
  double integral = 0;
  for(unsigned k = 0; k<getNumberOfBins(); ++k) integral += value[k];
  return integral;
  
}

unsigned Hist::getIndexOfMaximum() const{
  
  return max_element(value.begin(), value.end()) - value.begin();

}

double Hist::getMaximumValue() const{
  
  return *max_element(value.begin(), value.end());

}

double Hist::getRangeWidth() const{

  return edge[getNumberOfBins()+1] - edge[0];
  
}

double Hist::getRandomFromRejection(TRandom3& ran) const{
  
  const double ymax = getMaximumValue();
  const double xmax = getRangeWidth();
  double y = ymax*ran.Rndm();
  double x = xmax*ran.Rndm();
  while(y > getValueAt(x)){
    
    y = ymax*ran.Rndm();
    x = xmax*ran.Rndm();
    
  }
  return x;

}

vector<double> Hist::getCumulativeIntegral() const{
  
  vector<double> cumulativeIntegral(getNumberOfBins() + 1);
  double currentIntegral = 0;

  for(unsigned k = 0; k<getNumberOfBins(); ++k){
    
    cumulativeIntegral[k] = currentIntegral;
    currentIntegral += value[k];
    
  }
  cumulativeIntegral[getNumberOfBins()] = currentIntegral;//add the last element which is the integral of all bins
  
  return cumulativeIntegral;
  
}

double Hist::getRandomFromTransform(TRandom3& ran) const{
  
  vector<double> cumulativeIntegral = getCumulativeIntegral();
  double y = cumulativeIntegral[getNumberOfBins()]*ran.Rndm();
  unsigned lowEdgeIndex = find_if(cumulativeIntegral.begin()+1, cumulativeIntegral.end(), [&y](const double currentIntegral){return y<currentIntegral;}) - 1 - cumulativeIntegral.begin();
  return edge[lowEdgeIndex] + (edge[lowEdgeIndex+1] - edge[lowEdgeIndex])*(y - cumulativeIntegral[lowEdgeIndex])/(cumulativeIntegral[lowEdgeIndex+1] - cumulativeIntegral[lowEdgeIndex]);

}

double Hist::getRandom(TRandom3& ran, randomType choice) const{
  
  if(choice == transform) return getRandomFromTransform(ran);
  else return getRandomFromRejection(ran);

}

void Hist::setBinContent(const double y, const unsigned int i){
  
  value.at(i) = y;

}


