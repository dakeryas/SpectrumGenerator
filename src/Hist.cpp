# include "Hist.hpp"

std::ostream& operator<<(std::ostream& output, const Hist& h){
  
  output<<h.getName()<<"\n";
  for(unsigned k = 0; k<h.getNumberOfBins(); ++k) output<<"["<<std::setw(6)<<std::internal<<h.getBinLowEdge(k)<<", "<<std::setw(6)<<std::internal<<h.getBinUpEdge(k)<<"]"<<std::setw(8)<<std::left<<" "<<"-->"<<std::setw(8)<<std::left<<" "<<std::setw(12)<<std::left<<h.getBinContent(k)<<"\n";
  return output;
  
}

Hist::Hist(const std::string& name, const std::vector<double>& edge, const std::vector<double>& value):name(name),edge(edge),value(value){
  
}

Hist::Hist(const std::string& name, const std::vector<double>& edge):Hist(name, edge, std::vector<double> (edge.size()-1, 0)){

}

Hist::Hist(const std::string& name, unsigned numberOfBins, double lowerEdge, double upperEdge):Hist(name, std::vector<double>(numberOfBins+1), std::vector<double>(numberOfBins, 0)){

  buildConstantBining(numberOfBins, lowerEdge, upperEdge);
  
}

void Hist::buildConstantBining(unsigned numberOfBins, double lowerEdge, double upperEdge){
  
  if(numberOfBins == getNumberOfBins() && upperEdge - lowerEdge >0) edge = std::vector<double>(numberOfBins+1, (upperEdge - lowerEdge)/numberOfBins);

}


const std::string& Hist::getName() const{
  
  return name;
  
}

const double& Hist::getBinLowEdge(unsigned i) const{
  
  return edge.at(i);

}

const double& Hist::getBinUpEdge(unsigned i) const{
  
  return edge.at(i+1);

}

double Hist::getBinWidth(unsigned i) const{
  
  return edge.at(i+1) -edge.at(i);

}

unsigned Hist::getNumberOfBins() const{
  
  return value.size();

}

const double& Hist::getBinContent(unsigned i) const{

  return value.at(i); 
  
}

unsigned Hist::getIndexOfBinContaining(double x) const{

  return std::find_if(edge.begin()+1, edge.end(), [&x](double upEdge){return x<upEdge;}) - 1 - edge.begin();
  
}

double Hist::getLowEdgeOfBinContaining(double x) const{
  
  return *(std::find_if(edge.begin()+1, edge.end(), [&x](double upEdge){return x<upEdge;}) - 1);

}

double Hist::getUpEdgeOfBinContaining(double x) const{

  return *std::find_if(edge.begin()+1, edge.end(), [&x](double upEdge){return x<upEdge;});
  
}

double Hist::getValueAt(double x) const{
  
  unsigned lowIndex = getIndexOfBinContaining(x);
  return value[lowIndex];

}

double Hist::getIntegral() const{
  
  double integral = 0;
  for(unsigned k = 0; k<getNumberOfBins(); ++k) integral += value[k];
  return integral;
  
}

unsigned Hist::getIndexOfMaximum() const{
  
  return std::max_element(value.begin(), value.end()) - value.begin();

}

double Hist::getMaximumValue() const{
  
  return *std::max_element(value.begin(), value.end());

}

double Hist::getRangeWidth() const{

  return edge[getNumberOfBins()+1] - edge[0];
  
}

double Hist::getRandomFromRejection(TRandom3& ran) const{
  
  double ymax = getMaximumValue();
  double xmax = getRangeWidth();
  double y = ymax*ran.Rndm();
  double x = xmax*ran.Rndm();
  while(y > getValueAt(x)){
    
    y = ymax*ran.Rndm();
    x = xmax*ran.Rndm();
    
  }
  return x;

}

std::vector<double> Hist::getCumulativeIntegral() const{
  
  std::vector<double> cumulativeIntegral(getNumberOfBins() + 1);
  double currentIntegral = 0;

  for(unsigned k = 0; k<getNumberOfBins(); ++k){
    
    cumulativeIntegral[k] = currentIntegral;
    currentIntegral += value[k];
    
  }
  cumulativeIntegral[getNumberOfBins()] = currentIntegral;//add the last element which is the integral of all bins
  
  return cumulativeIntegral;
  
}

double Hist::getRandomFromTransform(TRandom3& ran) const{
  
  std::vector<double> cumulativeIntegral = getCumulativeIntegral();
  double y = cumulativeIntegral[getNumberOfBins()]*ran.Rndm();
  unsigned lowEdgeIndex = std::find_if(cumulativeIntegral.begin()+1, cumulativeIntegral.end(), [&y](double currentIntegral){return y<currentIntegral;}) - 1 - cumulativeIntegral.begin();
  return edge[lowEdgeIndex] + (edge[lowEdgeIndex+1] - edge[lowEdgeIndex])*(y - cumulativeIntegral[lowEdgeIndex])/(cumulativeIntegral[lowEdgeIndex+1] - cumulativeIntegral[lowEdgeIndex]);

}

double Hist::getRandom(TRandom3& ran, randomType choice) const{
  
  if(choice == transform) return getRandomFromTransform(ran);
  else return getRandomFromRejection(ran);

}

void Hist::setBinContent(double y, unsigned i){
  
  value.at(i) = y;

}


