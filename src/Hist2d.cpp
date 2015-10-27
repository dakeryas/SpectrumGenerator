# include "Hist2d.hpp"

std::ostream& operator<<(std::ostream& output, const Hist2d& h){
  
  output<<h.getName()<<"\n";
  for(unsigned i = 0; i<h.getNumberOfBinsX(); ++i)
    for(unsigned j = 0; j<h.getNumberOfBinsY(); ++j)
      output<<"["<<std::setw(6)<<std::internal<<h.getBinXLowEdge(i)<<", "<<std::setw(6)<<std::internal<<h.getBinXUpEdge(i)<<"]"
       <<"  x  "<<"["<<std::setw(6)<<std::internal<<h.getBinYLowEdge(j)<<", "<<std::setw(6)<<std::internal<<h.getBinYUpEdge(j)<<"]"
      <<std::setw(8)<<std::left<<" "<<"-->"<<std::setw(8)<<std::left<<" "<<std::setw(12)<<std::left<<h.getBinContent(i, j)<<"\n";
  return output;
  
}

Hist2d::Hist2d():name("dummy"){

}

Hist2d::Hist2d(const std::string& name, const std::vector<double>& xedge, const std::vector<double>& yedge, const std::vector<double>& value):name(name),xedge(xedge),yedge(yedge),value(value){
  
}

Hist2d::Hist2d(const std::string& name, const std::vector<double>& xedge, const std::vector<double>& yedge):Hist2d(name, xedge, yedge, std::vector<double> ((xedge.size()-1)*(yedge.size()-1), 0)){

}

Hist2d::Hist2d(const std::string& name, unsigned numberOfBinsX, double lowerEdgeX, double upperEdgeX, unsigned numberOfBinsY, double lowerEdgeY, double upperEdgeY):Hist2d(name, std::vector<double>(numberOfBinsX+1), std::vector<double>(numberOfBinsY+1, 0), std::vector<double>(numberOfBinsX*numberOfBinsY, 0)){

  buildConstantBiningX(numberOfBinsX, lowerEdgeX, upperEdgeX);
  buildConstantBiningY(numberOfBinsY, lowerEdgeY, upperEdgeY);
  
}

Hist2d::Hist2d(const TH2F* hROOT):name(hROOT->GetName()),xedge(hROOT->GetNbinsX()+1), yedge(hROOT->GetNbinsY()+1),value((hROOT->GetNbinsX())*(hROOT->GetNbinsY())){

  double xarray[xedge.size() - 1];
  hROOT->GetXaxis()->GetLowEdge(xarray);
  xedge = std::vector<double>(xarray, xarray + xedge.size() - 1);
  xedge.back() = hROOT->GetXaxis()->GetBinLowEdge(xedge.size());//the low edge of the overflow bin should be the upper edge of the previous bin
  
  double yarray[yedge.size() - 1];
  hROOT->GetYaxis()->GetLowEdge(yarray);
  yedge = std::vector<double>(yarray, yarray + yedge.size() - 1);
  yedge.back() = hROOT->GetYaxis()->GetBinLowEdge(yedge.size());
  
  for(unsigned i = 0; i<getNumberOfBinsX(); ++i)
    for(unsigned j = 0; j<getNumberOfBinsY(); ++j)
      setBinContent(hROOT->GetBinContent(i+1, j+1), i, j);
  
}

void Hist2d::buildConstantBiningX(unsigned numberOfBinsX, double lowerEdgeX, double upperEdgeX){
  
  if(numberOfBinsX == getNumberOfBins() && upperEdgeX - lowerEdgeX >0) xedge = std::vector<double>(numberOfBinsX+1, (upperEdgeX - lowerEdgeX)/numberOfBinsX);

}

void Hist2d::buildConstantBiningY(unsigned numberOfBinsY, double lowerEdgeY, double upperEdgeY){
  
  if(numberOfBinsY == getNumberOfBins() && upperEdgeY - lowerEdgeY >0) yedge = std::vector<double>(numberOfBinsY+1, (upperEdgeY - lowerEdgeY)/numberOfBinsY);

}

const std::vector<double>& Hist2d::getXedge() const{

  return xedge;
  
}

const std::vector<double>& Hist2d::getYedge() const{

  return yedge;
  
}

const std::string& Hist2d::getName() const{
  
  return name;
  
}

const double& Hist2d::getBinXLowEdge(unsigned i) const{
  
  return xedge.at(i);

}

double Hist2d::getBinXCenter(unsigned int i) const{
  
  return 0.5*(xedge.at(i)+xedge.at(i+1));

}

const double& Hist2d::getBinXUpEdge(unsigned i) const{
  
  return xedge.at(i+1);

}

double Hist2d::getBinXWidth(unsigned int i) const{
  
  return xedge.at(i+1) -xedge.at(i);

}

const double& Hist2d::getBinYLowEdge(unsigned i) const{
  
  return yedge.at(i);

}

double Hist2d::getBinYCenter(unsigned int i) const{
  
  return 0.5*(yedge.at(i)+yedge.at(i+1));

}

const double& Hist2d::getBinYUpEdge(unsigned i) const{
  
  return yedge.at(i+1);

}

double Hist2d::getBinYWidth(unsigned int i) const{
  
  return yedge.at(i+1) -yedge.at(i);

}

unsigned Hist2d::getNumberOfBinsX() const{
  
  return xedge.size()-1;

}

unsigned Hist2d::getNumberOfBinsY() const{
  
  return yedge.size()-1;

}

unsigned Hist2d::getNumberOfBins() const{
  
  return value.size();

}

const double& Hist2d::getBinContent(unsigned i, unsigned j) const{

  return value.at(i*getNumberOfBinsY() + j); 
  
}

unsigned Hist2d::getIndexOfBinXContaining(double x) const{

  return std::find_if(xedge.begin()+1, xedge.end(), [&x](double upEdge){return x<upEdge;}) - 1 - xedge.begin();
  
}

double Hist2d::getLowEdgeOfBinXContaining(double x) const{
  
  return *(std::find_if(xedge.begin()+1, xedge.end(), [&x](double upEdge){return x<upEdge;}) - 1);

}

double Hist2d::getUpEdgeOfBinXContaining(double x) const{

  return *std::find_if(xedge.begin()+1, xedge.end(), [&x](double upEdge){return x<upEdge;});
  
}

unsigned Hist2d::getIndexOfBinYContaining(double y) const{

  return std::find_if(yedge.begin()+1, yedge.end(), [&y](double upEdge){return y<upEdge;}) - 1 - yedge.begin();
  
}

double Hist2d::getLowEdgeOfBinYContaining(double y) const{
  
  return *(std::find_if(yedge.begin()+1, yedge.end(), [&y](double upEdge){return y<upEdge;}) - 1);

}

double Hist2d::getUpEdgeOfBinYContaining(double y) const{

  return *std::find_if(yedge.begin()+1, yedge.end(), [&y](double upEdge){return y<upEdge;});
  
}

double Hist2d::getValueAt(double x, double y) const{
  
  unsigned lowXIndex = getIndexOfBinXContaining(x);
  unsigned lowYIndex = getIndexOfBinYContaining(y);
  return value[lowXIndex*getNumberOfBinsY() + lowYIndex];

}

double Hist2d::getIntegral() const{
  
  double integral = 0;
  for(unsigned k = 0; k<getNumberOfBins(); ++k) integral += value[k];
  return integral;
  
}

std::pair<unsigned, unsigned> Hist2d::getIndicesOfMaximum() const{
  
  std::pair<unsigned, unsigned> indices;
  unsigned linearIndex =  std::max_element(value.begin(), value.end()) - value.begin();
  indices.first = linearIndex/getNumberOfBinsY();
  indices.second = linearIndex - indices.first*getNumberOfBinsY();
  return indices;

}

double Hist2d::getMaximumValue() const{
  
  return *std::max_element(value.begin(), value.end());

}

double Hist2d::getXRangeWidth() const{

  return xedge[getNumberOfBinsX()+1] - xedge[0];
  
}

double Hist2d::getYRangeWidth() const{
  
  return yedge[getNumberOfBinsY()+1] - yedge[0];

}

void Hist2d::getInterpolatedSliceAtX(double x, Hist& interpolated) const{
    
  unsigned index = getIndexOfBinXContaining(x);
  
  if(x > getBinXCenter(index) && index != getNumberOfBinsX() - 1){
    
    double correction = std::abs(x-getBinXCenter(index))/(getBinXCenter(index+1) - getBinXCenter(index));
    for(unsigned j = 0;j<getNumberOfBinsY(); ++j)
      interpolated.setBinContent(getBinContent(index,j)+correction*(getBinContent(index+1,j) - getBinContent(index,j)), j);
    
  }
  else if(x < getBinXCenter(index) && index != 0){
    
    double correction = std::abs(x-getBinXCenter(index))/(getBinXCenter(index) - getBinXCenter(index-1));
    for(unsigned j = 0;j<getNumberOfBinsY(); ++j)
      interpolated.setBinContent(getBinContent(index,j)+correction*(getBinContent(index-1,j) - getBinContent(index,j)), j);
    
  }
  else {//if we only have the first or the last bin, we cannot interpolate
    
    for(unsigned j = 0;j<getNumberOfBinsY(); ++j)
      interpolated.setBinContent(getBinContent(index,j), j);
    
  }
  
}

void Hist2d::getInterpolatedSliceAtY(double y, Hist& interpolated) const{
    
  unsigned index = getIndexOfBinYContaining(y);
  
  if(y > getBinYCenter(index) && index != getNumberOfBinsY() - 1){
    
    double correction = std::abs(y-getBinYCenter(index))/(getBinYCenter(index+1) - getBinYCenter(index));
    for(unsigned i = 0;i<getNumberOfBinsX(); ++i)
      interpolated.setBinContent(getBinContent(i,index)+correction*(getBinContent(i,index+1) - getBinContent(i,index)), i);
    
  }
  else if(y < getBinYCenter(index) && index != 0){
    
    double correction = std::abs(y-getBinYCenter(index))/(getBinYCenter(index) - getBinYCenter(index-1));
    for(unsigned i = 0;i<getNumberOfBinsX(); ++i)
      interpolated.setBinContent(getBinContent(i,index)+correction*(getBinContent(i,index-1) - getBinContent(i,index)), i);
    
  }
  else {//if we only have the first or the last bin, we cannot interpolate
    
    for(unsigned i = 0;i<getNumberOfBinsX(); ++i)
      interpolated.setBinContent(getBinContent(i,index), i);
    
  }
  
}

double Hist2d::getRandomFromInterpolatedSliceAtX(TRandom3& ran, double x) const{
  
  Hist interpolated(name+"_sliced_at"+std::to_string(x), getYedge());
  getInterpolatedSliceAtX(x, interpolated);
  return interpolated.getRandom(ran);

}

double Hist2d::getRandomFromInterpolatedSliceAtY(TRandom3& ran, double y) const{
  
  Hist interpolated(name+"_sliced_at"+std::to_string(y), getXedge());
  getInterpolatedSliceAtY(y, interpolated);
  return interpolated.getRandom(ran);

}

void Hist2d::setBinContent(double z, unsigned int i, unsigned j){
  
  value.at(i*getNumberOfBinsY() + j) = z;

}


