# include "Hist2d.hpp"

ostream& operator<<(ostream& output, const Hist2d& h){
  
  output<<h.getName()<<"\n";
  for(unsigned i = 0; i<h.getNumberOfBinsX(); ++i)
    for(unsigned j = 0; j<h.getNumberOfBinsY(); ++j)
      output<<"["<<setw(6)<<internal<<h.getBinXLowEdge(i)<<", "<<setw(6)<<internal<<h.getBinXUpEdge(i)<<"]"
       <<"  x  "<<"["<<setw(6)<<internal<<h.getBinYLowEdge(j)<<", "<<setw(6)<<internal<<h.getBinYUpEdge(j)<<"]"
      <<setw(8)<<left<<" "<<"-->"<<setw(8)<<left<<" "<<setw(12)<<left<<h.getBinContent(i, j)<<"\n";
  return output;
  
}

Hist2d::Hist2d():name("dummy"){

}

Hist2d::Hist2d(const string& name, const vector<double>& xedge, const vector<double>& yedge, const vector<double>& value):name(name),xedge(xedge),yedge(yedge),value(value){
  
}

Hist2d::Hist2d(const string& name, const vector<double>& xedge, const vector<double>& yedge):Hist2d(name, xedge, yedge, vector<double> ((xedge.size()-1)*(yedge.size()-1), 0)){

}

Hist2d::Hist2d(const string& name, const unsigned numberOfBinsX, const double lowerEdgeX, const double upperEdgeX, const unsigned numberOfBinsY, const double lowerEdgeY, const double upperEdgeY):Hist2d(name, vector<double>(numberOfBinsX+1), vector<double>(numberOfBinsY+1, 0), vector<double>(numberOfBinsX*numberOfBinsY, 0)){

  buildConstantBiningX(numberOfBinsX, lowerEdgeX, upperEdgeX);
  buildConstantBiningY(numberOfBinsY, lowerEdgeY, upperEdgeY);
  
}

Hist2d::Hist2d(const TH2F* hROOT):name(hROOT->GetName()),xedge(hROOT->GetNbinsX()+1), yedge(hROOT->GetNbinsY()+1),value((hROOT->GetNbinsX())*(hROOT->GetNbinsY())){

  double xarray[xedge.size() - 1];
  hROOT->GetXaxis()->GetLowEdge(xarray);
  xedge = vector<double>(xarray, xarray + xedge.size() - 1);
  xedge.back() = hROOT->GetXaxis()->GetBinLowEdge(xedge.size());//the low edge of the overflow bin should be the upper edge of the previous bin
  
  double yarray[yedge.size() - 1];
  hROOT->GetYaxis()->GetLowEdge(yarray);
  yedge = vector<double>(yarray, yarray + yedge.size() - 1);
  yedge.back() = hROOT->GetYaxis()->GetBinLowEdge(yedge.size());
  
  for(unsigned i = 0; i<getNumberOfBinsX(); ++i)
    for(unsigned j = 0; j<getNumberOfBinsY(); ++j)
      setBinContent(hROOT->GetBinContent(i+1, j+1), i, j);
  
}

void Hist2d::buildConstantBiningX(const unsigned int numberOfBinsX, const double lowerEdgeX, const double upperEdgeX){
  
  if(numberOfBinsX == getNumberOfBins() && upperEdgeX - lowerEdgeX >0) xedge = vector<double>(numberOfBinsX+1, (upperEdgeX - lowerEdgeX)/numberOfBinsX);

}

void Hist2d::buildConstantBiningY(const unsigned int numberOfBinsY, const double lowerEdgeY, const double upperEdgeY){
  
  if(numberOfBinsY == getNumberOfBins() && upperEdgeY - lowerEdgeY >0) yedge = vector<double>(numberOfBinsY+1, (upperEdgeY - lowerEdgeY)/numberOfBinsY);

}

const vector<double>& Hist2d::getXedge() const{

  return xedge;
  
}

const vector<double>& Hist2d::getYedge() const{

  return yedge;
  
}

const string& Hist2d::getName() const{
  
  return name;
  
}

const double& Hist2d::getBinXLowEdge(const unsigned i) const{
  
  return xedge.at(i);

}

double Hist2d::getBinXCenter(const unsigned int i) const{
  
  return 0.5*(xedge.at(i)+xedge.at(i+1));

}

const double& Hist2d::getBinXUpEdge(const unsigned i) const{
  
  return xedge.at(i+1);

}

double Hist2d::getBinXWidth(const unsigned int i) const{
  
  return xedge.at(i+1) -xedge.at(i);

}

const double& Hist2d::getBinYLowEdge(const unsigned i) const{
  
  return yedge.at(i);

}

double Hist2d::getBinYCenter(const unsigned int i) const{
  
  return 0.5*(yedge.at(i)+yedge.at(i+1));

}

const double& Hist2d::getBinYUpEdge(const unsigned i) const{
  
  return yedge.at(i+1);

}

double Hist2d::getBinYWidth(const unsigned int i) const{
  
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

const double& Hist2d::getBinContent(const unsigned i, const unsigned j) const{

  return value.at(i*getNumberOfBinsY() + j); 
  
}

unsigned Hist2d::getIndexOfBinXContaining(const double x) const{

  return find_if(xedge.begin()+1, xedge.end(), [&x](const double upEdge){return x<upEdge;}) - 1 - xedge.begin();
  
}

double Hist2d::getLowEdgeOfBinXContaining(const double x) const{
  
  return *(find_if(xedge.begin()+1, xedge.end(), [&x](const double upEdge){return x<upEdge;}) - 1);

}

double Hist2d::getUpEdgeOfBinXContaining(const double x) const{

  return *find_if(xedge.begin()+1, xedge.end(), [&x](const double upEdge){return x<upEdge;});
  
}

unsigned Hist2d::getIndexOfBinYContaining(const double y) const{

  return find_if(yedge.begin()+1, yedge.end(), [&y](const double upEdge){return y<upEdge;}) - 1 - yedge.begin();
  
}

double Hist2d::getLowEdgeOfBinYContaining(const double y) const{
  
  return *(find_if(yedge.begin()+1, yedge.end(), [&y](const double upEdge){return y<upEdge;}) - 1);

}

double Hist2d::getUpEdgeOfBinYContaining(const double y) const{

  return *find_if(yedge.begin()+1, yedge.end(), [&y](const double upEdge){return y<upEdge;});
  
}

double Hist2d::getValueAt(const double x, const double y) const{
  
  unsigned lowXIndex = getIndexOfBinXContaining(x);
  unsigned lowYIndex = getIndexOfBinYContaining(y);
  return value[lowXIndex*getNumberOfBinsY() + lowYIndex];

}

double Hist2d::getIntegral() const{
  
  double integral = 0;
  for(unsigned k = 0; k<getNumberOfBins(); ++k) integral += value[k];
  return integral;
  
}

pair<unsigned, unsigned> Hist2d::getIndicesOfMaximum() const{
  
  pair<unsigned, unsigned> indices;
  unsigned linearIndex =  max_element(value.begin(), value.end()) - value.begin();
  indices.first = linearIndex/getNumberOfBinsY();
  indices.second = linearIndex - indices.first*getNumberOfBinsY();
  return indices;

}

double Hist2d::getMaximumValue() const{
  
  return *max_element(value.begin(), value.end());

}

double Hist2d::getXRangeWidth() const{

  return xedge[getNumberOfBinsX()+1] - xedge[0];
  
}

double Hist2d::getYRangeWidth() const{
  
  return yedge[getNumberOfBinsY()+1] - yedge[0];

}

void Hist2d::getInterpolatedSliceAtX(const double x, Hist& interpolated) const{
    
  unsigned index = getIndexOfBinXContaining(x);
  
  if(x > getBinXCenter(index) && index != getNumberOfBinsX() - 1){
    
    double correction = abs(x-getBinXCenter(index))/(getBinXCenter(index+1) - getBinXCenter(index));
    for(unsigned j = 0;j<getNumberOfBinsY(); ++j)
      interpolated.setBinContent(getBinContent(index,j)+correction*(getBinContent(index+1,j) - getBinContent(index,j)), j);
    
  }
  else if(x < getBinXCenter(index) && index != 0){
    
    double correction = abs(x-getBinXCenter(index))/(getBinXCenter(index) - getBinXCenter(index-1));
    for(unsigned j = 0;j<getNumberOfBinsY(); ++j)
      interpolated.setBinContent(getBinContent(index,j)+correction*(getBinContent(index-1,j) - getBinContent(index,j)), j);
    
  }
  else {//if we only have the first or the last bin, we cannot interpolate
    
    for(unsigned j = 0;j<getNumberOfBinsY(); ++j)
      interpolated.setBinContent(getBinContent(index,j), j);
    
  }
  
}

void Hist2d::getInterpolatedSliceAtY(const double y, Hist& interpolated) const{
    
  unsigned index = getIndexOfBinYContaining(y);
  
  if(y > getBinYCenter(index) && index != getNumberOfBinsY() - 1){
    
    double correction = abs(y-getBinYCenter(index))/(getBinYCenter(index+1) - getBinYCenter(index));
    for(unsigned i = 0;i<getNumberOfBinsX(); ++i)
      interpolated.setBinContent(getBinContent(i,index)+correction*(getBinContent(i,index+1) - getBinContent(i,index)), i);
    
  }
  else if(y < getBinYCenter(index) && index != 0){
    
    double correction = abs(y-getBinYCenter(index))/(getBinYCenter(index) - getBinYCenter(index-1));
    for(unsigned i = 0;i<getNumberOfBinsX(); ++i)
      interpolated.setBinContent(getBinContent(i,index)+correction*(getBinContent(i,index-1) - getBinContent(i,index)), i);
    
  }
  else {//if we only have the first or the last bin, we cannot interpolate
    
    for(unsigned i = 0;i<getNumberOfBinsX(); ++i)
      interpolated.setBinContent(getBinContent(i,index), i);
    
  }
  
}

double Hist2d::getRandomFromInterpolatedSliceAtX(TRandom3& ran, const double x) const{
  
  Hist interpolated(name+"_sliced_at"+to_string(x), getYedge());
  getInterpolatedSliceAtX(x, interpolated);
  return interpolated.getRandom(ran);

}

double Hist2d::getRandomFromInterpolatedSliceAtY(TRandom3& ran, const double y) const{
  
  Hist interpolated(name+"_sliced_at"+to_string(y), getXedge());
  getInterpolatedSliceAtY(y, interpolated);
  return interpolated.getRandom(ran);

}

void Hist2d::setBinContent(const double z, const unsigned int i, const unsigned j){
  
  value.at(i*getNumberOfBinsY() + j) = z;

}


