#ifndef HIST2D_H
#define HIST2D_H

#include <cmath>
#include <TH2F.h>
#include "Hist.hpp"

class Hist2d { 
    
  string name;
  vector<double> xedge;
  vector<double> yedge;
  vector<double> value;
  void buildConstantBiningX(const unsigned numberOfBinsX, const double lowerEdgeX, const double upperEdgeX);
  void buildConstantBiningY(const unsigned numberOfBinsY, const double lowerEdgeY, const double upperEdgeY);
  const vector<double>& getXedge() const;
  const vector<double>& getYedge() const;
  
public:
  Hist2d();
  Hist2d(const string& name, const vector<double>& xedge, const vector<double>& yedge);
  Hist2d(const string& name, const vector<double>& xedge, const vector<double>& yedge, const vector<double>& value);
  Hist2d(const string& name, const unsigned numberOfBinsX, const double lowerEdgeX, const double upperEdgeX, const unsigned numberOfBinsY, const double lowerEdgeY, const double upperEdgeY);
  Hist2d(const TH2F* hROOT);//copy the contents of the ROOT histogram
  const string& getName() const;
  const double& getBinXLowEdge(const unsigned i) const;
  double getBinXCenter(const unsigned i) const;
  const double& getBinXUpEdge(const unsigned i) const;
  double getBinXWidth(const unsigned i) const; 
  const double& getBinYLowEdge(const unsigned i) const;
  double getBinYCenter(const unsigned i) const;
  const double& getBinYUpEdge(const unsigned i) const;
  double getBinYWidth(const unsigned i) const; 
  unsigned getNumberOfBinsX() const;
  unsigned getNumberOfBinsY() const;
  unsigned getNumberOfBins() const;
  const double& getBinContent(const unsigned i, const unsigned j) const;
  unsigned getIndexOfBinXContaining(const double x) const;
  double getLowEdgeOfBinXContaining(const double x) const;
  double getUpEdgeOfBinXContaining(const double x) const;
  unsigned getIndexOfBinYContaining(const double y) const;
  double getLowEdgeOfBinYContaining(const double y) const;
  double getUpEdgeOfBinYContaining(const double y) const;
  double getValueAt(const double x, const double y) const;
  double getIntegral() const;
  pair<unsigned, unsigned> getIndicesOfMaximum() const;
  double getMaximumValue() const;
  double getXRangeWidth() const;
  double getYRangeWidth() const;
  void getInterpolatedSliceAtX(const double x, Hist& interpolated) const; //the result is stored in interpolated and assumes that interpolated.getNumberOfBins() = getNumberOfBinsY()
  void getInterpolatedSliceAtY(const double y, Hist& interpolated) const;
  double getRandomFromInterpolatedSliceAtX(TRandom3& ran, const double x) const;
  double getRandomFromInterpolatedSliceAtY(TRandom3& ran, const double y) const;
  void setBinContent(const double z, const unsigned i, const unsigned j);

};

ostream& operator<<(ostream& output, const Hist2d& h); 

#endif