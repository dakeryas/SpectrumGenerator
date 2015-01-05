#ifndef HIST_H
#define HIST_H

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <TRandom3.h>

using namespace::std;

class Hist { 
  
  string name;
  vector<double> edge;
  vector<double> value;
  enum randomType{rejection, transform};
  void buildConstantBining(const unsigned numberOfBins, const double lowerEdge, const double upperEdge);
  double getRandomFromRejection(TRandom3& ran) const;
  double getRandomFromTransform(TRandom3& ran) const;
  vector<double> getCumulativeIntegral() const;
  
public:
  Hist(const string& name, const vector<double>& edge);
  Hist(const string& name, const vector<double>& edge, const vector<double>& value);
  Hist(const string& name, const unsigned numberOfBins, const double lowerEdge, const double upperEdge);
  const string& getName() const;
  const double& getBinLowEdge(const unsigned i) const;
  const double& getBinUpEdge(const unsigned i) const;
  double getBinWidth(const unsigned i) const; 
  unsigned getNumberOfBins() const;
  const double& getBinContent(const unsigned i) const;
  unsigned getIndexOfBinContaining(const double x) const;
  double getLowEdgeOfBinContaining(const double x) const;
  double getUpEdgeOfBinContaining(const double x) const;
  double getValueAt(const double x) const;
  double getIntegral() const;
  unsigned getIndexOfMaximum() const;
  double getMaximumValue() const;
  double getRangeWidth() const;
  double getRandom(TRandom3& ran, randomType choice = transform) const;
  void setBinContent(const double y, const unsigned i);

};

ostream& operator<<(ostream& output, const Hist& h); 

#endif