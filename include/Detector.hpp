#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include "Volume.hpp"

class Detector {
  
public:
  
friend ostream& operator<<(ostream& output, const Detector& detector);

  vector<Volume*> volumes;
  Detector();
  void InsertVolume(const Volume* volume);

};

#endif