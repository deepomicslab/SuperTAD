//
// Created by wang mengbo on 2019-09-01.
//

#ifndef PROGRAM_DATA_H
#define PROGRAM_DATA_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <map>
#include "params.h"
#include "inputAndOutput.h"


class Data {
private:
  Eigen::MatrixXd _contactMat;
  
  // upper tri is intra; lower tri is inter
  Eigen::MatrixXd _edgeCount;
  
  Reader *_reader;
  
public:
  Data (std::string fileName);
  
  ~Data ();
  
  void init ();
  
  Eigen::MatrixXd &edgeCount () { return _edgeCount; }
  
  double getVol (int s, int e);
  
  double edgeSum () { return _edgeCount.coeff(0, _N-1); }
};

#endif //PROGRAM_DATA_H
