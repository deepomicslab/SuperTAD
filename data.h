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


namespace data {
  class Reader {
  private:
    std::ifstream _infile;
    std::string _fileName;
  
  public:
    Reader (std::string fileName);
    
    ~Reader () {};
    
    int parse (Eigen::MatrixXd &contactMat);
    
  };
  
  
  class Edges {
  private:
    std::map<int, std::map<int, double>> _contacts;
  public:
    Edges ();
    
    ~Edges ();
    
    void init ();
  };
  
  
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
    
    double edgeSum () { return _edgeCount.coeff(0, _N-1); }
  };
  
}


#endif //PROGRAM_DATA_H
