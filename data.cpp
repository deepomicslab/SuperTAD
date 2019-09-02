//
// Created by wang mengbo on 2019-09-01.
//

#include "data.h"


namespace data {
  
  Reader::Reader (std::string fileName)
  {
    _fileName = fileName;
  }
  
  
  int Reader::parse (Eigen::MatrixXd &contactMat)
  {
    std::string line;
    double c;
    
    _infile.open (_fileName);
    bool init = false;
    int count = 0;
    int pos1 = 0;
    while (getline (_infile, line)) {
      std::istringstream iss (line);
      if (!init) {
        std::string lineBack = line;
        std::istringstream issBack (lineBack);
        while (issBack >> c) {
            count++;
        }
        contactMat.resize (count, count);
        init = true;
      }
      int pos2 = 0;
      while (iss >> c) {
        contactMat(pos1, pos2) = c;
        pos2++;
      }
      pos1++;
    }
    _infile.close ();
    
    return count;
  }
  
  
  Data::Data (std::string fileName)
  {
    _reader = new Reader (fileName);
    _N = _reader->parse (_contactMat);
  }
  
  
  Data::~Data ()
  {
    delete _reader;
  }
  
  
  void Data::init ()
  {
    _edgeCount.resize (_N, _N);
    for (int i = 0; i < _N; i++) {
      for (int j = i; j < _N; j++) {
        int n = j - i + 1;
        Eigen::MatrixXd currentMat = _contactMat.block (i, j, n, n);
        int countIntra = (currentMat.sum () - currentMat.diagonal().sum ()) * .5;
        int countInter = _contactMat.block (0, i, i, j-i+1).sum () + _contactMat.block (i, j+1, j-i+1, n-j).sum ();
        _edgeCount(i, j) = countIntra;
        _edgeCount(j, i) = countInter;
      }
    }
  }
}
