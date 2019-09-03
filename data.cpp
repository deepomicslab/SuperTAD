//
// Created by wang mengbo on 2019-09-01.
//

#include "data.h"


Data::Data (std::string fileName)
{
  _reader = new Reader (fileName);
  _N = _reader->parse (_contactMat);
  std::cout << "N=" << _N << std::endl;
  if (_K < 0)
    _K = _N / 5;
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
//      std::cout << "i=" << i << ", j=" << j << ", n=" << n << std::endl;
      Eigen::MatrixXd currentMat = _contactMat.block (i, i, n, n);
      int countIntra = (currentMat.sum () - currentMat.diagonal().sum ()) * .5;
      
      double sum1 = _contactMat.block (0, i, i, j-i+1).sum ();
//      std::cout << "sum1=" << sum1 << std::endl;
  
      double sum2;
      if (j + 1 < _N)
        sum2 = _contactMat.block (i, j+1, j-i+1, _N-1-j).sum ();
      else
        sum2 = 0;
//      std::cout << "sum2=" << sum2 << std::endl;
      int countInter = sum1 + sum2;
      _edgeCount(i, j) = countIntra;
      _edgeCount(j, i) = countInter;
    }
  }
//  std::cout << "_edgeCount:\n" << _edgeCount << std::endl;
}
