//
// Created by wang mengbo on 2019-09-02.
//

#ifndef PROGRAM_UTILS_H
#define PROGRAM_UTILS_H

#include <cmath>
#include <utility>


namespace utils {
  
  int randInt (int low=0, int high=10);
  
  double randDouble (double low=0., double high=1.);
  
  typedef std::pair<int, int> boundary;
  
  inline bool cmpBoundary (const boundary &p1, const boundary &p2) { return p1.first < p2.first; }
  
  typedef std::pair<int, double> normLeaf;
  
  inline bool cmpNormLeaf (const normLeaf &p1, const normLeaf &p2) { return p1.second < p2.second; }
  
  bool doubleArrayEqual (double a1[], double a2[], int n);
  
  void copyDoubleArray (double from[], double to[], int n);
  
}

#endif //PROGRAM_UTILS_H
