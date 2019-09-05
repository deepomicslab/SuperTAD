//
// Created by wang mengbo on 2019-09-03.
//

#ifndef PROGRAM_DETECTORMULTI_H
#define PROGRAM_DETECTORMULTI_H

#include <limits>
#include "params.h"
#include "data.h"
#include "inputAndOutput.h"
#include "multiTree.h"
#include <map>


bool cmpBoundary (const std::pair<int, int> &p1, const std::pair<int, int> &p2);

class DetectorMulti {
private:
  Data *_data;
  Eigen::MatrixXd *_edgeCount;
  Writer _writer;
  multi::Tree _multiTree;
  std::vector<multi::TreeNode *> *_nodeList;
  double *****_table;
  int *****_minIndexArray;
  int *****_leftKArray;
  std::vector<std::pair<int, int>> _boundary;
  std::map<int, int> _kToIdx;
  
public:
  DetectorMulti (Data &data);

  ~DetectorMulti ();
  
  int indexK (int k) { return _kToIdx.find(k)->second; }
  
  double getSE (int x, int y, double a, double b);
  
  double getSE (double x, double a, double b);
  
  void execute ();
  
  void fillTable ();
  
  void initH (int h);
  
  void backTrace (int k, int h, bool add=false);
  
  void multiSplit (int start, int end, int k, int h, int parentEnd, bool add=false);
};

#endif //PROGRAM_DETECTORMULTI_H
