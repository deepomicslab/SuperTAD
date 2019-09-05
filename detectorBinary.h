//
// Created by wang mengbo on 2019-09-01.
//

#ifndef PROGRAM_DETECTORBINARY_H
#define PROGRAM_DETECTORBINARY_H

#include <limits>
#include <algorithm>
#include <set>
#include "params.h"
#include "data.h"
#include "binaryTree.h"
#include "inputAndOutput.h"
#include "utils.h"


typedef std::numeric_limits<double> infDouble;

bool cmpBoundary (const std::pair<int, int> &p1, const std::pair<int, int> &p2);

bool doubleArrayEqual (double a1[], double a2[], int n);

void copyDoubleArray (double from[], double to[], int n);

class DetectorBinary {
private:
  Data *_data;
  Eigen::MatrixXd *_edgeCount;
  Writer _writer;
  binary::Tree *_binaryTree;
  std::vector<binary::TreeNode *> *_nodeList;
  double ***_table;
  int ***_minIndexArray;
  int ***_leftKArray;
  std::vector<std::pair<int, int>> _boundary;
  std::set<binary::TreeNode *> _trueNodeList;
  std::map<int, int> _kToIdx;

public:
  DetectorBinary (Data &data);
  
  ~DetectorBinary ();
  
  void execute ();
  
  void init ();
  
  void fillTable ();
  
  int indexK (int k) { return _kToIdx.find(k)->second; };
  
  void backTrace (int k, bool add=false);
  
  void binarySplit (int start, int end, int k, bool add=false);
  
  void calculateD (binary::TreeNode &node);
  
  void calculateDensity (binary::TreeNode &node);
  
  double minusParent (double d, binary::TreeNode &node);
  
  void filterNodes ();
  
  double getX (binary::TreeNode &node);
  
  double getY (binary::TreeNode &node);
  
  bool simpleLinearRegression (std::vector<std::pair<int, binary::TreeNode *>> &nodeList, double ab[]);
};


#endif //PROGRAM_DETECTORBINARY_H
