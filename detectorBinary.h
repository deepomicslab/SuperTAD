//
// Created by wang mengbo on 2019-09-01.
//

#ifndef PROGRAM_DETECTORBINARY_H
#define PROGRAM_DETECTORBINARY_H

#include <limits>
#include <algorithm>
#include "params.h"
#include "data.h"
#include "binaryTree.h"
#include "inputAndOutput.h"
#include "utils.h"


typedef std::numeric_limits<double> infDouble;

bool cmpBoundary (const std::pair<int, int> &p1, const std::pair<int, int> &p2) {
  return p1.first < p2.first;
}

class DetectorBinary {
private:
  data::Data *_data;
  
  Eigen::MatrixXd *_edgeCount;
  
  Writer *_writer;
  
  BinaryTree *_binaryTree;
  
  std::vector<TreeNode *> *_nodeList;
  
  double ***_table;
  
  int ***_minIndexArray;
  
  int ***_leftKArray;
  
  std::vector<std::pair<int, int>> _boundary;
  
//  std::vector<double> _nodeSize;
  
  std::vector<std::pair<int, TreeNode *>> _trueNodes;
  
public:
  DetectorBinary (data::Data &data);
  
  ~DetectorBinary ();
  
  void execute ();
  
  void init ();
  
  void fillTable ();
  
  void backTrace (int k, bool add=false);
  
  void binarySplit (int start, int end, int k, bool add=false);
  
  void calculateD (TreeNode &node);
  
  void calculateDensity (TreeNode &node);
  
  double minusParent (double d, TreeNode &node);
  
  void filterNodes ();
  
  double getX (TreeNode &node);
  
  double getY (TreeNode &node);
  
  void simpleLinearRegression (std::vector<std::pair<int, TreeNode *>> &nodeList, double *ab);
};


#endif //PROGRAM_DETECTORBINARY_H
