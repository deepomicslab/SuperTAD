//
// Created by wang mengbo on 2019-09-01.
//

#include "detectorBinary.h"


DetectorBinary::DetectorBinary (data::Data &data)
{
  _data = &data;
  _edgeCount = &data.edgeCount ();
  _table = new double **[_N];
  _minIndexArray = new int **[_N];
  _leftKArray = new int **[_N];
  for (int i = 0; i < _N; i++) {
    _table[i] = new double *[_N];
    _minIndexArray[i] = new int *[_N];
    _leftKArray[i] = new int *[_N];
    for (int j = 0; j < _N; j++) {
      _table[i][j] = new double [_K];
      _minIndexArray[i][j] = new int [_K];
      _leftKArray[i][j] = new int [_K];
    }
  }
  _boundary.reserve (_K);
  _binaryTree = new BinaryTree ();
  
}


DetectorBinary::~DetectorBinary ()
{

}


void DetectorBinary::execute ()
{
  fillTable ();
  
  // determine K
  std::vector <double> sumOfEntropy;
  std::vector <double> sumOfLeaves;
  std::vector <std::pair<int, double>> normLeaves;
  struct cmpNormLeaf {
    bool operator ()(const std::pair<int, double> &p1, const std::pair<int, double> &p2) {
      return p1.second < p2.second;
    }
  };
  
  for (int num = 1; num < _K; num++) {
    sumOfEntropy.emplace_back (_table[0][_N-1][num]);
    init ();
    backTrace (num);
    double leafSum = 0;
    for (int leaf=0; leaf < _boundary.size (); leaf++) {
      int currentStart = _boundary[leaf].first;
      int currentEnd = _boundary[leaf].second;
      double currentVolumn;
      if (currentStart == currentEnd)
        currentVolumn = 2 * _edgeCount[currentStart][currentEnd] + _edgeCount[currentEnd][currentStart];
      else
        currentVolumn = _edgeCount[currentEnd][currentStart];
      leafSum += _edgeCount[currentEnd][currentStart] / (2 * _data->edgeSum ()) * log2(2 * _data->edgeSum () / currentVolumn);
      leafSum += _table[currentStart][currentEnd][0];
    }
    sumOfLeaves.emplace_back (leafSum);
    normLeaves.emplace_back (num, leafSum / (log2(_N / num) + _N * (num - 1) / (num * (_N - 1)) * log2(num)));
  }
  
  sort(normLeaves.begin (), normLeaves.end (), cmpNormLeaf ());
  int index = normLeaves[0].first;
  backTrace (index, true);
  
  _nodeList = &_binaryTree->nodeList ();
  _writer->writeTree (_WORK_DIR, "original_boundaries.txt", *_nodeList);
  
  // filtering
  if (_FILTERING) {
    calculateD (_binaryTree->root ());
    calculateDensity (_binaryTree->root ());
    
  }
  
}


void DetectorBinary::init ()
{
  _boundary.clear ();
}


void DetectorBinary::fillTable ()
{
  for (int start = 0; start < _N; start++) {
    for (int end = start; end < _N; end++) {
      double currentVolumn;
      if (start != end)
        currentVolumn = 2 * _data->edgeCount ().coeff (start, end) + _data->edgeCount ().coeff (end, start);
      else
        currentVolumn = _data->edgeCount ().coeff (end, start);
      for (int leaf = start; leaf < end + 1; leaf++) {
        double leafDegree = _data->edgeCount ().coeff (leaf, leaf);
        if (leafDegree != 0) {
          _table[start][end][0] += (leafDegree / (2. * _data->edgeSum ())) * log2(currentVolumn / leafDegree);
        }
      }
    }
  }
  
  for (int a = 1; a < _K; a++) {
    for (int start = 0; start < _N; start++) {
      for (int end = start; end < _N; end++) {
        double minTmp = infDouble::infinity ();
        int minIdx = 0;
        int leftK = 0;
        for (int binaryK = 1; binaryK < a; binaryK++) {
          for (int mid = start; mid < end; mid++) {
            double tmp = _table[start][mid][binaryK] + _table[mid + 1][end][a - binaryK];
            double volumnParent, currentVolumn1, currentVolumn2;
            if (start != end)
              volumnParent = 2 * _edgeCount->coeff(start, end) + _edgeCount->coeff(end, start);
            else
              volumnParent = _edgeCount->coeff(end, start);
            
            if (start != mid)
              currentVolumn1 = 2 * _edgeCount->coeff(start, mid) + _edgeCount->coeff(mid, start);
            else
              currentVolumn1 = _edgeCount->coeff(mid, start);
            
            if (mid +1 != end)
              currentVolumn2 = 2 * _edgeCount->coeff(mid + 1, end) + _edgeCount->coeff(end, mid+1);
            else
              currentVolumn2 = _edgeCount->coeff(end, mid+1);
            
            if (currentVolumn1 != 0)
              tmp += _edgeCount->coeff(mid, start) / (2. * _data->edgeSum ()) * log2(volumnParent / currentVolumn1);
            if (currentVolumn2 != 0)
              tmp += _edgeCount->coeff(end, mid+1) / (2. * _data->edgeSum ()) * log2(volumnParent / currentVolumn2);
            if (tmp < minTmp) {
              minTmp = tmp;
              minIdx = mid;
              leftK = binaryK;
            }
          }
        }
        _minIndexArray[start][end][a] = minIdx;
        _table[start][end][a] = minTmp;
        _leftKArray[start][end][a] = leftK;
      }
    }
  }
}


void DetectorBinary::backTrace (int k, bool add)
{
  binarySplit (0, _N-1, k-1, add);
//  _boundary[0].emplace_back (0);
//  std::sort(_boundary[0].begin (), _boundary[0].end (), xxx);
  _boundary.emplace_back (0, 0);
  sort(_boundary.begin (), _boundary.end (), cmpBoundary);
  for (int i = 0; i < _boundary.size (); i++) {
    if (i == _boundary.size () - 1) {
//      _boundary[1].emplace_back (_N - 1);
      _boundary[i].second = _N - 1;
    }
    else {
//      _boundary[1].emplace_back (_boundary[0][i + 1] - 1);
      _boundary[i].second = _boundary[i+1].first - 1;
    }
  }
}


void DetectorBinary::binarySplit (int start, int end, int k, bool add)
{
  if (add)
    _binaryTree->add (start, end, k);
  
  if (k == 0)
    return;
  else {
    int midPos = _minIndexArray[start][end][k];
    int leftK = _leftKArray[start][end][k];
    _boundary.emplace_back (midPos+1, -1);
    binarySplit (start, midPos, leftK, add);
    binarySplit (midPos+1, end, k-leftK, add);
  }
}



void DetectorBinary::calculateD (TreeNode &node)
{
  int start = node._val[0];
  int end = node._val[1];
  node._D = (double) _edgeCount->coeff(start, end) / ((end - start + 1) * (end - start) * .5);
  if (node._left != NULL)
    calculateD (*node._left);
  if (node._right != NULL)
    calculateD (*node._left);
}


void DetectorBinary::calculateDensity (TreeNode &node)
{
  int start = node._val[0];
  int end = node._val[1];
  
  if (node._left == NULL && node._right==NULL) {
    node._info = minusParent (node._D, node);
  }
  else {
    int leftStart = node._left->_val[0];
    int leftEnd = node._left->_val[1];
    int rightStart = node._right->_val[0];
    int rightEnd = node._right->_val[1];
    int delta = (end - start + 1) * (end - start) * .5;
    int deltaLeft = (leftEnd - leftStart + 1) * (leftEnd - leftStart) * .5;
    int deltaRight = (rightEnd - rightStart + 1) * (rightEnd - rightStart) * .5;
    double densitySum = (delta * node._D - deltaLeft * node._left->_D - deltaRight * node._right->_D) / (double) (delta - deltaLeft - deltaRight);
    node._info = minusParent (densitySum, node);
  }
  if (node._left != NULL)
    calculateDensity (*node._left);
  if (node._right != NULL)
    calculateDensity (*node._right);
}


double DetectorBinary::minusParent (double d, TreeNode &node)
{
  TreeNode *currentNode = &node;
  while (!(*currentNode == _binaryTree->root ())) {
    currentNode = currentNode->_parent;
    d -= currentNode->_info;
  }
  return d;
}


void DetectorBinary::filterNodes ()
{
//  for (int i = 0; i < _nodeList->size (); i++) {
//    (*_nodeList)[i]->_size = (*_nodeList)[i]->_val[1] - (*_nodeList)[i]->_val[0] + 1;
//  }
  Eigen::MatrixXi scoreMat(_nodeList->size (), _nodeList->size ());
  scoreMat.setZero ();
  
  int time = 0;
  std::vector<std::pair<int, TreeNode *>> nodeList1, nodeList2;
  nodeList1.reserve (_nodeList->size () / 2);
  nodeList2.reserve (_nodeList->size () / 2);
  
  double *ab1 = new double[2];
  double *ab2 = new double[2];
  
  while (time < 1000) {
    nodeList1.clear ();
    nodeList2.clear ();
    double oldAB1[2] = {0, 0};
    double oldAB2[2] = {0, 0};
    bool converged = false;
    
    // init filter
    for (int i = 0; i < _nodeList->size (); i++) {
      if (utils::randInt (0, 2) == 0)
        nodeList1.emplace_back (i, (*_nodeList)[i]);
      else
        nodeList2.emplace_back (i, (*_nodeList)[i]);
    }
    
    while (!converged) {
      simpleLinearRegression (nodeList1, ab1);
      simpleLinearRegression (nodeList2, ab2);
      if (ab1 == oldAB1 && ab2 == oldAB2)
        converged = true;
      else {
        nodeList1.clear ();
        nodeList2.clear ();
        for (int i = 0; i < _nodeList->size (); i++) {
          TreeNode *nodeTmp = (*_nodeList)[i];
          double dist1 = pow (ab1[0] * getX(*nodeTmp) + ab1[1] - nodeTmp->_info, 2);
          double dist2 = pow (ab2[0] * getX(*nodeTmp) + ab2[1] - nodeTmp->_info, 2);
          if (dist1 < dist2)
            nodeList1.emplace_back (i, nodeTmp);
          else
            nodeList2.emplace_back (i, nodeTmp);
        }
      }
    }
    
    if (ab1[0] < ab2[0])
      _trueNodes = nodeList1;
    else
      _trueNodes = nodeList2;
    
    if (abs(ab1[0] - ab2[0]) > .5) {
      time++;
      for (int m = 0; m < _trueNodes.size()-1; m++) {
        for (int n = m+1; n < _trueNodes.size (); n++) {
          scoreMat(_trueNodes[m].first, _trueNodes[n].first) ++;
        }
      }
    }
  }
}


double DetectorBinary::getX (TreeNode &node)
{
  double size = node._val[1] - node._val[0] + 1;
  return 1. / 3. * (size + 1);
}


double DetectorBinary::getY (TreeNode &node)
{
  return node._info;
}


void DetectorBinary::simpleLinearRegression (std::vector<std::pair<int, TreeNode *>> &nodeList, double *ab)
{
  double sumX = 0;
  double sumY = 0;
  for (int i = 0; i < nodeList.size (); i++) {
    sumX += getX (*nodeList[i].second);
    sumY += getY (*nodeList[i].second);
  }
  double meanX = sumX / nodeList.size ();
  double meanY = sumY / nodeList.size ();
  double covXY = 0;
  double varX = 0;
  for (int i = 0; i < nodeList.size (); i++) {
    covXY += (getX(*nodeList[i].second) - meanX) * (getY(*nodeList[i]) - meanY);
    varX += pow(getX(*nodeList[i].second) - meanX, 2);
  }
  
  ab[0] = covXY / varX;
  ab[1] = meanY - ab[0] * meanX;
}