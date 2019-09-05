//
// Created by wang mengbo on 2019-09-03.
//

#include "multiTree.h"


namespace multi {
  bool operator<(const TreeNode &t1, const TreeNode &t2)
  {
    return t1._val[1] < t2._val[0];
  }
  
  
  Tree::Tree ()
  {
    _root = new TreeNode (0, _N);
    _nodeList.emplace_back (_root);
  }
  
  
  Tree::~Tree ()
  {
    for (int i = 0; i < _nodeList.size (); i++)
      delete _nodeList[i];
  }
  
  bool Tree::add (multi::TreeNode &parentNode, multi::TreeNode &newNode)
  {
    if (parentNode._val[0] <= newNode._val[0] && parentNode._val[1] >= newNode._val[1] && parentNode._children.emplace (&newNode).second) {
      _nodeList.emplace_back (&newNode);
      return true;
    }
    else {
      for (auto it: parentNode._children) {
        if (add (*it, newNode))
          return true;
      }
    }
    return false;
  }
  
  
  void Tree::insert (int start, int end)
  {
    TreeNode *treeNode = new TreeNode (start, end);
    add(*treeNode, *_root);
  }
}