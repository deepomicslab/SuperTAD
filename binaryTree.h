//
// Created by wang mengbo on 2019-09-02.
//

#ifndef PROGRAM_BINARYTREE_H
#define PROGRAM_BINARYTREE_H


#include <vector>
#include <stack>
#include <iostream>
#include "params.h"

namespace binary {
  
  struct TreeNode {
    int _val[1];
    TreeNode *_left;
    TreeNode *_right;
    double _se;
    double _info;
    TreeNode *_parent;
    double _D;
    //  double _size;
    
    TreeNode (int start, int end)
    {
      _val[0] = start;
      _val[1] = end;
      _se = 0;
      _left = NULL;
      _right = NULL;
      _info = 0;
      _parent = NULL;
      _D = 0;
    }
    
    TreeNode &operator= (const TreeNode &copy)
    {
      _val[0] = copy._val[0];
      _val[1] = copy._val[1];
      _se = copy._se;
      _left = copy._left;
      _right = copy._right;
      _info = copy._info;
      _parent = copy._parent;
      _D = copy._D;
      return *this;
    }
    
    bool operator== (const TreeNode &t) const
    {
      return _val[0] == t._val[0] && _val[1] == t._val[1];
    }
  };
  
  std::ostream& operator<< (std::ostream &os, const TreeNode &node);
  
  class Tree {
  private:
    TreeNode *_root;
    std::stack<TreeNode *> _t;
    std::vector<TreeNode *> _nodeList;
  
  public:
    Tree ();
    
    ~Tree ();
    
    void add (int &start, int &end, int &k);
    
    std::vector<TreeNode *> &nodeList () { return _nodeList; }
    
    TreeNode &root () { return *_root; }
  };
}

#endif //PROGRAM_BINARYTREE_H
