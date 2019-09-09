//
// Created by wang mengbo on 2019-09-03.
//

#ifndef PROGRAM_MULTITREE_H
#define PROGRAM_MULTITREE_H

#include <vector>
#include <set>
#include <iostream>
#include "params.h"


namespace multi {
  
  struct TreeNode {
    int _val[2];
    double _info=0;
    std::set<TreeNode *> _children;
    TreeNode *_parent=NULL;
    
    TreeNode (int start, int end) {
      _val[0] = start;
      _val[1] = end;
    }
    
    TreeNode& operator=(const TreeNode &copy) {
      _val[0] = copy._val[0];
      _val[1] = copy._val[1];
      _info = copy._info;
      _children = copy._children;
      _parent = copy._parent;
      return *this;
    }
  
    bool operator==(const TreeNode &t) const {
      return _val[0] == t._val[0] && _val[1] == t._val[1];
    }
  };
  
  
  bool operator<(const TreeNode &t1, const TreeNode &t2);
  
  std::ostream& operator<< (std::ostream &os, const TreeNode &node);
  
  
  class Tree {
  private:
    TreeNode *_root;
    std::vector<TreeNode *> _nodeList;
    
  public:
    Tree ();
    
    ~Tree ();
    
    bool add (TreeNode &parentNode, TreeNode &newNode);
  
    std::vector<TreeNode *> &nodeList () { return _nodeList; }
    
    void insert (int start, int end);
  };
}


#endif //PROGRAM_MULTITREE_H
