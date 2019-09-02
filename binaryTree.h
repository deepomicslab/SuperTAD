//
// Created by wang mengbo on 2019-09-02.
//

#ifndef PROGRAM_BINARYTREE_H
#define PROGRAM_BINARYTREE_H


#include <vector>
#include <stack>


struct TreeNode {
  int _val[2];
  TreeNode *_left;
  TreeNode *_right;
  double _info;
  TreeNode *_parent;
  double _D;
//  double _size;
  
  TreeNode (int start, int end) {
    _val[0] = start;
    _val[1] = end;
    _left = NULL;
    _right = NULL;
    _info = 0;
    _parent = NULL;
    _D = 0;
//    _size = 0;
  }
  
  TreeNode & operator=(const TreeNode &copy) {
    _val[0] = copy._val[0];
    _val[1] = copy._val[1];
    _left = copy._left;
    _right = copy._right;
    _info = copy._info;
    _parent = copy._parent;
    _D = copy._D;
//    _size = copy._size;
  }
  
  bool operator==(const TreeNode &t) const{
    return _val[0] == t._val[0] && _val[1] == t._val[1];
  }
};


class BinaryTree {
private:
  TreeNode *_root;
  std::stack<TreeNode *> _t;
  std::vector<TreeNode *> _nodeList;

public:
  BinaryTree ();
  
  ~BinaryTree ();
  
  void add (int start, int end, int k);
  
  std::vector<TreeNode *> &nodeList () { return _nodeList; }
  
  TreeNode &root () { return *_root; }
};

#endif //PROGRAM_BINARYTREE_H
