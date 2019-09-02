//
// Created by wang mengbo on 2019-09-02.
//

#include "binaryTree.h"


BinaryTree::BinaryTree ()
{
  _root = NULL;
  _nodeList.reserve (1000);
}


BinaryTree::~BinaryTree ()
{
  for (int i = 0; i < _nodeList.size (); i++) {
    delete _nodeList[i];
  }
}


void BinaryTree::add (int start, int end,  int k)
{
  TreeNode *treeNode = new TreeNode(start, end);
  if (k == 0) {
    TreeNode *treeExistNode = _t.top();
    if (treeExistNode->_left == NULL) {
      treeExistNode->_left = treeNode;
      treeNode->_parent = treeExistNode;
    }
    else {
      treeExistNode->_right = treeNode;
      treeNode->_parent = treeExistNode;
      _t.pop ();
    }
  }
  else {
    if (_root == NULL) {
      _root = treeNode;
      _t.push (_root);
    }
    else {
      TreeNode *treeExistNode = _t.top ();
      if (treeExistNode->_left == NULL) {
        treeExistNode->_left = treeNode;
        treeNode->_parent = treeExistNode;
        _t.push (treeExistNode->_left);
      }
      else {
        treeExistNode->_right = treeNode;
        treeNode->_parent = treeExistNode;
        _t.pop ();
        _t.push (treeExistNode->_right);
      }
    }
  }
  _nodeList.emplace_back (treeNode);
}
