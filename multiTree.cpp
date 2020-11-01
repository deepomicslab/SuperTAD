//
// Created by wang mengbo on 2019-09-03.
//

#include "multiTree.h"
#include "params.h"


namespace multi {

//  bool operator<(const TreeNode &t1, const TreeNode &t2)
//  {
//    return t1._val[1] < t2._val[0];
//  }


//  std::ostream& operator<< (std::ostream &os, const TreeNode &node)
//  {
//    os << "self=(" << node._val[0] << ", " << node._val[1] << ")";
//    os << ", info=" << node._info << ", len(children)=" << node._children.size ();
//    if (node._parent) {
//        os << ", parent=(" << node._parent->_val[0] << ", " << node._parent->_val[1] << ")";
//    }
//    else {
//        os << ", no parent";
//    }
//    return os;
//  }


    Tree::Tree()
    {
        _root = new TreeNode(0, SuperTAD::_N_);
//    _nodeList.emplace_back (_root);
    }


    Tree::~Tree()
    {
//        for (int i = 0; i < _nodeList.size (); i++)
//            delete _nodeList[i];
        delete _root;
    }


    void Tree::insert(int start, int end)
    {
        TreeNode *treeNode = new TreeNode(start, end);
        if (!add(*treeNode, *_root)) {
            fprintf(stderr, "cannot add node: %s\n", treeNode->verbose().c_str());
        }
    }


    bool Tree::add(multi::TreeNode &newNode, multi::TreeNode &parentNode)
    {
        if (parentNode._val[0] <= newNode._val[0]
            && parentNode._val[1] >= newNode._val[1]
            && !parentNode.hasBiggerChild(newNode)
            && parentNode._children.emplace(&newNode).second) {
            newNode._parent = &parentNode;
//            _nodeList.emplace_back(&newNode);

//            std::cout << "parentNode=\n";
//            treeNodeVerbose(parentNode);
//            std::cout << "newNode=" << newNode << "\n";
//            std::cout << "\n";
            return true;
        }
        else {
            for (auto it: parentNode._children) {
                if (add(newNode, *it))
                    return true;
            }
        }
        return false;
    }


    void Tree::getNodeList(std::vector<TreeNode *> &nl)
    {
//        nl.emplace_back(_root);
        _root->getChildren(nl);
    }

}