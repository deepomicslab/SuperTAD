//
// Created by wang mengbo on 2019-09-03.
//

#include "multiTree.h"


namespace SuperTAD::multi {

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
        _root = new TreeNode(0, SuperTAD::_N_-1);
//    _nodeList.emplace_back (_root);
    }


    Tree::Tree(Data &d) : Tree()
    {
        _data = &d;
        _root->setVol(d);
    }


    Tree::~Tree()
    {
        for (int i = 0; i < _nodeList.size (); i++)
            delete _nodeList[i];
        delete _root;
    }


    void Tree::setData(Data &d)
    {
        _data = &d;
        _root->setVol(d);
    }


    TreeNode* Tree::add(int start, int end)
    {
        TreeNode *node = new TreeNode(start, end);
        if (!insert(*node, *_root)) {
            fprintf(stderr, "cannot add node: %s\n", node->verbose().c_str());
        }
        return node;
    }


//    bool Tree::insert(multi::TreeNode &newNode, multi::TreeNode &parentNode)
//    {
//        if (parentNode._val[0] <= newNode._val[0]
//            && parentNode._val[1] >= newNode._val[1]
//            && !parentNode.hasBiggerChild(newNode)
//            && parentNode._children.emplace(&newNode).second) {
//            newNode._parent = &parentNode;
//            _nodeList.emplace_back(&newNode);
//
////            std::cout << "parentNode=\n";
////            treeNodeVerbose(parentNode);
//            std::cout << "emplaced newNode: " << newNode << "\n";
//            return true;
//        }
//        else {
//            for (auto it: parentNode._children) {
//                if (insert(newNode, *it))
//                    return true;
//            }
//        }
//        return false;
//    }


    bool Tree::insert(multi::TreeNode &newNode, multi::TreeNode &parent)
    {
        if (parent._children.size() > 0) {
            for (auto it: parent._children) {
                if (insert(newNode, *it)) {
                    return true;
                }
            }
        }
        else if (parent._val[0] <= newNode._val[0] && parent._val[1] >= newNode._val[0]) {
            newNode._parent = &parent;
            newNode.setVol(*_data);
            newNode.setSE(*_data);
//            printf("vol=%f, se=%f\n", newNode._vol, newNode._se);
            _nodeList.emplace_back(&newNode);
//            std::cout << "emplaced node: " << newNode << "\n";
            return true;
        } else {
            return false;
        }
    }


    void Tree::getNodeList(std::vector<TreeNode *> &nl)
    {
//        nl.emplace_back(_root);
        _root->getChildren(nl);
    }


//    double Tree::getSE()
//    {
//        double se = 0;
//        for (auto it=_root->_children.begin(); it!=_root->_children.end(); it++) {
//            se +=
//        }
//    }
}