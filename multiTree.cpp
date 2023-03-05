//
// Created by wang mengbo on 2019-09-03.
//

#include "multiTree.h"


namespace SuperTAD::multi {

    Tree::Tree()
    {
        _root = new TreeNode(0, SuperTAD::_N_-1);
//        _nodeList.emplace_back (_root);
    }


    Tree::Tree(Data &d) : Tree()
    {
        _data = &d;
        _root->setVol(d);
    }


    Tree::~Tree()
    {
        for (auto & i : _nodeList)
            delete i;
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
        if (! insert(*node, *_root)) {
            fprintf(stderr, "cannot add node: %s\n", node->verbose().c_str());
        }
        return node;
    }


    bool Tree::insert(multi::TreeNode &newNode, multi::TreeNode &parent)
    {
        if (!parent._children.empty()) {
            for (auto & i : parent._children)
            {
                if (i->_val[0] <= newNode._val[0] && i->_val[1] >= newNode._val[1])
                    if(insert(newNode, *i))
                        return true;
            }
        }
        if (parent._val[0] <= newNode._val[0] && parent._val[1] >= newNode._val[0]) {
            newNode._parent = &parent;
            parent._children.emplace_back(&newNode);
            newNode.setVol(*_data);
            newNode.setSE(*_data);
//            printf("start=%d, end=%d, vol=%f, se=%f\n", newNode._val[0], newNode._val[1], newNode._vol, newNode._se);
            _nodeList.emplace_back(&newNode);
//            std::cout << "emplaced node: " << newNode << "\n";
            return true;
        }
        return false;
    }

    void Tree::preOrderTraversal(multi::TreeNode &node)
    {
        _nodeList.emplace_back(&node);
//        printf("start %d, end %d, %d children\n", node._val[0], node._val[1], node._children.size());
        for (int i=node._children.size()-1; i>=0; i--)
        {
            preOrderTraversal(*node._children[i]);
        }
    }

    void Tree::getPreOrder()
    {
        _nodeList.clear();
        preOrderTraversal(*_root);
    }
}