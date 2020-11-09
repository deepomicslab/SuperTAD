//
// Created by wang mengbo on 2019-09-03.
//

#ifndef PROGRAM_MULTITREE_H
#define PROGRAM_MULTITREE_H

#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <sstream>
#include "params.h"
#include "data.h"


namespace SuperTAD::multi {

    struct TreeNode {
        int _val[2];
        double _info=0;
        std::set<TreeNode*> _children;
        TreeNode *_parent=NULL;

        TreeNode(int start, int end) {
            _val[0] = start;
            _val[1] = end;
        }

        ~TreeNode() {
            while (!_children.empty()) {
                auto it = _children.begin();
                delete *it;
                _children.erase(it);
            }
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

        bool hasOverLapChild(const TreeNode &t) const {
            for (auto it=_children.begin(); it!=_children.end(); it++) {
                if ((_val[0]-t._val[1])*(_val[1]-t._val[0]) < 0)
                    return true;
            }
            return false;
        }

        bool hasBiggerChild(const TreeNode &t) const {
            for (auto it=_children.begin(); it!=_children.end(); it++) {
                if ((*it)->_val[0] <= t._val[0] && (*it)->_val[1] >= t._val[1])
                    return true;
            }
            return false;
        }

        std::string verbose(int numHeadingSpace=0) const {
            std::stringstream iss;
            iss << "self=(" << _val[0] << ", " << _val[1] << ")";
            iss << ", info=" << _info << ", len(children)=" << _children.size ();
            if (_parent) {
                iss << ", parent=(" << _parent->_val[0] << ", " << _parent->_val[1] << ")";
            }
            else {
                iss << ", no parent";
            }
            return iss.str();
        }

        void getChildren(std::vector<TreeNode*> &nl) const {
            for (auto it=_children.begin(); it!=_children.end(); it++) {
                nl.emplace_back(*it);
                (*it)->getChildren(nl);
            }
        }

    };

    inline bool operator<(const TreeNode &t1, const TreeNode &t2)
    {
        return t1._val[1] < t2._val[0];
    }

    inline std::ostream& operator<<(std::ostream &os, const TreeNode &node)
    {
        os << node.verbose();
        return os;
    }

    inline void treeNodeVerbose(multi::TreeNode &node, int numHeadingSpace=0)
    {
        for (int i=0; i<numHeadingSpace; i++)
            std::cout << " ";
        std::cout << node << "\n";
        for (auto it=node._children.begin(); it!=node._children.end(); it++)
            treeNodeVerbose(**it, numHeadingSpace+2);
    }


    class Tree {
    private:
//        std::vector<TreeNode*> _nodeList;

    public:
        TreeNode *_root = NULL;
        std::vector<TreeNode*> _nodeList;

        Tree();

        ~Tree();

        void add(int start, int end);

        bool insert(TreeNode &newNode, TreeNode &parentNode);

//        std::vector<TreeNode*> &nodeList() { return _nodeList; }

        void getNodeList(std::vector<TreeNode *> &nl);

        double getSE();
    };
}


#endif //PROGRAM_MULTITREE_H
