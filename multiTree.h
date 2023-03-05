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
        double _info=0, _g=0, _vol=0, _se=0;
        int _index = 0;
        std::vector<TreeNode*> _children;
        TreeNode *_parent=NULL;

        TreeNode(int start, int end) {
            _val[0] = start;
            _val[1] = end;
        }

        ~TreeNode() {
            for (auto i : _children)
                delete i;
            _children.clear();
        }

        TreeNode& operator=(const TreeNode &copy) {
            _val[0] = copy._val[0];
            _val[1] = copy._val[1];
            _info = copy._info;
            _g = copy._g;
            _vol = copy._vol;
            _children = copy._children;
            _parent = copy._parent;
            return *this;
        }

        bool operator==(const TreeNode &t) const {
            return _val[0] == t._val[0] && _val[1] == t._val[1];
        }

        bool hasOverLapChild(const TreeNode &t) const {
            for (auto i : _children)
            {
                if ((i->_val[0] - t._val[1]) * (i->_val[1] - t._val[0]) < 0)
                    return true;
            }
            return false;
        }

        bool hasBiggerChild(const TreeNode &t) const {
            for (auto i : _children)
            {
                if (i->_val[0] <= t._val[0] && i->_val[1] >= t._val[0])
                    return true;
            }
            return false;
        }

        std::string verbose(int numHeadingSpace=0) const {
            std::stringstream iss;
            iss << "self=(" << _val[0] << ", " << _val[1] << ")";
            iss << ", info=" << _info << ", len(children)=" << _children.size();
            if (_parent) {
                iss << ", parent=(" << _parent->_val[0] << ", " << _parent->_val[1] << ")";
            }
            else {
                iss << ", no parent";
            }
            iss << ", vol=" << _vol;
            iss << ", se=" << _se;
            return iss.str();
        }

        void setG(Data &data)
        {
            _g = data.getG(_val[0], _val[1]);
        }

        void setVol(Data &data)
        {
            _vol = data.getVol(_val[0], _val[1]);
        }

        void setSE(Data &data)
        {
            _se = data.getSE(_val[0], _val[1], _parent->_vol, _vol);
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
        std::reverse(node._children.begin(), node._children.end());
    }


    class Tree {
    private:
//        std::vector<TreeNode*> _nodeList;

    public:
        Data *_data;
        TreeNode *_root = NULL;
        std::vector<TreeNode*> _nodeList;  // Not containing root

        Tree();

        Tree(Data &d);

        ~Tree();

        void setData(Data &d);

        TreeNode *add(int start, int end);

        bool insert(TreeNode &newNode, TreeNode &parent);

        void preOrderTraversal(TreeNode &node);

        void getPreOrder();

//        std::vector<TreeNode*> &nodeList() { return _nodeList; }

    };
}


#endif //PROGRAM_MULTITREE_H
