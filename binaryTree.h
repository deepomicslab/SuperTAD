//
// Created by wang mengbo on 2019-09-02.
//

#ifndef PROGRAM_BINARYTREE_H
#define PROGRAM_BINARYTREE_H

#include <vector>
#include <stack>
#include <iostream>
#include <limits>
#include <cmath>
#include "data.h"
#include "multiTree.h"
#include "params.h"
#include "utils.h"


namespace SuperTAD::binary {

    struct TreeNode {
        TreeNode *_left=NULL, *_right=NULL, *_parent=NULL;
        double _se=0, _info=0, _D=0, _vol=0;
        int _val[2] = {0, 0}, _idx;

        TreeNode(int start, int end)
        {
            _val[0] = start;
            _val[1] = end;
        }

        TreeNode(int start, int end, Data &data) : TreeNode(start, end)
        {
            this->setVol(data);
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
            _vol = copy._vol;
            _idx = copy._idx;
            return *this;
        }

        void setIdx(int idx)
        {
            _idx = idx;
        }

        // calculate structure entropy as node consists of only leaves, given parent node and return
        double getSEasLeaf(Data &data, TreeNode &pNode)
        {
            double se = data.getSE(_val[0], _val[1], pNode._val[0], pNode._val[1]);

            se += (data.getVol(_val[0], _val[1]) * data._logVolTable[_val[0]][_val[1]-_val[0]]) / data._doubleEdgeSum;

            if (_val[0] == 0)
                se -= data._sumOfGtimesLogG[_val[1]] / data._doubleEdgeSum;
            else
                se -= (data._sumOfGtimesLogG[_val[1]] - data._sumOfGtimesLogG[_val[0]-1]) / data._doubleEdgeSum;

            return se;
        }

        // calculate structure entropy, given parent node and return
        double getSE(Data &data, TreeNode &pNode)
        {
            double se = data.getSE(_val[0], _val[1], pNode._val[0], pNode._val[1]);
            if (_left)
                se += _left->getSE(data, *this);
            if (_right)
                se += _right->getSE(data, *this);
            return se;
        }

        // set volume
        void setVol(Data &data)
        {
            _vol = data.getVol(_val[0], _val[1]);
        }

        bool operator==(const TreeNode &t) const
        {
            return _val[0] == t._val[0] && _val[1] == t._val[1];
        }
    };

    inline std::ostream& operator<<(std::ostream &os, const TreeNode &node)
    {
        os << "idx=" << node._idx << ", ";
        os << "self=(" << node._val[0] << ", " << node._val[1] << ")";
        if (node._left != NULL) {
            os << ", left=(" << node._left->_val[0] << ", " << node._left->_val[1] << ")";
        }
        else {
            os << ", no left child";
        }

        if (node._right != NULL) {
            os << ", right=(" << node._right->_val[0] << ", " << node._right->_val[1] << ")";
        }
        else {
            os << ", no right child";
        }

        os << ", vol=" << node._vol;
        return os;
    }


    class Tree {
    private:
        std::stack<TreeNode*> _t;

    public:
        Data *_data;
        TreeNode *_root;
        std::vector<TreeNode*> _nodeList;

        Tree();

        ~Tree();

        void setData(Data &d);

        // BUG; deprecated; use insert instead
        void add(int start, int end, int k);

        void insert(TreeNode *newNode, TreeNode *parentNode);

        TreeNode &root() { return *_root; }

//        double getSE(TreeNode &child, TreeNode &parent);
    };

    static const int PruneMethod1 = 1;
    static const int PruneMethod2 = 2;

    class BasePruner {
    public:
        Data *_data;
        double **_minHtable, _optimalSE;
        int **_minIdxTable, _K, _mu, _optimalK;
        binary::Tree *_tree;
        multi::Tree _prunedTree;

        BasePruner(Tree &tree);

        virtual ~BasePruner();

        virtual void execute() { fprintf(stderr, "execute() in BasePruner should not be called\n"); };
    };


    // legacy!!! prune deep binary tree
    class Pruner1 : public BasePruner {
    public:
        Pruner1(Tree &tree, int k=10);

        ~Pruner1() override;

        void execute() override;

        double getH(TreeNode &node, int k);

        void backTrace(TreeNode &node, int k);
    };


    // prune deep binary tree
    // H(mu, k)=min_k1{H(mu_l, k1)+H(mu_r, k-k2)}
    class Pruner2 : public BasePruner {
    public:
        Pruner2(Tree &tree);

        ~Pruner2() override;

//        void init();

        void execute() override;

        // recursive approach for optimal H given k
        double getH(TreeNode &node, int k);

        void backTrace(TreeNode &node, int k);
    };

}

#endif //PROGRAM_BINARYTREE_H
