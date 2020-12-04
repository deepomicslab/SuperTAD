//
// Created by wang mengbo on 2019-09-02.
//

#include "binaryTree.h"


namespace SuperTAD::binary
{

//    std::ostream& operator<<(std::ostream &os, const TreeNode &node)
//    {
//        os << "idx=" << node._idx << ", ";
//        os << "self=(" << node._val[0] << ", " << node._val[1] << ")";
//        if (node._left != NULL && node._right != NULL) {
//            os << ", left=(" << node._left->_val[0] << ", " << node._left->_val[1] << ")";
//            os << ", right=(" << node._right->_val[0] << ", " << node._right->_val[1] << ")";
//        }
//        return os;
//    }


    Tree::Tree()
    {
        _root = NULL;

        // reserve 1000 pointers of nodes
        _nodeList.reserve(1000);
    }


//    Tree::Tree(Data &d)
//    {
//        Tree();
//        _data = &d;
//    }


    Tree::~Tree()
    {
        for (int i = 0; i < _nodeList.size(); i++) {
            delete _nodeList[i];
        }
    }


    void Tree::setData(Data &d)
    {
        _data = &d;
    }


    void Tree::add(int start, int end, int k)
    {
        TreeNode *treeNode = new TreeNode(start, end, *_data);

        // leaf node (only 1 bin)
        if (k == 0) {
            if (SuperTAD::_DEBUG_)
                printf("leaf node: (%d, %d)\n", start, end);

            TreeNode *treeExistNode = _t.top();
            if (treeExistNode->_left == NULL) {
                treeExistNode->_left = treeNode;
                treeNode->_parent = treeExistNode;
            } else {
                treeExistNode->_right = treeNode;
                treeNode->_parent = treeExistNode;
                _t.pop();
            }
        }
        // other nodes
        else {
            if (_root == NULL) {
                _root = treeNode;
                _t.push(_root);
                _root->setIdx(_nodeList.size());
                _nodeList.emplace_back(_root);
            } else {
                TreeNode *treeExistNode = _t.top();
                if (treeExistNode->_left == NULL) {
                    treeExistNode->_left = treeNode;
                    treeNode->_parent = treeExistNode;
                    _t.push(treeExistNode->_left);
                } else {
                    treeExistNode->_right = treeNode;
                    treeNode->_parent = treeExistNode;
                    _t.pop();
                    _t.push(treeExistNode->_right);
                }
            }
        }

        if (treeNode != _root) {
            treeNode->setIdx(_nodeList.size());
//            treeNode->setVol(*_data);
            _nodeList.emplace_back(treeNode);
        }
    }


    void Tree::insert(TreeNode *treeNode, TreeNode *parentNode)
    {
        if (parentNode->_left == NULL) {
            parentNode->_left = treeNode;
            treeNode->_parent = parentNode;
            treeNode->setIdx(_nodeList.size());
            treeNode->setVol(*_data);
            _nodeList.emplace_back(treeNode);
        } else if (treeNode->_val[1] <= parentNode->_left->_val[1]) {
            insert(treeNode, parentNode->_left);
        } else if (parentNode->_right == NULL) {
            parentNode->_right = treeNode;
            treeNode->_parent = parentNode;
            treeNode->setIdx(_nodeList.size());
            _nodeList.emplace_back(treeNode);
        } else {
            insert(treeNode, parentNode->_right);
        }
    }


//    TreeNode* Tree::getNode(int idx)
//    {
//        if (idx==0)
//            return _root;
//        int i=0;
//        TreeNode *node = _root;
//        while (node->_left || node->_right) {
//            i++;
//            if (i==idx)
//                return node->_left;
//            i++;
//            if (i==idx)
//                return node->_right;
//
//        }
//    }


    BasePruner::BasePruner(Tree &tree)
    {
        _data = tree._data;
        _tree = &tree;
        _mu = _tree->_nodeList.size();
    }


    Pruner1::Pruner1(Tree &tree, int k) : BasePruner(tree)
    {
        _K = k;
        _minHtable = new double *[_K];
        _minIdxTable = new int *[_K];
        for (int i = 0; i < _K; i++) {
            _minHtable[i] = new double[_mu]{};
            std::fill(&_minHtable[i][0], &_minHtable[i][_mu], std::numeric_limits<double>::infinity());
            _minIdxTable[i] = new int[_mu]{};
            std::fill(&_minIdxTable[i][0], &_minIdxTable[i][_mu], -1);
        }
    }


    Pruner1::~Pruner1()
    {
        for (int i = 0; i < _K; i++) {
            delete[] _minHtable[i];
            delete[] _minIdxTable[i];
        }
        delete[] _minHtable;
        delete[] _minIdxTable;
    }


    void Pruner1::execute()
    {
//        printf("before: minIdxtable:\n");
//        utils::print2Darray(_minIdxTable, _K, _Mu);
        getH(*_tree->_root, _K);
//        printf("minHtable:\n");
//        utils::print2Darray(_minHtable, _K, _Mu);
//        printf("minIdxtable:\n");
//        utils::print2Darray(_minIdxTable, _K, _Mu);
        backTrace(*_tree->_root, _K);
    }


    double Pruner1::getH(TreeNode &node, int k)
    {
        if (k == 1) {
            double tmp = node.getSE(*(_tree->_data), *(_tree->_root));
            // if leaf node is bin
            if (node._left == NULL || node._right == NULL) {
                // no extra op
            } else { // if leaf node is bin (NOT node contains bin(s))
                for (int i = node._val[0]; i <= node._val[1]; i++) {
                    tmp += _data->getSE(i, i, node._vol);
                }
//                std::cout << "k=" << k << ", current node chosen:" << no  de << ", se=" << tmp << "\n";
                _minHtable[k-1][node._idx] = tmp;
            }
            return tmp;
        } else {
            if (node._left == NULL || node._right == NULL) {
//                std::cout << "k=" << k << ", but current node is leaf node:" << node << "\n";
                return std::numeric_limits<double>::infinity();
            } else {
                double minH = std::numeric_limits<double>::infinity();
                int minK1 = -1;
                for (int k1 = 1; k1 < k; k1++) {
                    double tmp = getH(*node._left, k1) + getH(*node._right, k - k1);
//                    if (k==10) {
//                        std::cout << "split node: " << node << " into " << k1 << " and " << k-k1 << "\nse=" << tmp << "\n";
//                    }
                    if (tmp < minH) {
                        minH = tmp;
                        minK1 = k1;
                    }
                }

                if (minH < std::numeric_limits<double>::infinity()) {
//                    if (k==10) {
//                        std::cout << "k=" << k << ", current node:" << node << "\n";
//                        printf("minH=%f, minK1=%d\n", minH, minK1);
//                    }
                    _minHtable[k - 1][node._idx] = minH;
                    _minIdxTable[k - 1][node._idx] = minK1;
                    return minH;
                } else {
//                    std::cout << "k=" << k << ", current node cannot be splited into " << k << " nodes:" << node << "\n";
                    return std::numeric_limits<double>::infinity();
                }
            }
        }
    }


    void Pruner1::backTrace(TreeNode &node, int k)
    {
//        std::cout << "root: " << *(_tree->_root) << "\n";
//        _prunedTree.add(_tree->_root->_val[0], _tree->_root->_val[1]);

        if (k == 1) {
//            std::cout << "selected node:" << node << "\n";
            _prunedTree.add(node._val[0], node._val[1]);
            return;
        }

        int k1 = _minIdxTable[k - 1][node._idx];
//        printf("k=%d, k1=%d, k-k1=%d\n", k, k1, k-k1);
        if (node._left)
            backTrace(*node._left, k1);
        if (node._right)
            backTrace(*node._right, k - k1);
    }


    Pruner2::Pruner2(Tree &tree) : BasePruner(tree)
    {
        _minHtable = new double *[_N_];
        _minIdxTable = new int *[_N_];
        for (int i = 0; i < _N_; i++) {
            _minHtable[i] = new double[_mu]{};
            std::fill(&_minHtable[i][0], &_minHtable[i][_mu], std::numeric_limits<double>::infinity());
            _minIdxTable[i] = new int[_mu]{};
            std::fill(&_minIdxTable[i][0], &_minIdxTable[i][_mu], -1);
        }
//        init();
    }


    Pruner2::~Pruner2()
    {
        for (int i = 0; i < _N_; i++) {
            delete[] _minHtable[i];
            delete[] _minIdxTable[i];
        }
        delete[] _minHtable;
        delete[] _minIdxTable;
    }


    void Pruner2::init()
    {
        for (int i = 0; i < _N_; i++) {
            std::fill(&_minHtable[i][0], &_minHtable[i][_mu], std::numeric_limits<double>::infinity());
            std::fill(&_minIdxTable[i][0], &_minIdxTable[i][_mu], -1);
        }
    }


    void Pruner2::execute()
    {
        double minSE = std::numeric_limits<double>::infinity();
        int minIdx;
        for (int k=1; k<_N_; k++) {
            getH(*_tree->_root, k);
            double tmpMinSE = _minHtable[k - 1][_tree->_root->_idx];
            printf("========\nk=%d, minSE=%f\n", k, tmpMinSE);
            if (minSE > tmpMinSE) {
                minSE = tmpMinSE;
                minIdx = k;
            }
            else {
                printf("met potential inflection\n");
                break;
            }
        }
        _optimalK = minIdx;
        printf("optimal k is %d with minial se %f\n", minIdx, minSE);
        backTrace(*_tree->_root, minIdx);
    }


    double Pruner2::getH(TreeNode &node, int k)
    {
        if (k==1) {
            double tmp = node.getSEasLeaf(*(_tree->_data), *(_tree->_root));
            _minHtable[k-1][node._idx] = tmp;
            return tmp;
        } else {
            if (node._left==NULL || node._right==NULL) {
                return std::numeric_limits<double>::infinity();
            } else {
                double minH = std::numeric_limits<double>::infinity();
                int minK1;
                for (int k1 = 1; k1 < k; k1++) {
                    double tmp = getH(*node._left, k1) + getH(*node._right, k - k1);
                    if (tmp < minH) {
                        minH = tmp;
                        minK1 = k1;
                    }
                }
                if (minH < std::numeric_limits<double>::infinity()) {
                    _minHtable[k - 1][node._idx] = minH;
                    _minIdxTable[k - 1][node._idx] = minK1;
                    return minH;
                } else {
                    return std::numeric_limits<double>::infinity();
                }
            }
        }
    }


    void Pruner2::backTrace(TreeNode &node, int k)
    {
//        _prunedTree.add(_tree->_root->_val[0], _tree->_root->_val[1]);

        if (k == 1) {
            std::cout << "selected node:" << node << "\n";
            _prunedTree.add(node._val[0], node._val[1]);
            return;
        }

        int k1 = _minIdxTable[k - 1][node._idx];
//        printf("k=%d, k1=%d, k-k1=%d\n", k, k1, k-k1);
        if (node._left)
            backTrace(*node._left, k1);
        if (node._right)
            backTrace(*node._right, k-k1);
    }

}
