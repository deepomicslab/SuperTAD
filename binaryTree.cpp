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
        delete _root;
    }


    void Tree::setData(Data &d)
    {
        _data = &d;
    }


    void Tree::add(int start, int end, int k)
    {
        if (_VERBOSE_)
            fprintf(stderr, "binary::Tree::add() is deprecated and going to be removed. Please use insert instead.\n");

        TreeNode *treeNode = new TreeNode(start, end, *_data);

        // add leaf nodes (only 1 bin)
        if (k == 0) {
            if (SuperTAD::_DEBUG_)
                printf("leaf node: (%d, %d)\n", start, end);

            TreeNode *treeExistNode = _t.top();
            if (treeExistNode->_left == NULL) {
                treeExistNode->_left = treeNode;
                treeNode->_parent = treeExistNode;
            }
            else {
                treeExistNode->_right = treeNode;
                treeNode->_parent = treeExistNode;
                _t.pop();
            }
        } else {    // add regular nodes
            if (_root == NULL) {
                _root = treeNode;
                _t.push(_root);
                _root->setIdx(_nodeList.size());
//                _nodeList.emplace_back(_root);
            }
            else {
                TreeNode *treeExistNode = _t.top();
                if (treeExistNode->_left == NULL) {
                    treeExistNode->_left = treeNode;
                    treeNode->_parent = treeExistNode;
                    _t.push(treeExistNode->_left);
                }
                else {
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


    void Tree::insert(TreeNode *newNode, TreeNode *parentNode)
    {
        if (parentNode->_left == NULL) {
            parentNode->_left = newNode;
            newNode->_parent = parentNode;
            newNode->setIdx(_nodeList.size());
            newNode->setVol(*_data);
            _nodeList.emplace_back(newNode);
        }
        else if (newNode->_val[1] <= parentNode->_left->_val[1]) {
            insert(newNode, parentNode->_left);
        }
        else if (parentNode->_right == NULL) {
            parentNode->_right = newNode;
            newNode->_parent = parentNode;
            newNode->setIdx(_nodeList.size());
            _nodeList.emplace_back(newNode);
        }
        else {
            insert(newNode, parentNode->_right);
        }
    }

    BasePruner::BasePruner(Tree &tree)
    {
        _data = tree._data;
        _tree = &tree;
        _mu = _tree->_nodeList.size();
        _prunedTree.setData(*_data);
    }

    BasePruner::~BasePruner()
    {
        _data = NULL;
        _tree = NULL;
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


//    void Pruner2::init()
//    {
//        for (int i = 0; i < _N_; i++) {
//            std::fill(&_minHtable[i][0], &_minHtable[i][_mu], std::numeric_limits<double>::infinity());
//            std::fill(&_minIdxTable[i][0], &_minIdxTable[i][_mu], -1);
//        }
//    }


    void Pruner2::execute()
    {
        double minSE = std::numeric_limits<double>::infinity();
        int minIdx;

//        // test
//        int k=63;
//        getH(*_tree->_root, k);
//        double tmpMinSE = _minHtable[k - 1][_tree->_root->_idx];
//        printf("========\nk=%d, minSE=%f\n", k, tmpMinSE);
//        backTrace(*_tree->_root, k);
//        _prunedTree._root->setVol(*_data);
//        std::cout << "prunedTree.root:" << *_prunedTree._root << "\n";
//        double se=0;
//        printf("root vol=%f\n", _prunedTree._root->_vol);
//        std::ofstream file;
//        file.open("../data/test/edge_graph_3421.txt.deepbinary.pruned.tsv.test.log.tsv");
//        file << "start\tend\tvol\tse\n";
//        for (auto it:_prunedTree._nodeList) {
//            file << it->_val[0] << "\t" << it->_val[1] << "\t" << it->_vol << "\t" << it->_se << "\n";
//            std::cout << *it << "\n";
//            double tmp=_data->getSE(it->_val[0],it->_val[1], _prunedTree._root->_vol, it->_vol);
//            printf("se=%f\n", tmp);
//            se += tmp;
//            printf("\n");
//        }
//        printf("\n");
//        for (auto it:_prunedTree._nodeList) {
////            file << it->_val[0] << "\t" << it->_val[1] << "\t" << it->_vol << "\t" << it->_se << "\n";
////            std::cout << *it << "\n";
////            double tmp=_data->getSE(it->_val[0],it->_val[1], _prunedTree._root->_vol, it->_vol);
////            printf("se=%f\n", tmp);
////            se += tmp;
//            double tmp;
//            for (int i=it->_val[0]; i<=it->_val[1]; i++) {
//                tmp = _data->getSE(i, i, it->_vol);
//                printf("se of bin %d=%f\n", i, tmp);
//                file << i << "\t" << i << "\t" << it->_vol << "\t" << tmp << "\n";
//                se += tmp;
//            }
//        }
//        file.close();
//        printf("final se=%f\n", se);
//        _optimalK = k;
//        _optimalSE = se;
////        exit(0);

        // aggressive prune mode
        if (_TURBO_PRUNE_) {
//            std::clock_t tTmp = std::clock();
            for (int k = 1; k <= _N_; k++) {
                getH(*_tree->_root, k);
                double tmpMinSE = _minHtable[k - 1][_tree->_root->_idx];
                printf("========\nk=%d, minSE=%f\n", k, tmpMinSE);
                if (minSE > tmpMinSE) {
                    printf("update min SE\n");
                    minSE = tmpMinSE;
                    minIdx = k;
                } else {
                    printf("met potential inflection\n");
                    break;
                }
            }
//            printf("getH for %d consumes %fs\n", _N_, (float)(std::clock() - tTmp)/CLOCKS_PER_SEC);
        } else {    // general prune mode
            std::clock_t tTmp = std::clock();
            for (int k=_N_; k>0; k--) {
                getH(*_tree->_root, k);
                double tmpMinSE = _minHtable[k - 1][_tree->_root->_idx];
                printf("========\nk=%d, minSE=%f\n", k, tmpMinSE);
                if (minSE >= tmpMinSE) {
                    printf("update min SE\n");
                    minSE = tmpMinSE;
                    minIdx = k;
                }
            }
//            printf("getH for %d consumes %fs\n", _N_, (float)(std::clock() - tTmp)/CLOCKS_PER_SEC);
        }
        _optimalK = minIdx;
        _optimalSE = minSE;
        printf("optimal k is %d with minial se %f\n", _optimalK, _optimalSE);
        backTrace(*_tree->_root, _optimalK);
    }


    double Pruner2::getH(TreeNode &node, int k)
    {
//        printf("getH(node=%d, k=%d)\n", node._idx, k);
        if (_minHtable[k-1][node._idx] != std::numeric_limits<double>::infinity()) {
//            printf("already obtained optimal H for node=%d and k=%d;optimal H=%f, k1=%d\n", node._idx, k, _minHtable[0][node._idx], _minIdxTable[k-1][node._idx]);
            return _minHtable[k-1][node._idx];
        } else {
            if (k == 1) {
                double tmp = node.getSEasLeaf(*(_tree->_data), *(_tree->_root));
                _minHtable[0][node._idx] = tmp;
//                std::cout << "for node: " << node;
//                printf(", getH(k=1, id=%d)=%f\n", node._idx, tmp);
                return tmp;
            } else {
                if (node._left == NULL || node._right == NULL) {
//                    printf("getH(node=%d, k=%d): no child exists\n", node._idx, k);
//                    std::cout << node << "\n";
                    return std::numeric_limits<double>::infinity();
                } else if (node._val[1]-node._val[0]+1<k) {
//                    printf("getH(node=%d, k=%d): k>#bins\n", node._idx, k);
//                    std::cout << node << "\n";
                    return std::numeric_limits<double>::infinity();
                }
                else {
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
//                        printf("getH(k=%d, id=%d)=%f, k1=%d\n", k - 1, node._idx, minH, minK1);
                    }
                    return minH;
                }
            }
        }
    }


    void Pruner2::backTrace(TreeNode &node, int k)
    {
        if (k == 1) {
//            std::cout << "selected node:" << node << "\n";
            multi::TreeNode *tn = _prunedTree.add(node._val[0], node._val[1]);
//            std::cout << "inserted node:" << *tn << "\n\n";
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
