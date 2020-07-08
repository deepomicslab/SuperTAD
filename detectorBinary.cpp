//
// Created by wang mengbo on 2019-09-01.
//

#include "detectorBinary.h"


namespace binary {

    Detector::Detector (Data &data)
    {
        _data = &data;
        _edgeCount = &data.edgeCount();
        _table = new double **[_N_];
        _minIndexArray = new int **[_N_];
        _leftKArray = new int **[_N_];
        _minIndexTableForBold = new int **[_N_];
        int size = 0;
        for (int s=0; s < _N_; s++) {
            _table[s] = new double *[_N_];
            _minIndexArray[s] = new int *[_N_];
            _leftKArray[s] = new int *[_N_];
            _minIndexTableForBold[s] = new int *[_N_];
            int k;
            for (int e=s; e < _N_; e++) {
                k = (e-s+1 < _K_ ? e-s+1 : _K_);
                size += k;
                _table[s][e] = new double[k]{};
                _minIndexArray[s][e] = new int[k]{};
                _leftKArray[s][e] = new int[k]{};
                if (_FAST_) {
                    _minIndexTableForBold[s][e] = new int[k]{};
                    memset(_minIndexTableForBold[s][e], -1, k * sizeof(int));
                }
            }
        }
        printf("db table size=%d\n", size);
        _boundary.reserve(_K_);
        _binaryTree = new binary::Tree();
        _numBins = new int(0);
        _kTmpIdx = new int(0);
        _kMinusKtmpIdx = new int(0);
    }


    Detector::~Detector ()
    {
        delete _binaryTree;
        for (int s = 0; s < _N_; s++) {
            for (int e = s; e < _N_; e++) {
                delete _table[s][e];
                delete _minIndexArray[s][e];
                delete _leftKArray[s][e];

                if (_FAST_) {
                    delete _minIndexTableForBold[s][e];
                }
            }
            delete _table[s];
            delete _minIndexArray[s];
            delete _leftKArray[s];
        }
        delete _table;
        delete _minIndexArray;
        delete _leftKArray;

        delete _numBins;
        delete _kTmpIdx;
        delete _kMinusKtmpIdx;
    }


    void Detector::execute ()
    {
        std::clock_t tTmp;

        fillTable();

        std::vector<double> sumOfEntropy;
        std::vector<double> sumOfLeaves;
        std::vector<intDoublePair> normLeaves;



        int index;
        if (_DETERMINE_K_) {
            if (_VERBOSE_) {
                printf("start determine k\n========\n");
                tTmp = std::clock();
            }

            double entropy;

            // determine K
            for (int k = 2; k <= _K_; k++) {
                if (_VERBOSE_)
                    std::cout << "k=" << k << std::endl;
                indexKtmp(k);

                entropy = _table[0][_N_ - 1][*_kTmpIdx];
                if (_VERBOSE_)
                    printf("min structure entropy=%f\n", entropy);
                sumOfEntropy.emplace_back(entropy);

                backTrace(k);
                double leafSum = 0;
                for (int leaf = 0; leaf < _boundary.size(); leaf++) {
                    leafSum += _data->getSE(_boundary[leaf].first, _boundary[leaf].second, 2.*_data->_edgeSum);
                    leafSum += _table[_boundary[leaf].first][_boundary[leaf].second][0];
                }
                sumOfLeaves.emplace_back(leafSum);
                double divisor = log2(_N_ / (double) k) + (_N_ * (k - 1) / (double) (k * (_N_ - 1))) * log2(k);
                normLeaves.emplace_back(k, leafSum / divisor);
                if (_VERBOSE_)
                    std::cout << "========\n";
            }

            sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePairBySecond);

            if (_VERBOSE_) {
                printf("finish determine k\n");
                printf("determination consumes %fs\n", (float)(std::clock()-tTmp)/CLOCKS_PER_SEC);
            }
            else
                printf("determine k\n");
        }
        else {

            int kTmp = _K_;
            indexKtmp(kTmp);

            double entropy = _table[0][_N_ - 1][*_kTmpIdx];
            sumOfEntropy.emplace_back(entropy);

            backTrace(kTmp);
            double leafSum = 0;
            int currentStart, currentEnd;
            for (int leaf = 0; leaf < _boundary.size(); leaf++) {
                currentStart = _boundary[leaf].first;
                currentEnd = _boundary[leaf].second;
                leafSum += _data->getSE(currentStart, currentEnd, 2. * _data->_edgeSum);
                leafSum += _table[currentStart][currentEnd][0];
            }
            sumOfLeaves.emplace_back(leafSum);
            double divisor = log2(_N_ / (double) kTmp) + (_N_ * (kTmp - 1) / (double) (kTmp * (_N_ - 1))) * log2(kTmp);
            normLeaves.emplace_back(kTmp, leafSum/divisor);
        }

        index = normLeaves[0].first;
        printf("obtain optimal structure at k=%d\n", index);
        fflush(stdout);
        backTrace(index, true);


        _nodeList = &_binaryTree->nodeList();
        Writer::writeTree(_OUTPUT_ + ".original_boundaries", *_nodeList);

        // filter
        if (_FILTERING_) {
            if (_VERBOSE_) {
                printf("start filtering\n");
                tTmp = std::clock();
            }
            else
                printf("filter nodes\n");

            calculateD (_binaryTree->root());
            calculateDensity(_binaryTree->root());
            filterNodes();
            std::vector<binary::TreeNode *> trueNodes;
            for (auto it = _trueNodeList.begin(); it != _trueNodeList.end(); it++) {
                trueNodes.emplace_back((*it));
            }

            Writer::writeTree(_OUTPUT_ + ".filter_boundaries", trueNodes);
            if (_VERBOSE_)
                printf("filtering consumes %fs\n", (float)(std::clock() - tTmp) / CLOCKS_PER_SEC);
        }

    }


    void Detector::fillTable()
    {
        std::clock_t t, tTmp;
        if (_VERBOSE_) {
            printf("start filling dp table\n");
            t = std::clock();
        }
        else
            printf("fill dp table\n");

        // process k=1
        if (_VERBOSE_) {
            printf("start filling base case\n");
            tTmp = std::clock();
        }

        for (int s = 0; s < _N_; s++) {
            for (int e = s; e < _N_; e++) {
                double currentVol = _data->getVol(s, e);
                double binSum;
                if (s == 0) {
                    binSum = _data->getGtimesLogG(currentVol) - _data->_sumOfGtimesLogG[e];
                } else {
                    binSum = _data->getGtimesLogG(currentVol) - (_data->_sumOfGtimesLogG[e] - _data->_sumOfGtimesLogG[s - 1]);
                }
                _table[s][e][0] = binSum / (2. * _data->_edgeSum);
            }
        }

        if (_VERBOSE_) {
            printf("finish filling base case where k=1, _table[0][%d][0]=%fs\nfilling base case consumes %fs\n",
                   _N_ - 1, _table[0][_N_ - 1][0], (float) (std::clock() - tTmp)/CLOCKS_PER_SEC);
        }

        // process k>=2
        if (_VERBOSE_) {
            printf("start filling upper cases\n");
            tTmp = std::clock();
        }

        int kIdx, k, s, e;
        bool breakFlag = false;

        for (k = 2; k <= _K_; k++) {

            if (breakFlag) {
                if (_VERBOSE_)
                    printf("break at loop k; s=%d, e=%d, k=%d\n", s, e, k);
                break;
            }

            indexK(k, kIdx);

            for (s = 0; s < _N_; s++) {
                if (breakFlag) {
                    if (_VERBOSE_)
                        printf("break at loop s; s=%d, e=%d, k=%d\n", s, e, k);
                    break;
                }

                for (e = s; e < _N_; e++) {
                    if (breakFlag) {
                        if (_VERBOSE_)
                            printf("break at loop e; s=%d, e=%d, k=%d\n", s, e, k);
                        break;
                    }

                    // skip cases when #numBins<#numLeves
                    numBins(s, e);
                    if (*_numBins < k) {
                        continue;
                    }

                    double minSE, se, parentVol;
                    int leftI, leftK;

                    minSE = std::numeric_limits<double>::infinity();
                    leftI = 0;
                    leftK = 0;

                    /*
                     * find min{S(s,i,kTmp)+H_l(s,e,i)+S(i+1,e,k-kTmp)+H_r(s,e,i)}
                     * loop all meaningful comb of i, kTmp
                     */
                    for (int kTmp=1; kTmp<k; kTmp++) {
                        indexKtmp(kTmp);
                        indexK(k-kTmp, *_kMinusKtmpIdx);

                        double minSE2 = std::numeric_limits<double>::infinity();
                        int leftI2 = 0;

                        int endTmp = (_minIndexTableForBold[s][e][*_kTmpIdx] == -1 ?
                                      e : _minIndexTableForBold[s][e][*_kTmpIdx] + _PENALTY_);
                        if (endTmp > e)
                            endTmp = e;

                        for (int i=s; i<endTmp; i++) {
                            numBins(s, i);
                            if (*_numBins < kTmp) {
//                                printf("s=%d, i=%d, e=%d, #numBins=%d, k=%d, kTmp=%d; skip\n",
//                                    s, i, e, *_numBins, k, kTmp);
                                continue;
                            }
                            numBins(i+1, e);
                            if (*_numBins < k-kTmp) {
//                                printf("s=%d, i+1=%d, e=%d, #numBins=%d, k=%d, kTmp=%d, k-kTmp=%d; skip\n",
//                                       s, i + 1, e, *_numBins, k, kTmp, k - kTmp);
                                continue;
                            }

                            se = _table[s][i][*_kTmpIdx];
                            se += _table[i+1][e][*_kMinusKtmpIdx];
                            parentVol = _data->getVol(s, e);
                            se += _data->getSE(s, i, parentVol);
                            se += _data->getSE(i+1, e, parentVol);
                            if (se < minSE) {
                                minSE = se;
                                leftI = i;
                                leftK = kTmp;
                            }
                            if (se < minSE2) {
                                minSE2 = se;
                                leftI2 = i;
                            }
                        }

                        _minIndexTableForBold[s][e][*_kTmpIdx] = leftI2;

                    }
                    _minIndexArray[s][e][kIdx] = leftI;
                    _table[s][e][kIdx] = minSE;
                    _leftKArray[s][e][kIdx] = leftK;

                    if (k == _K_ && s == 0 && e == _N_ - 1) {
                        if (_VERBOSE_)
                            printf("already obtain optimal solution for minS(1, %d, %d); stop here\n", _N_, _K_);
                        breakFlag=true;
                    }
                }
            }

            if (_VERBOSE_) {
                printf("finishing filling upper case where k=%d, _table[0][%d][%d]=%fs\n",
                       k, _N_ - 1, kIdx, _table[0][_N_ - 1][kIdx]);
            }
        }

        if (_VERBOSE_)
            printf("finish filling table\nfilling table consumes %fs\n", (float)(std::clock()-tTmp)/CLOCKS_PER_SEC);

        return;
    }


    void Detector::backTrace(int k, bool add)
    {
        std::clock_t tTmp;
        if (_VERBOSE_) {
            printf("start backtrace\n");
            tTmp = std::clock();
        }
        else
            printf("backtrace k=%d\n", k);

        init();

        binarySplit(0, _N_ - 1, k, add);

        _boundary.emplace_back(0, 0);

        sort(_boundary.begin(), _boundary.end(), utils::cmpBoundary);
        for (int i = 0; i < _boundary.size(); i++) {
            if (i == _boundary.size() - 1)
                _boundary[i].second = _N_ - 1;
            else
                _boundary[i].second = _boundary[i + 1].first - 1;
        }
        if (_VERBOSE_) {
            printf("finish backtrace\nbacktrace consumes %fs\n", (float)(std::clock()-tTmp)/CLOCKS_PER_SEC);
        }
    }


    void Detector::init()
    {
        _boundary.clear();
    }


    void Detector::binarySplit(int s, int e, int k, bool add, int lv)
    {
        indexKtmp(k);

        if (add)
            _binaryTree->add(s, e, *_kTmpIdx);

        if (k == 1)
            return;
        else {
            int i = _minIndexArray[s][e][*_kTmpIdx];

            int kTmp = _leftKArray[s][e][*_kTmpIdx];

            _boundary.emplace_back(i+1, -1);

            binarySplit(s, i, kTmp, add, lv+1);

            binarySplit(i+1, e, k-kTmp, add, lv+1);
        }
    }


    void Detector::calculateD (binary::TreeNode &node)
    {
        int start = node._val[0];
        int end = node._val[1];
        node._D = (double) _edgeCount->coeff(start, end) / ((end - start + 1) * (end - start) * .5);
        if (node._left != NULL)
            calculateD (*node._left);
        if (node._right != NULL)
            calculateD (*node._right);
    }


    void Detector::calculateDensity (binary::TreeNode &node)
    {
        int start = node._val[0];
        int end = node._val[1];

        if (node._left == NULL && node._right == NULL) {
            node._info = minusParent (node._D, node);
        }
        else {
            int leftStart = node._left->_val[0];
            int leftEnd = node._left->_val[1];
            int rightStart = node._right->_val[0];
            int rightEnd = node._right->_val[1];
            int delta = (end - start + 1) * (end - start) * .5;
            int deltaLeft = (leftEnd - leftStart + 1) * (leftEnd - leftStart) * .5;
            int deltaRight = (rightEnd - rightStart + 1) * (rightEnd - rightStart) * .5;
            double densitySum = (delta * node._D - deltaLeft * node._left->_D - deltaRight * node._right->_D) /
                                (double) (delta - deltaLeft - deltaRight);
            node._info = minusParent (densitySum, node);
        }

        if (node._left != NULL)
            calculateDensity(*node._left);
        if (node._right != NULL)
            calculateDensity(*node._right);
    }


    double Detector::minusParent(double d, binary::TreeNode &node)
    {
        binary::TreeNode *currentNode = &node;
        while (!(*currentNode == _binaryTree->root())) {
            currentNode = currentNode->_parent;
            d -= currentNode->_info;
        }
        return d;
    }


    void Detector::filterNodes()
    {
        Eigen::MatrixXi scoreMat(_nodeList->size(), _nodeList->size());
        scoreMat.setZero();

        std::vector<std::pair<int, binary::TreeNode *>> nodeList1, nodeList2, trueNodeList;
        nodeList1.reserve(_nodeList->size() / 2);
        nodeList2.reserve(_nodeList->size() / 2);
        trueNodeList.reserve(_nodeList->size() / 2);

        double ab1[2]{};
        double ab2[2]{};

        int totalItr = 1000;
        int threshold = 900;
        int time = 0;
        while (time < totalItr) {
            _trueNodeList.clear();
            double oldAB1[2]{};
            double oldAB2[2]{};
            bool converged = false;

            // init filter
            while (true) {
                nodeList1.clear();
                nodeList2.clear();
                for (int i = 0; i < _nodeList->size(); i++) {
                    if (utils::randInt(0, 2) == 0)
                        nodeList1.emplace_back(i, (*_nodeList)[i]);
                    else
                        nodeList2.emplace_back(i, (*_nodeList)[i]);
                }

                if (nodeList1.size() > 1 && nodeList2.size() > 1) {
                    break;
                }
            }

            int countTmp = 0;
            while (!converged and countTmp < totalItr) {
                countTmp++;

                utils::copyDoubleArray(ab1, oldAB1, 2);
                utils::copyDoubleArray(ab2, oldAB2, 2);

                if (!simpleLinearRegression(nodeList1, ab1))
                    break;

                if (!simpleLinearRegression(nodeList2, ab2))
                    break;

                if (utils::equalDoubleArrays(ab1, oldAB1, 2) && utils::equalDoubleArrays(ab2, oldAB2, 2))
                    converged = true;
                else {
                    nodeList1.clear();
                    nodeList2.clear();
                    for (int i = 0; i < _nodeList->size(); i++) {
                        binary::TreeNode *nodeTmp = (*_nodeList)[i];
                        double dist1 = pow(ab1[0] * getX(*nodeTmp) + ab1[1] - nodeTmp->_info, 2);
                        double dist2 = pow(ab2[0] * getX(*nodeTmp) + ab2[1] - nodeTmp->_info, 2);
                        if (dist1 < dist2)
                            nodeList1.emplace_back(i, nodeTmp);
                        else
                            nodeList2.emplace_back(i, nodeTmp);
                    }
                }
            }

            if (ab1[0] < ab2[0])
                trueNodeList = nodeList1;
            else
                trueNodeList = nodeList2;

            if (abs (ab1[0] - ab2[0]) > .5) {
                time++;
                for (int m = 0; m < trueNodeList.size() - 1; m++) {
                    for (int n = m + 1; n < trueNodeList.size(); n++) {
                        scoreMat(trueNodeList[m].first, trueNodeList[n].first)++;
                    }
                }
            }
        }

        for (int i = 0; i < _nodeList->size(); i++) {
            for (int j = i + 1; j < _nodeList->size(); j++) {
                if (scoreMat.coeff(i, j) > threshold) {
                    _trueNodeList.emplace((*_nodeList)[i]);
                    _trueNodeList.emplace((*_nodeList)[j]);
                }
            }
        }
        return;
    }


    double Detector::getX(binary::TreeNode &node)
    {
        double size = node._val[1] - node._val[0] + 1;
        return 1. / 3. * (size + 1);
    }


    double Detector::getY(binary::TreeNode &node)
    {
        return node._info;
    }


    bool Detector::simpleLinearRegression(std::vector<std::pair<int, binary::TreeNode *>> &nodeList, double ab[])
    {
        double sumX = 0;
        double sumY = 0;
        for (int i = 0; i < nodeList.size(); i++) {
            sumX += getX(*nodeList[i].second);
            sumY += getY(*nodeList[i].second);
        }
        double meanX = sumX / nodeList.size();
        double meanY = sumY / nodeList.size();
        double covXY = 0;
        double varX = 0;
        for (int i = 0; i < nodeList.size(); i++) {
            covXY += (getX(*nodeList[i].second) - meanX) * (getY(*nodeList[i].second) - meanY);
            varX += pow (getX(*nodeList[i].second) - meanX, 2);
        }
        ab[0] = covXY / varX;
        ab[1] = meanY - ab[0] * meanX;
        if (std::isnan (ab[0]) || std::isnan (ab[1]))
            return false;

        return true;
    }

}
