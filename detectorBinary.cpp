//
// Created by wang mengbo on 2019-09-01.
//

#include "detectorBinary.h"


namespace SuperTAD::binary
{

    Detector::Detector(SuperTAD::Data &data)
    {
        _data = &data;
        _table = new double **[SuperTAD::_N_];
        _minIndexTable = new int **[SuperTAD::_N_];
        _leftKtable = new int **[SuperTAD::_N_];
        _minIndexTableForBold = new int **[SuperTAD::_N_];
        int size = 0;
        int k, s, e;
        for (s = 0; s < SuperTAD::_N_; s++)
        {
            _table[s] = new double *[SuperTAD::_N_];
            _minIndexTable[s] = new int *[SuperTAD::_N_];
            _leftKtable[s] = new int *[SuperTAD::_N_];
            _minIndexTableForBold[s] = new int *[SuperTAD::_N_];
            for (e = s; e < SuperTAD::_N_; e++)
            {
                k = (e - s + 1 < SuperTAD::_K_ ? e - s + 1 : SuperTAD::_K_);
                size += k;
                _table[s][e] = new double[k]{};
                _minIndexTable[s][e] = new int[k]{};
                _leftKtable[s][e] = new int[k]{};
                if (SuperTAD::_FAST_)
                {
                    _minIndexTableForBold[s][e] = new int[k]{};
                    memset(_minIndexTableForBold[s][e], -1, k * sizeof(int));
                }
            }
        }
        _boundaries.reserve(SuperTAD::_K_);
        _binaryTree = new Tree();
        _binaryTree->setData(*_data);
    }


    Detector::~Detector()
    {
        delete _binaryTree;
        for (int s = 0; s < SuperTAD::_N_; s++)
        {
            for (int e = s; e < SuperTAD::_N_; e++)
            {
                delete _table[s][e];
                delete _minIndexTable[s][e];
                delete _leftKtable[s][e];

                if (SuperTAD::_FAST_)
                {
                    delete _minIndexTableForBold[s][e];
                }
            }
            delete _table[s];
            delete _minIndexTable[s];
            delete _leftKtable[s];
        }
        delete _table;
        delete _minIndexTable;
        delete _leftKtable;
    }


    void Detector::execute()
    {
        std::clock_t tTmp;

        fillTable();

        std::vector<double> sumOfEntropy;
        std::vector<double> sumOfLeaves;
        std::vector<IntDoublePair> normLeaves;

        double binSum, se;
        if (SuperTAD::_DETERMINE_K_)
        {
            if (SuperTAD::_OPTIMAL_K_ < SuperTAD::_K_)
                SuperTAD::_K_ = SuperTAD::_OPTIMAL_K_;

            printf("Determine the optimal K as %d\n", SuperTAD::_K_);

            backTrace(SuperTAD::_K_, true);
        } else
        {
            printf("K=%d\n", SuperTAD::_K_);

            setIndex(SuperTAD::_K_, _kTmpIdx);

            se = _table[0][SuperTAD::_N_ - 1][_kTmpIdx];
            sumOfEntropy.emplace_back(se);

            backTrace(SuperTAD::_K_, true);
            binSum = 0;

            for (int leaf = 0; leaf < _boundaries.size(); leaf++)
            {
                binSum += _data->getSE(_boundaries[leaf].first, _boundaries[leaf].second, 0, SuperTAD::_N_-1);
                binSum += _table[_boundaries[leaf].first][_boundaries[leaf].second][0];
            }
            sumOfLeaves.emplace_back(binSum);
            double divisor = log2(SuperTAD::_N_ / (double) SuperTAD::_K_) + (SuperTAD::_N_ * (SuperTAD::_K_ - 1) / (double) (
                            SuperTAD::_K_ * (SuperTAD::_N_ - 1))) * log2(SuperTAD::_K_);
            normLeaves.emplace_back(SuperTAD::_K_, binSum / divisor);
        }
        _nodeList = &_binaryTree->_nodeList;
        if (SuperTAD::_VERBOSE_)
        {
            printf("nodes:");
            for (int i = 0; i < _nodeList->size(); i++)
            {
                printf("(%d, %d)", (*_nodeList)[i]->_val[0], (*_nodeList)[i]->_val[1]);
                if (i < _nodeList->size() - 1)
                    printf(", ");
            }
            printf("\n");
        }

        SuperTAD::Writer::writeTree(SuperTAD::_OUTPUT_ + ".binary.original", *_nodeList);

        if (SuperTAD::_PRUNE_)
            filter();
    }


    bool Detector::sortByStart(Boundary a, Boundary b)
    {
        if (a.first == b.first)
            return a.size > b.size;
        else
            return a.first < b.first;
    }

    void Detector::executeFilter(std::string result)
    {
        SuperTAD::Reader::parseBoundariesIn8ColsFormat(_boundaries, result);
        sort(_boundaries.begin(), _boundaries.end(),
             Detector::sortByStart);    //sort boundaries in the increasing of start pos
        binary::TreeNode *newNode;
        int s, e;
        s = 0;
        e = SuperTAD::_N_ - 1;

        // add root
        newNode = new binary::TreeNode(s, e);
        _binaryTree->_root = newNode;

        // construct coding tree
        for (int i = 0; i < _boundaries.size(); i++)
        {
            newNode = new binary::TreeNode(_boundaries[i].first - 1, _boundaries[i].second - 1);
            _binaryTree->insert(newNode, &_binaryTree->root());
        }
        _nodeList = &_binaryTree->_nodeList;

        filter();
    }


    void Detector::fillTable()
    {
        std::clock_t t;
        if (SuperTAD::_VERBOSE_)
        {
            printf("start filling dp table\n");
            t = std::clock();
        }

        double binSum;
        for (int s = 0; s < SuperTAD::_N_; s++)
        {
            for (int e = s; e < SuperTAD::_N_; e++)
            {
                binSum = _data->getVol(s, e) * _data->_logVolTable[s][e-s];
                if (s == 0)
                {
                    binSum -= _data->_sumOfGtimesLogG[e];
                } else
                {
                    binSum -= _data->_sumOfGtimesLogG[e] - _data->_sumOfGtimesLogG[s - 1];
                }
                _table[s][e][0] = binSum / (_data->_doubleEdgeSum);
            }
        }

        int kIdx, k, s, e, leftI, leftK, leftI2, endTmp;
        double minSE, tmpSE, minSE2;
        _breakFlag = false; // for restricting the searching area based on identified boundaries, under testing
        double normOfLeavesTmp = std::numeric_limits<double>::infinity();

        for (k = 2; k <= SuperTAD::_K_; k++)
        {
            if (_breakFlag)
            {}

            setIndex(k, kIdx);

            for (s = 0; s < SuperTAD::_N_; s++)
            {
                if (_breakFlag)
                {}

                for (e = s; e < SuperTAD::_N_; e++)
                {
                    if (_breakFlag)
                    {}

                    setNumBins(s, e);
                    if (_numBins < k)
                        continue;

                    minSE = std::numeric_limits<double>::infinity();
                    leftI = 0;
                    leftK = 0;

                    /*
                     * find min{S(s,i,k')+H_l(s,e,i)+S(i+1,e,k-k')+H_r(s,e,i)}
                     * loop all meaningful comb of i, k'
                     */
                    for (int kTmp = 1; kTmp < k; kTmp++)
                    {
                        setIndex(kTmp, _kTmpIdx);
                        setIndex(k - kTmp, _kMinusKtmpIdx);

//                        if (SuperTAD::_FAST_)
//                        {
//                            minSE2 = std::numeric_limits<double>::infinity();
//                            leftI2 = 0;
////                            endTmp = (_minIndexTableForBold[s][e][*_kTmpIdx] == -1 ? e : _minIndexTableForBold[s][e][*_kTmpIdx] + SuperTAD::_PENALTY_);
//                            endTmp = (_minIndexTableForBold[s][e][_kTmpIdx] == -1 ? e :
//                                      _minIndexTableForBold[s][e][_kTmpIdx] + SuperTAD::_PENALTY_);
//                            if (endTmp > e) endTmp = e;
//                        } else
                        endTmp = e;

                        for (int i = s; i < endTmp; i++)
                        {
                            setNumBins(s, i);
                            if (_numBins < kTmp)
                                continue;
                            setNumBins(i + 1, e);
                            if (_numBins < k - kTmp)
                                continue;

                            tmpSE = _table[s][i][_kTmpIdx] + _table[i + 1][e][_kMinusKtmpIdx];

                            // H_l(s,e,i)
                            tmpSE += _data->getSE(s, i, s, e);

                            // H_r(s,e,i)
                            tmpSE += _data->getSE(i + 1, e, s, e);

//                            if (tmpSE - minSE < -SuperTAD::_THRESHOLD_)
                            if (tmpSE < minSE)
                            {
                                minSE = tmpSE;
                                leftI = i;
                                leftK = kTmp;
                            }

//                            if (SuperTAD::_FAST_ && tmpSE < minSE2)
//                            {
//                                minSE2 = tmpSE;
//                                leftI2 = i;
//                            }
                        }

//                        if (SuperTAD::_FAST_)
//                            _minIndexTableForBold[s][e][_kTmpIdx] = leftI2;
                    }
                    _minIndexTable[s][e][kIdx] = leftI;
                    _table[s][e][kIdx] = minSE;
                    _leftKtable[s][e][kIdx] = leftK;

//                    if (k == SuperTAD::_K_ && s == 0 && e == SuperTAD::_N_ - 1)
//                    {
//                        if (SuperTAD::_DEBUG_)
//                            printf("already obtain optimal solution for minS(1, %d, %d); stop here\n",
//                                   SuperTAD::_N_, SuperTAD::_K_);
//                        _breakFlag = true;
//                    }
                }
            }

            if (SuperTAD::_VERBOSE_)
            {
                printf("finish filling upper case where k=%d, _table[0][%d][%d]=%f\n",
                       k, SuperTAD::_N_ - 1, kIdx, _table[0][SuperTAD::_N_ - 1][kIdx]);
            }

            if (SuperTAD::_DETERMINE_K_)
            {
                backTrace(k);

                double leafSum = 0;

                for (int leaf = 0; leaf < _boundaries.size(); leaf++)
                {
                    leafSum += _data->getSE(_boundaries[leaf].first, _boundaries[leaf].second, 0, SuperTAD::_N_-1);
                    leafSum += _table[_boundaries[leaf].first][_boundaries[leaf].second][0];
                }
                double divisor = log2(SuperTAD::_N_ / (double) k) + (SuperTAD::_N_ * (k - 1) / (double) (k * (
                                SuperTAD::_N_ - 1))) * log2(k);
                double Tmp = leafSum / divisor;
                if (SuperTAD::_VERBOSE_)
                    printf("--------\noptimalK=%d, normLeaves=%f, se=%f\n", k, Tmp, _table[0][SuperTAD::_N_ - 1][kIdx]);
                if (Tmp - normOfLeavesTmp < - SuperTAD::_THRESHOLD_)
                {
                    SuperTAD::_OPTIMAL_K_ = k;
                    normOfLeavesTmp = Tmp;
                } else
                    break;
            }
        }

        if (SuperTAD::_VERBOSE_)
            printf("finish filling dp table, consumes %fs\n", (float) (std::clock() - t) / CLOCKS_PER_SEC);

    }


    void Detector::backTrace(int k, bool add)
    {
        init();

        if (add){
            binary::TreeNode *newNode = new binary::TreeNode(0, SuperTAD::_N_-1);
            _binaryTree->_root = newNode;
        }

        setIndex(k, _kTmpIdx);
        int i = _minIndexTable[0][SuperTAD::_N_-1][_kTmpIdx];
        int kTmp = _leftKtable[0][SuperTAD::_N_-1][_kTmpIdx];

        _boundaries.emplace_back(i + 1, -1);
        binarySplit(0, i, kTmp, add, 1);

        binarySplit(i + 1, SuperTAD::_N_-1, k - kTmp, add, 1);

        _boundaries.emplace_back(0, 0);

        sort(_boundaries.begin(), _boundaries.end(), utils::cmpBoundary);

        for (int i = 0; i < _boundaries.size(); i++)
        {
            if (i == _boundaries.size() - 1)
                _boundaries[i].second = SuperTAD::_N_ - 1;
            else
                _boundaries[i].second = _boundaries[i + 1].first - 1;
        }

    }


    void Detector::init()
    {
        _boundaries.clear();
    }


    void Detector::binarySplit(int s, int e, int k, bool add, int lv)
    {
        setIndex(k, _kTmpIdx);

        if (add)
        {
            binary::TreeNode *newNode = new binary::TreeNode(s, e);
            _binaryTree->insert(newNode, &_binaryTree->root());
        }

        if (k == 1)
            return;
        else
        {
            int i = _minIndexTable[s][e][_kTmpIdx];
            int kTmp = _leftKtable[s][e][_kTmpIdx];

            _boundaries.emplace_back(i + 1, -1);

            binarySplit(s, i, kTmp, add, lv + 1);

            binarySplit(i + 1, e, k - kTmp, add, lv + 1);
        }
    }


    void Detector::calculateD(binary::TreeNode &node)
    {
        int s = node._val[0];   //  start from index 0
        int e = node._val[1];
        //        node._D = (double) _edgeCount->coeff(start, end) / ((end - start + 1) * (end - start) * .5);

        if (s != e)
            node._D = (double) _data->_edgeCountTable[s][e] / ((e - s + 1) * (e - s) * .5);
        else
            node._D = 0;

        if (node._parent != NULL)   //  calculate se
        {
            node._se = _data->getSE(s, e,node._parent->_val[0], node._parent->_val[1]);
        }

        if (node._left != NULL)
            calculateD(*node._left);
        if (node._right != NULL)
            calculateD(*node._right);
    }


    void Detector::calculateDensity(binary::TreeNode &node)
    {
        int start = node._val[0];
        int end = node._val[1];

        if (node._left == NULL && node._right == NULL)
        {
            node._info = minusParent(node._D, node);
        } else
        {
            int leftStart = node._left->_val[0];
            int leftEnd = node._left->_val[1];
            int rightStart = node._right->_val[0];
            int rightEnd = node._right->_val[1];
            int delta = (end - start + 1) * (end - start) * .5;
            int deltaLeft = (leftEnd - leftStart + 1) * (leftEnd - leftStart) * .5;
            int deltaRight = (rightEnd - rightStart + 1) * (rightEnd - rightStart) * .5;
            double densitySum =
                    (delta * node._D - deltaLeft * node._left->_D - deltaRight * node._right->_D) /
                    (double) (delta - deltaLeft - deltaRight);
            node._info = minusParent(densitySum, node);
        }

        if (node._left != NULL)
            calculateDensity(*node._left);
        if (node._right != NULL)
            calculateDensity(*node._right);
    }


    double Detector::minusParent(double d, binary::TreeNode &node)
    {
        binary::TreeNode *currentNode = &node;
        while (!(*currentNode == _binaryTree->root()))
        {
            currentNode = currentNode->_parent;
            d -= currentNode->_info;
        }
        return d;
    }


    void Detector::filter()
    {
        std::clock_t tTmp;
        if (SuperTAD::_VERBOSE_)
        {
            printf("start filtering\n");
            tTmp = std::clock();
        }
        calculateD(_binaryTree->root());
        calculateDensity(_binaryTree->root());

        int *label1;    //  label from 1000 times random experiments
        label1 = filterNodes();
        for (int i = 0; i < _nodeList->size(); i++)
        {
            int size = (*_nodeList)[i]->_val[1] - (*_nodeList)[i]->_val[0] + 1;
            int parent_size = (*_nodeList)[i]->_parent->_val[1] - (*_nodeList)[i]->_parent->_val[0] + 1;
            double size_diff = std::abs(sqrt(size * (parent_size - size)) - size);
            parent_size = parent_size * 0.04;
            if (size_diff <= parent_size)
            {
                //                printf("size_diff <= threshold parentsize*0.04, %d, %d \n", (*_nodeList)[i]->_val[0], (*_nodeList)[i]->_val[1]);
                if (label1[i] > 0 and (*_nodeList)[i]->_se > (*_nodeList)[i]->_parent->_se)
                {
                    _trueNodeList.emplace_back((*_nodeList)[i]);
                    //                        printf("both labels are 1, %d, %d, \n", (*_nodeList)[i]->_val[0], (*_nodeList)[i]->_val[1]);
                }
            } else
            {
                _trueNodeList.emplace_back((*_nodeList)[i]);
                //                    printf("size_diff > threshold parentsize*0.04, %d, %d \n", (*_nodeList)[i]->_val[0], (*_nodeList)[i]->_val[1]);
            }
        }
        if (SuperTAD::_VERBOSE_)
            printf("filtering consumes %fs\n", (float) (std::clock() - tTmp) / CLOCKS_PER_SEC);
        printf("%d TAD candidates filtered into %d true TADs\n", _nodeList->size(), _trueNodeList.size());

        SuperTAD::Writer::writeTree(SuperTAD::_OUTPUT_ + ".binary.filter", _trueNodeList);
    }


    int *Detector::filterNodes()
    {
        _scoreTable = new float[_nodeList->size()];
        for (int i = 0; i < _nodeList->size(); i++)
        {
            _scoreTable[i] = 0;
        }

        std::vector<std::pair<int, binary::TreeNode *>> nodeList1, nodeList2, trueNodeList;
        nodeList1.reserve(_nodeList->size() / 2);
        nodeList2.reserve(_nodeList->size() / 2);
        trueNodeList.reserve(_nodeList->size() / 2);

        double ab1[2]{}, ab2[2]{};

        int totalItr = 1000;
        int time = 0;
        while (time < totalItr)
        {
            _trueNodeList.clear();
            double oldAB1[2]{}, oldAB2[2]{};
            bool converged = false;

            // init filter
            while (true)
            {
                nodeList1.clear();
                nodeList2.clear();
                for (int i = 0; i < _nodeList->size(); i++)
                {
                    if (utils::randInt(0, 2) == 0)
                        nodeList1.emplace_back(i, (*_nodeList)[i]); //  index & Treenode
                    else
                        nodeList2.emplace_back(i, (*_nodeList)[i]);
                }

                if (nodeList1.size() > 1 && nodeList2.size() > 1)
                {
                    break;
                }
            }

            int countTmp = 0;
            while (!converged and countTmp < totalItr)
            {
                memcpy(oldAB1, ab1, sizeof(double) * 2);
                memcpy(oldAB2, ab2, sizeof(double) * 2);

                //                printf("list1's size=%d, list2's size=%d\n", nodeList1.size(), nodeList2.size());
                if (!simpleLinearRegression(nodeList1, ab1))
                    break;
                if (!simpleLinearRegression(nodeList2, ab2))
                    break;
                if (utils::equalArrays(ab1, oldAB1, 2) && utils::equalArrays(ab2, oldAB2, 2))
                    converged = true;
                else
                {
                    nodeList1.clear();
                    nodeList2.clear();
                    for (int i = 0; i < _nodeList->size(); i++)
                    {
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
            if (SuperTAD::_VERBOSE_)
            {
                printf("finish convergence, %f, %f, %f, %f\n", ab1[0], ab1[1], ab2[0], ab2[1]);
                if (std::isnan(ab1[0]))
                { break; }
            }
            if (ab1[0] < ab2[0])
                trueNodeList = nodeList1;
            else
                trueNodeList = nodeList2;

            if (ab1[0] < 0 and ab2[0] < 0)
            {
                time++;
                for (int m = 0; m < trueNodeList.size(); m++)
                {
                    _scoreTable[trueNodeList[m].first]++;
                }
            }
        }

        double sum = 0;
        for (int i = 0; i < _nodeList->size(); i++)
        {
            _scoreTable[i] = _scoreTable[i] / 1000.0;
            sum += _scoreTable[i];
        }
        double mean = sum / _nodeList->size();

        for (int i = 0; i < _nodeList->size(); i++)
        {
            //            printf("%d, %f, %d \n", i, _scoreTable[i], _scoreTable[i] > mean);
            if (_scoreTable[i] > mean)
                _scoreTable[i] = 1;
            else
                _scoreTable[i] = 0;
        }
        return reinterpret_cast<int *>(_scoreTable);
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
        double sumX = 0;    // sum(Xi)
        double sumY = 0;    // sum(Yi)
        double sumXsquare = 0;  // sum(Xi*Xi)
        double sumXY = 0;   // sum(Xi*Yi)
        for (int i = 0; i < nodeList.size(); i++)
        {
            //            printf("node_info=%f\n", getY(*nodeList[i].second));
            sumX += getX(*nodeList[i].second);
            sumY += getY(*nodeList[i].second);
            sumXsquare += pow(getX(*nodeList[i].second), 2);
            sumXY += getX(*nodeList[i].second) * getY(*nodeList[i].second);
        }
        double temp = (nodeList.size() * sumXsquare - sumX * sumX);
        //        printf("simplelinearrepression, sumX=%f,sumY=%f, sumXsquare=%f,sumXY=%f, temp=%f\n", sumX, sumY, sumXsquare, sumXY, temp);
        if (temp)
        {
            ab[0] = (nodeList.size() * sumXY - sumX * sumY) / temp;
            ab[1] = (sumXsquare * sumY - sumX * sumXY) / temp;
        } else
        {
            ab[0] = 1;
            ab[1] = 0;
        }
        if (std::isnan(ab[0]) || std::isnan(ab[1]))
        {
            return false;
        }

        return true;
    }

}
