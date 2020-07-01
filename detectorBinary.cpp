//
// Created by wang mengbo on 2019-09-01.
//

#include "detectorBinary.h"


namespace binary {

    Detector::Detector (Data &data)
    {
        _data = &data;
        _edgeCount = &data.edgeCount ();
        _table = new double **[_N];
        _minIndexArray = new int **[_N];
        _leftKArray = new int **[_N];
        for (int s=0; s<_N; s++) {
            _table[s] = new double *[_N];
            _minIndexArray[s] = new int *[_N];
            _leftKArray[s] = new int *[_N];
            int k;
            for (int e=s; e<_N; e++) {
                k = (e-s+1 < _K ? e-s+1 : _K);
                _table[s][e] = new double[k]{};
                _minIndexArray[s][e] = new int[k]{};
                _leftKArray[s][e] = new int[k]{};
            }
        }
        _boundary.reserve(_K);
        _binaryTree = new binary::Tree();

//        int k = 1;
//        for (int i = 0; i < _K; i++)
//            _kToIdx.emplace (k++, i);
//        std::cout << "_kToIndex.size=" << _kToIdx.size () << "\n";

        _numBins = new int(0);
        _kTmpIdx = new int(0);
        _kMinusTmpIdx = new int(0);
    }


    Detector::~Detector ()
    {
        delete _binaryTree;
        for (int s = 0; s < _N; s++) {
            for (int e = s; e < _N; e++) {
                delete _table[s][e];
                delete _minIndexArray[s][e];
                delete _leftKArray[s][e];
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
        delete _kMinusTmpIdx;
    }


    void Detector::execute ()
    {
        fillTable();

        std::vector<double> sumOfEntropy;
        std::vector<double> sumOfLeaves;
        std::vector<utils::intDoublePair> normLeaves;

        int index = -1;
        if (_DETERMINE_K) {
            double entropy;

            // determine K
            for (int kTmp = 2; kTmp <= _K; kTmp++) {
                std::cout << "kTmp=" << kTmp << std::endl;
                indexKtmp(kTmp);

                entropy = _table[0][_N - 1][*_kTmpIdx];
                printf("_table[0][_N-1][%d]=%f\n", *_kTmpIdx, entropy);
                sumOfEntropy.emplace_back(entropy);

                backTrace(kTmp);
                double leafSum = 0;
//                int currentS, currentE;
                for (int leaf = 0; leaf < _boundary.size(); leaf++) {
//                    currentS = _boundary[leaf].first;
//                    currentE = _boundary[leaf].second;
                    leafSum += _data->getSE(_boundary[leaf].first, _boundary[leaf].second, 2.*_data->_edgeSum);
                    leafSum += _table[_boundary[leaf].first][_boundary[leaf].second][0];
                }
                sumOfLeaves.emplace_back(leafSum);
                //    std::cout << "log2((double)_N / (double)num)=" << log2((double)_N / (double)num) << std::endl;
                double divisor = log2(_N / (double) kTmp) + (_N * (kTmp - 1) / (double) (kTmp * (_N - 1))) * log2(kTmp);
                std::cout << "leafSum=" << leafSum << ", divisor=" << divisor << std::endl;
                normLeaves.emplace_back(kTmp, leafSum / divisor);
                std::cout << "========\n\n";
            }
            for (int i = 0; i < normLeaves.size(); i++) {
                std::cout << normLeaves[i].first << ", " << normLeaves[i].second << std::endl;
            }
            sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePairBySecond);
            index = normLeaves[0].first;
            std::cout << "k chosen=" << index << std::endl;

        } else {

            int kTmp = _K;
            std::cout << "kTmp=" << kTmp << std::endl;
            indexKtmp(kTmp);

            double entropy = _table[0][_N - 1][*_kTmpIdx];
//            std::cout << "_table[0][" << _N - 1 << "][" << num << "]=" << entropy << std::endl;
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
            //    std::cout << "log2((double)_N / (double)num)=" << log2((double)_N / (double)num) << std::endl;
            double divisor = log2(_N / (double) kTmp) + (_N * (kTmp - 1) / (double) (kTmp * (_N - 1))) * log2(kTmp);
            std::cout << "leafSum=" << leafSum << ", divisor=" << divisor << std::endl;
            normLeaves.emplace_back(kTmp, leafSum/divisor);
            std::cout << "========\n\n";
        }

        index = normLeaves[0].first;
        backTrace(index, true);

        _nodeList = &_binaryTree->nodeList();
        _writer.writeTree(_INPUT+".original_boundaries.txt", *_nodeList);

        // filtering
        std::clock_t t = std::clock();
        if (_FILTERING) {
            calculateD (_binaryTree->root());
            calculateDensity(_binaryTree->root());
            filterNodes();
            std::vector<binary::TreeNode *> trueNodes;
            for (auto it = _trueNodeList.begin(); it != _trueNodeList.end(); it++) {
                trueNodes.emplace_back((*it));
            }
            _writer.writeTree(_INPUT+".filter_boundaries.txt", trueNodes);
        }
        float time = (float)(std::clock() - t) / CLOCKS_PER_SEC;
        printf("filtering consumes %fs\n", time);
    }


    void Detector::fillTable()
    {
//        str_2_i2dMap map;
        if (_VERBOSE)
            printf("filling dp table\n");

        // process k=1
        if (_VERBOSE)
            printf("filling base cases\n");
        for (int s = 0; s < _N; s++) {
            for (int e = s; e < _N; e++) {
                double currentVol = _data->getVol(s, e);
                double binSum;
                if (s == 0) {
                    binSum = _data->getGtimesLogG(currentVol) - _data->_sumOfGtimesLogG[e];
                } else {
                    binSum = _data->getGtimesLogG(currentVol) - (_data->_sumOfGtimesLogG[e] - _data->_sumOfGtimesLogG[s - 1]);
                }
//                indexK(1);
//                printf("s=%d, e=%d\t", s, e);
//                printf("_table[%d][%d][0]=%f\n", s, e, _table[s][e][0]);
                _table[s][e][0] = binSum / (2. * _data->_edgeSum);
            }
        }
        if (_VERBOSE)
            printf("finish filling base case where k=1, _table[0][%d][0]=%f\n", _N-1, _table[0][_N-1][0]);

        // process k>=2
        double minSe, seTmp, parentVol;
        int minIdx, leftK, kIdx;
        bool breakFlag = false;
        for (int k = 2; k <= _K; k++) {
            if (breakFlag)
                break;
            indexK(k, kIdx);

            for (int s = 0; s < _N; s++) {
                if (breakFlag)
                    break;

                for (int e = s; e < _N; e++) {
                    if (breakFlag)
                        break;

                    minSe = std::numeric_limits<double>::infinity();
                    minIdx = 0;
                    leftK = 0;

                    // skip cases when #numBins<#numLeves
                    numBins(s, e);
                    if (*_numBins < k) {
//                        if (_VERBOSE)
//                            printf("s=%d, e=%d, k=%d, #numBins(%d)<#numLeaves(%d); meaningless; skip\n", s, e, k, *_numBins, k);
                        continue;
                    }

                    // find min{S(s,i,k_t)+H_l(s,e,i)+S(i+1,e,k-k_t)+H_r(s,e,i)}
                    for (int kTmp=1; kTmp<k; kTmp++)
                    {
                        indexKtmp(kTmp);
                        indexK(k-kTmp, *_kMinusTmpIdx);

//                        std::string key = std::to_string(s) + "_" + std::to_string(e) + "_" + std::to_string(k) + "_" + std::to_string(kTmp);
//                        i2dMap map2;
                        for (int i=s; i<e; i++)
                        {
                            if (i-s+1 < kTmp || e-(i+1)+1 < k-kTmp) {
//                                if (_VERBOSE) {
//                                    numBins(s, i);
//                                    printf("s=%d, i=%d, #numBins=%d, kTmp=%d;%s\n",
//                                        s, i, *_numBins, kTmp, (*_numBins<kTmp?"#numBins<#numLeaves; skip":""));
//                                    numBins(i+1, e);
//                                    printf("i+1=%d, e=%d, #numBins=%d, k-kTmp=%d;%s\n",
//                                        i+1, e, *_numBins, k-kTmp, (*_numBins<k-kTmp?"#numBins<#numLeaves; skip":""));
//                                }
                                continue;
                            }

                            seTmp = _table[s][i][*_kTmpIdx];
                            seTmp += _table[i+1][e][*_kMinusTmpIdx];
                            parentVol = _data->getVol(s, e);
                            seTmp += _data->getSE(s, i, parentVol);
                            seTmp += _data->getSE(i + 1, e, parentVol);
                            if (seTmp < minSe)
                            {
                                minSe = seTmp;
                                minIdx = i;
                                leftK = kTmp;
                            }
//                            map2.emplace(mid, tmp);
                        }
//                        map.emplace(key, map2);
                    }
                    _minIndexArray[s][e][kIdx] = minIdx;
                    _table[s][e][kIdx] = minSe;
                    _leftKArray[s][e][kIdx] = leftK;

                    if (k==_K && s==0 && e==_N-1) {
                        if (_VERBOSE)
                            printf("already obtain optimal solution for minS(1, %d, %d); stop here\n", _N, _K);
                        breakFlag=true;
                    }
                }
            }
            printf("finishing filling upper events where k=%d, _table[0][%d][%d]=%f\n", k, _N-1, kIdx, _table[0][_N-1][kIdx]);
        }
//        if (_TMP_PATH_!="") {
//            Writer::writeListOfCoordinates(map, _TMP_PATH_);
//        }
    }


    void Detector::backTrace (int k, bool add)
    {
        printf("start backtrace\n");
        init();

        binarySplit(0, _N-1, k, add);

        _boundary.emplace_back(0, 0);

        sort(_boundary.begin(), _boundary.end(), utils::cmpBoundary);
        for (int i = 0; i < _boundary.size(); i++) {
            if (i == _boundary.size() - 1) {
                _boundary[i].second = _N - 1;
            } else {
                _boundary[i].second = _boundary[i + 1].first - 1;
            }
        }
        printf("finish backtrace\n");
    }


    void Detector::init()
    {
        _boundary.clear();
    }


    void Detector::binarySplit(int s, int e, int k, bool add, int lv)
    {
        std::string gap = "";
        for (int i=0; i<lv; i++)
            gap += " ";

        printf("%s----binarySplit k=%d----\n", gap.c_str(), k);

        indexKtmp(k);

        if (add) {
            _binaryTree->add(s, e, *_kTmpIdx);
        }

        if (k == 1) {
            printf("%sk==1, return\n\n", gap.c_str());
            return;
        }
        else {
            printf("%s_minIndexArray[%d][%d][%d]=%f\n", gap.c_str(), s, e, *_kTmpIdx, _minIndexArray[s][e][*_kTmpIdx]);
            printf("%s_leftKArray[%d][%d][%d]=%f\n", gap.c_str(), s, e, *_kTmpIdx, _leftKArray[s][e][*_kTmpIdx]);

            int i = _minIndexArray[s][e][*_kTmpIdx];

            int kTmp = _leftKArray[s][e][*_kTmpIdx];

            printf("%si=%d, kTmp=%d\n", gap.c_str(), i, kTmp);

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
//      std::cout << "node=" << node << std::endl;
//      std::cout << "node._left=" << *node._left << std::endl;
            int leftStart = node._left->_val[0];
            int leftEnd = node._left->_val[1];
//      std::cout << "node._right=" << *node._right << std::endl;
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
            calculateDensity (*node._left);
        if (node._right != NULL)
            calculateDensity (*node._right);
    }


    double Detector::minusParent (double d, binary::TreeNode &node)
    {
        binary::TreeNode *currentNode = &node;
        while (!(*currentNode == _binaryTree->root ())) {
            currentNode = currentNode->_parent;
            d -= currentNode->_info;
        }
        return d;
    }


    void Detector::filterNodes ()
    {
        Eigen::MatrixXi scoreMat (_nodeList->size (), _nodeList->size ());
        scoreMat.setZero ();

        std::vector<std::pair<int, binary::TreeNode *>> nodeList1, nodeList2, trueNodeList;
        nodeList1.reserve (_nodeList->size () / 2);
        nodeList2.reserve (_nodeList->size () / 2);
        trueNodeList.reserve (_nodeList->size () / 2);

        double ab1[2]{};
        double ab2[2]{};

        int totalItr = 1000;
        int threshold = 900;
        int time = 0;
        //  int count = 0;
        //  int countElse = 0;
        while (time < totalItr) {
            //    std::cout << "--------\ncount=" << ++count << std::endl;
            //    std::cout << "--------\ntime=" << time << std::endl;
            //    std::cout << "--------\ncountElse=" << countElse << std::endl;
            _trueNodeList.clear ();
            double oldAB1[2]{};
            double oldAB2[2]{};
            bool converged = false;

            // init filter
            while (true) {
                nodeList1.clear ();
                nodeList2.clear ();
                for (int i = 0; i < _nodeList->size (); i++) {
                    if (utils::randInt (0, 2) == 0)
                        nodeList1.emplace_back (i, (*_nodeList)[i]);
                    else
                        nodeList2.emplace_back (i, (*_nodeList)[i]);
                }

                if (nodeList1.size () > 1 && nodeList2.size () > 1) {
                    break;
                }
            }

            int countTmp = 0;
            while (!converged and countTmp < totalItr) {
                countTmp++;
                //      std::cout << "countTmp=" << ++countTmp << std::endl;

                utils::copyDoubleArray (ab1, oldAB1, 2);
                utils::copyDoubleArray (ab2, oldAB2, 2);

                if (!simpleLinearRegression (nodeList1, ab1))
                    break;
                //      std::cout << "ab1=(" << ab1[0] << ", " << ab1[1] << ")\n";

                if (!simpleLinearRegression (nodeList2, ab2))
                    break;
                //      std::cout << "ab2=(" << ab2[0] << ", " << ab2[1] << ")\n";

                if (utils::equalDoubleArrays(ab1, oldAB1, 2) && utils::equalDoubleArrays(ab2, oldAB2, 2))
                    converged = true;
                else {
                    nodeList1.clear();
                    nodeList2.clear();
                    for (int i = 0; i < _nodeList->size(); i++) {
                        binary::TreeNode *nodeTmp = (*_nodeList)[i];
                        double dist1 = pow(ab1[0] * getX (*nodeTmp) + ab1[1] - nodeTmp->_info, 2);
                        double dist2 = pow(ab2[0] * getX (*nodeTmp) + ab2[1] - nodeTmp->_info, 2);
                        if (dist1 < dist2)
                            nodeList1.emplace_back(i, nodeTmp);
                        else
                            nodeList2.emplace_back(i, nodeTmp);
                    }
                }
            }

//            if (!converged)
//                continue;

            if (ab1[0] < ab2[0])
                trueNodeList = nodeList1;
            else
                trueNodeList = nodeList2;

            if (abs (ab1[0] - ab2[0]) > .5) {
                time++;
                for (int m = 0; m < trueNodeList.size () - 1; m++) {
                    for (int n = m + 1; n < trueNodeList.size (); n++) {
                        scoreMat (trueNodeList[m].first, trueNodeList[n].first)++;
                    }
                }
            }
            //    else
            //      countElse++;
        }

        for (int i = 0; i < _nodeList->size (); i++) {
            for (int j = i + 1; j < _nodeList->size (); j++) {
                if (scoreMat.coeff (i, j) > threshold) {
                    _trueNodeList.emplace ((*_nodeList)[i]);
                    _trueNodeList.emplace ((*_nodeList)[j]);
                }
            }
        }
        return;
    }


    double Detector::getX (binary::TreeNode &node)
    {
        double size = node._val[1] - node._val[0] + 1;
        return 1. / 3. * (size + 1);
    }


    double Detector::getY (binary::TreeNode &node)
    {
        return node._info;
    }


    bool Detector::simpleLinearRegression (std::vector<std::pair<int, binary::TreeNode *>> &nodeList, double ab[])
    {
        double sumX = 0;
        double sumY = 0;
        for (int i = 0; i < nodeList.size (); i++) {
            //    std::cout << "x=" << getX (*nodeList[i].second);
            sumX += getX (*nodeList[i].second);
            //    std::cout << "_D=" << nodeList[i].second->_D << ", y=" << getY (*nodeList[i].second) << std::endl;
            sumY += getY (*nodeList[i].second);
        }
        double meanX = sumX / nodeList.size ();
        double meanY = sumY / nodeList.size ();
        double covXY = 0;
        double varX = 0;
        //  std::cout << "meanX=" << meanX << ", meanY=" << meanY << std::endl;
        for (int i = 0; i < nodeList.size (); i++) {
            covXY += (getX (*nodeList[i].second) - meanX) * (getY (*nodeList[i].second) - meanY);
            varX += pow (getX (*nodeList[i].second) - meanX, 2);
        }
        //  std::cout << "covXY=" << covXY << ", varX=" << varX << ", meanY=" << meanY << ", meanX=" << meanX << std::endl;
        ab[0] = covXY / varX;
        ab[1] = meanY - ab[0] * meanX;
        if (std::isnan (ab[0]) || std::isnan (ab[1]))
            return false;

        return true;
    }


    bool Detector::meaningfulComb(int s, int e, int k)
    {
        return e-s+1 >= k;
    }
}
