//
// Created by wang mengbo on 2019-09-01.
//

#include "detectorBinary.h"


namespace binary {

    Detector::Detector (Data &data)
    {
        _data = &data;
        _edgeCount = &data.edgeCount ();
//        if (_H==1)
//            _baseTable = new double *[_N];
//        else
        _table = new double **[_N];
        _minIndexArray = new int **[_N];
        _leftKArray = new int **[_N];
        for (int i = 0; i < _N; i++) {
            _table[i] = new double *[_N];
            _minIndexArray[i] = new int *[_N];
            _leftKArray[i] = new int *[_N];
            for (int j = 0; j < _N; j++) {
                _table[i][j] = new double[_K]{};
                _minIndexArray[i][j] = new int[_K]{};
                _leftKArray[i][j] = new int[_K]{};
            }
        }
        _boundary.reserve(_K);
        _binaryTree = new binary::Tree();

        int k = 1;
        for (int i = 0; i < _K; i++) {
            _kToIdx.emplace (k, i);
            k++;
        }
//        std::cout << "_kToIndex.size=" << _kToIdx.size () << "\n";
    }


    Detector::~Detector ()
    {
        delete _binaryTree;

//        if (_H==1) {
//            for (int i=0; i<_N; i++) {
//                delete _table[i];
//                delete _minIndexArray[i];
//                delete _leftKArray[i];
//            }
//        }
//        else {
        for (int i = 0; i < _N; i++) {
            for (int j = 0; j < _N; j++) {
                delete _table[i][j];
                delete _minIndexArray[i][j];
                delete _leftKArray[i][j];
            }
            delete _table[i];
            delete _minIndexArray[i];
            delete _leftKArray[i];
        }
//        }
        delete _table;
        delete _minIndexArray;
        delete _leftKArray;
    }


    void Detector::execute ()
    {
        fillTable();

        std::vector<double> sumOfEntropy;
        std::vector<double> sumOfLeaves;
        std::vector<utils::intDoublePair> normLeaves;
//        struct cmpNormLeaf {
//            bool operator() (const std::pair<int, double> &p1, const std::pair<int, double> &p2)
//            {
//                return p1.second < p2.second;
//            }
//        };

        int index = -1;
        if (_DETERMINE_K) {
            // determine K
            for (int num = 2; num < _K + 1; num++) {
                std::cout << "num=" << num << std::endl;
                double entropy = _table[0][_N - 1][indexK(num)];
                std::cout << "_table[0][" << _N - 1 << "][" << num << "]=" << entropy << std::endl;
                sumOfEntropy.emplace_back(entropy);

                backTrace(num);
                double leafSum = 0;
                for (int leaf = 0; leaf < _boundary.size(); leaf++) {
                    int currentStart = _boundary[leaf].first;
                    int currentEnd = _boundary[leaf].second;
                    leafSum += _data->getSE(currentStart, currentEnd, 2. * _data->_edgeSum);
                    leafSum += _table[currentStart][currentEnd][indexK(1)];
                }
                sumOfLeaves.emplace_back(leafSum);
                //    std::cout << "log2((double)_N / (double)num)=" << log2((double)_N / (double)num) << std::endl;
                double divisor = log2(_N / (double) num) + (_N * (num - 1) / (double) (num * (_N - 1))) * log2(num);
                std::cout << "leafSum=" << leafSum << ", divisor=" << divisor << std::endl;
                normLeaves.emplace_back(num, leafSum / divisor);
                std::cout << "--------\n";
            }
            for (int i = 0; i < normLeaves.size(); i++) {
                std::cout << normLeaves[i].first << ", " << normLeaves[i].second << std::endl;
            }
            sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePairBySecond);
            index = normLeaves[0].first;
            std::cout << "k chosen=" << index << std::endl;
        } else {
            int num = _K;
            std::cout << "num=" << num << std::endl;
            double entropy = _table[0][_N - 1][indexK(num)];
//            std::cout << "_table[0][" << _N - 1 << "][" << num << "]=" << entropy << std::endl;
            sumOfEntropy.emplace_back(entropy);

            backTrace(num);
            double leafSum = 0;
            for (int leaf = 0; leaf < _boundary.size(); leaf++) {
                int currentStart = _boundary[leaf].first;
                int currentEnd = _boundary[leaf].second;
                leafSum += _data->getSE(currentStart, currentEnd, 2. * _data->_edgeSum);
                leafSum += _table[currentStart][currentEnd][indexK(1)];
            }
            sumOfLeaves.emplace_back(leafSum);
            //    std::cout << "log2((double)_N / (double)num)=" << log2((double)_N / (double)num) << std::endl;
            double divisor = log2(_N / (double) num) + (_N * (num - 1) / (double) (num * (_N - 1))) * log2(num);
            std::cout << "leafSum=" << leafSum << ", divisor=" << divisor << std::endl;
            normLeaves.emplace_back(num, leafSum / divisor);
            std::cout << "--------\n";
        }
        index = normLeaves[0].first;
        backTrace(index, true);

        _nodeList = &_binaryTree->nodeList();
        _writer.writeTree(_INPUT+".original_boundaries.txt", *_nodeList);

        // filtering
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
    }


    void Detector::fillTable()
    {
//        if (_H==1) {
//            preSumGlog();
//            for (int i=0; i<_N; i++) {
//                double curVol = _edgeCount->coeff(i, 0) + _edgeCount->coeff(i, 0);
//                double binSum = getGlog(curVol) - _gLogSum[i];
//                _table[i][0] = _data->getSE(_edgeCount->coeff(i, 0), )
//
//            }
//        } else {
        for (int start = 0; start < _N; start++) {
            for (int end = start; end < _N; end++) {
                double currentVol = _data->getVol(start, end);
                for (int leaf = start; leaf < end + 1; leaf++) {
                    _table[start][end][indexK (1)] += _data->getSE(leaf, leaf, currentVol);
                }
            }
        }
        std::cout << "finishing filling the basic events in dp_table." << std::endl;
        for (int a = 2; a < _K + 1; a++) {
            for (int start = 0; start < _N; start++) {
                for (int end = start; end < _N; end++) {
                    double minTmp = std::numeric_limits<double>::infinity ();
                    int minIdx = 0;
                    int leftK = 0;
                    for (int binaryK = 1; binaryK < a; binaryK++)
                    {
                        for (int mid = start; mid < end; mid++)
                        {
                            double tmp = _table[start][mid][indexK(binaryK)] + _table[mid + 1][end][indexK(a - binaryK)];
                            double volParent;
                            volParent = _data->getVol(start, end);
                            tmp += _data->getSE(start, mid, volParent);
                            tmp += _data->getSE(mid+1, end, volParent);
                            if (tmp < minTmp)
                            {
                                minTmp = tmp;
                                minIdx = mid;
                                leftK = binaryK;
                            }
                        }
                    }
                    _minIndexArray[start][end][indexK(a)] = minIdx;
                    _table[start][end][indexK(a)] = minTmp;
                    _leftKArray[start][end][indexK(a)] = leftK;
                }
            }
            std::cout << "Finishing filling upper events where k = " << a << ", " << _table[0][_N-1][indexK(a)] << std::endl;
        }
//        }
    }


    void Detector::backTrace (int k, bool add)
    {
        init();

        binarySplit(0, _N - 1, k, add);
        _boundary.emplace_back(0, 0);
        sort(_boundary.begin(), _boundary.end(), utils::cmpBoundary);
        for (int i = 0; i < _boundary.size(); i++) {
            if (i == _boundary.size() - 1) {
                _boundary[i].second = _N - 1;
            } else {
                _boundary[i].second = _boundary[i + 1].first - 1;
            }
        }
    }


    void Detector::init()
    {
        _boundary.clear();
    }


    void Detector::binarySplit(int start, int end, int k, bool add)
    {
        //  std::cout << "----\n";
        //  std::cout << "k=" << k << std::endl;
        if (add)
            _binaryTree->add(start, end, indexK(k));

        if (k == 1)
            return;
        else {
//            std::cout << "start=" << start << ", end=" << end << ", k=" << k << std::endl;
//            printf("_minIndexArray[%d][%d][%d]=%e\n", start, end, k, _minIndexArray[start][end][k]);
            int midPos = _minIndexArray[start][end][indexK(k)];
//            int midPos = _minIndexArray[start][end][k];

//            printf("_leftKArray[%d][%d][%d]=%e\n", start, end, k, _leftKArray[start][end][k]);
            int leftK = _leftKArray[start][end][indexK(k)];
//            int leftK = _leftKArray[start][end][k];

//            std::cout << "midPos=" << midPos << ", leftK=" << leftK << std::endl;
            _boundary.emplace_back(midPos + 1, -1);
            binarySplit(start, midPos, leftK, add);

//            std::cout << "midPos+1=" << midPos+1 << ", k-leftK=" << k-leftK << std::endl;
            binarySplit(midPos + 1, end, k - leftK, add);
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


//    void Detector::preSumGlog()
//    {
//        _gLogSum[0] = getGlog(_edgeCount->coeff(0,0))
//        for (int i=1; i<_N; i++) {
//            _gLogSum[i] = _gLogSum[i-1] + getGlog(_edgeCount.coeff(i, i));
//        }
//    }
//
//
//    double Detector::getGlog(double binG)
//    {
//        if (binG==0)
//            return binG;
//        else
//            return binG * log2(binG);
//    }
}
