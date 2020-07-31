//
// Created by wang mengbo on 2019-09-03.
//

#include "detectorMulti.h"

namespace multi {

    Detector::Detector(Data &data)
    {
        _data = &data;
        _boundaries.reserve(_K_);
        _table = new double ****[_N_];
        _minIndexArray = new int ****[_N_];
        _leftKArray = new int ****[_N_];
        for (int s = 0; s < _N_; s++) {
            _table[s] = new double ***[_N_];
            _minIndexArray[s] = new int ***[_N_];
            _leftKArray[s] = new int ***[_N_];
            for (int e = s; e < _N_; e++) {
                _table[s][e] = new double **[_K_];
                _minIndexArray[s][e] = new int **[_K_];
                _leftKArray[s][e] = new int **[_K_];
                for (int k = 0; k < _K_; k++) {
                    _table[s][e][k] = new double *[_H_];
                    _minIndexArray[s][e][k] = new int *[_H_];
                    _leftKArray[s][e][k] = new int *[_H_];
                    for (int h = 0; h < _H_; h++) {
                        _table[s][e][k][h] = new double[_N_]{};
                        _minIndexArray[s][e][k][h] = new int[_N_]{};
                        _leftKArray[s][e][k][h] = new int[_N_]{};
                    }
                }
            }
        }
    }


    Detector::~Detector ()
    {
        for (int s = 0; s < _N_; s++) {
            for (int e = s; e < _N_; e++) {
                for (int k = 0; k < _K_; k++) {
                    for (int h = 0; h < _H_; h++) {
                        delete _table[s][e][k][h];
                        delete _minIndexArray[s][e][k][h];
                        delete _leftKArray[s][e][k][h];
                    }
                    delete _table[s][e][k];
                    delete _minIndexArray[s][e][k];
                    delete _leftKArray[s][e][k];
                }
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
    }

    void Detector::execute ()
    {
        std::clock_t tTmp;

        fillTable();

        std::vector<IntDoublePair> sumOfEntropy;
        std::vector<double> sumOfLeaves;
        std::vector<IntDoublePair> normLeaves;

        int kOpt = -1;
        double entropy;

        if (_DETERMINE_K_) {
            if (_VERBOSE_) {
                printf("start determine optimal K\n");
                tTmp = std::clock();
            }
            else
                printf("determine optimal K\n");

            for (int k=2; k<=_K_; k++) {

                if (_VERBOSE_)
                    printf("--------\nK=%d\n", k);


                entropy = _table[0][_N_ - 1][k - 1][_H_ - 1][_N_ - 1];
                if (_VERBOSE_)
                    printf("min structure entropy=%f\n", entropy);

//                sumOfEntropy.emplace_back(num, _table[0][_N_ - 1][indexK(num)][_H_ - 1][_N_ - 1]);
                sumOfEntropy.emplace_back(k, _table[0][_N_ - 1][k - 1][_H_ - 1][_N_ - 1]);

                backTrace(k, _H_);

                double leafSum = 0;
                int curS, curE;
                for (int leaf = 0; leaf < _boundaries.size(); leaf++) {
                    curS = _boundaries[leaf].first;
                    curE = _boundaries[leaf].second;
//                    leafSum += _data->getSE(curS, curE, 2 * _data->_edgeSum);
                    leafSum += _data->getSE(curS, curE, _data->_doubleEdgeSum);
                    leafSum += _table[curS][curE][0][0][curE];
                }
                sumOfLeaves.emplace_back(leafSum);
                double divisor =
                    log2(_N_ / (double) k) + (_N_ * (k - 1) / (double) (k * (_N_ - 1))) * log2((double) k);
                normLeaves.emplace_back(k, leafSum / divisor);
                fflush(stdout);

            }
            if (_VERBOSE_){
                printf("--------\n");
                printf("sumOfEntropy:\n");
                for (int i=0; i<sumOfEntropy.size(); i++)
                    printf("(%d, %f)\n", sumOfEntropy[i].first, sumOfEntropy[i].second);
                printf("normLeaves:\n");
                for (int i=0; i<normLeaves.size(); i++)
                    printf("(%d, %f)\n", normLeaves[i].first, normLeaves[i].second);
            }

            if (_H_ == 1 || _H_ == 2) {
                if (_optimalK_ < _K_)
                    kOpt = _optimalK_;
                else {
//                    sort(sumOfEntropy.begin(), sumOfEntropy.end(), utils::cmpIntDoublePairBySecond);
//                    kOpt = sumOfEntropy[0].first;
                    sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePairBySecond);
                    kOpt = normLeaves[0].first;
                }
            }
            else {
                sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePairBySecond);
                kOpt = normLeaves[0].first;
            }

            printf("optimal K is %d\n", kOpt);

            if (_VERBOSE_)
                printf("finish determine optimal K\n");

            if (_DEBUG_)
                printf("determining optimal K consumes %fs\n", (float)(std::clock()-tTmp)/CLOCKS_PER_SEC);

            backTrace(kOpt, _H_, true);

        }
        else {
            printf("K=%d\n", _K_);
//            sumOfEntropy.emplace_back(num, _table[0][_N_ - 1][indexK(num)][_H_ - 1][_N_ - 1]);
            sumOfEntropy.emplace_back(_K_, _table[0][_N_ - 1][_K_ - 1][_H_ - 1][_N_ - 1]);

            backTrace(_K_, _H_, true);

            double leafSum = 0;
            int curS, curE;
            for (int leaf = 0; leaf < _boundaries.size(); leaf++) {
                curS = _boundaries[leaf].first;
                curE = _boundaries[leaf].second;
//                leafSum += _data->getSE(currentStart, currentEnd, 2 * _data->_edgeSum);
                leafSum += _data->getSE(curS, curE, _data->_doubleEdgeSum);
                leafSum += _table[curS][curE][0][0][curE];
            }
            sumOfLeaves.emplace_back(leafSum);
            double divisor =
                log2(_N_ / (double) _K_) + (_N_ * (_K_ - 1) / (double) (_K_ * (_N_ - 1))) * log2((double) _K_);
            normLeaves.emplace_back(_K_, leafSum / divisor);
        }

        _nodeList = &_multiTree.nodeList();
        if (_VERBOSE_) {
            printf("nodes:");
            for (int i = 0; i < _nodeList->size(); i++) {
                printf("(%d, %d)", (*_nodeList)[i]->_val[0], (*_nodeList)[i]->_val[1]);
                if (i < _nodeList->size()-1)
                    printf(", ");
            }
            printf("\n");
        }

        _writer.writeTree(_OUTPUT_ + ".multi", *_nodeList);
    }


    void Detector::fillTable ()
    {
        std::clock_t t;
        if (_VERBOSE_) {
            if (_DEBUG_)
                t = std::clock();

            printf("start filling db table\n");
        }
        else
            printf("fill db table\n");

        for (int s = 0; s < _N_; s++) {
            for (int e = s; e < _N_; e++) {
                double currentVolume = _data->getVol(s, e);
                double binSum;
                if (s == 0) {
                    binSum = _data->getGtimesLogG(currentVolume) - _data->_sumOfGtimesLogG[e];
                } else {
                    binSum = _data->getGtimesLogG(currentVolume) - (_data->_sumOfGtimesLogG[e] - _data->_sumOfGtimesLogG[s - 1]);
                }
                _table[s][e][0][0][e] = binSum / (2. * _data->_edgeSum);
            }
        }

        for (int k = 1; k < _K_; k++) {
            for (int s = 0; s < _N_; s++) {
                for (int parentEnd = s; parentEnd < _N_; parentEnd++) {
                    for (int e = s; e < parentEnd + 1; e++) {

//                        if (e-s+1 < k)
//                            continue;

                        double parentVol = _data->getVol(s, parentEnd);
                        double minTmp = std::numeric_limits<double>::infinity();
                        int minIdx = 0;
                        for (int i = s; i < e; i++) {
                            double tmp;
                            if (k - 1 == 0) {
                                tmp = _table[s][i][k - 1][0][i];
                                tmp += _data->getSE(s, i, parentVol);
                            }
                            else {
                                tmp = _table[s][i][k - 1][0][parentEnd];
                            }

                            tmp += _data->getSE(i + 1, e, parentVol);

                            tmp += _table[i + 1][e][0][0][e];
                            if (tmp <= minTmp) {
                                minTmp = tmp;
                                minIdx = i;
                            }
                        }
                        _minIndexArray[s][e][k][0][parentEnd] = minIdx;
                        _table[s][e][k][0][parentEnd] = minTmp;
                    }
                }
            }
        }

        for (int h = 1; h < _H_; h++) {
            initH(h);
            for (int k = 2; k < _K_ + 1; k++) {
                for (int s = 0; s < _N_; s++) {
                    for (int parentEnd = s; parentEnd < _N_; parentEnd++) {
                        for (int e = s; e < parentEnd + 1; e++) {

//                            if (e-s+1 < k)
//                                continue;

                            double minTmp, currentVol, parentVol;
                            int minIdx, leftK;
                            if (e - s + 1 >= k) {
//                                minTmp = _table[start][end][indexK(cluster)][height - 1][end];
                                minTmp = _table[s][e][k - 1][h - 1][e];
                                parentVol = _data->getVol(s, parentEnd);
                                minTmp += _data->getSE(s, e, parentVol);
                                minIdx = s;
                                leftK = 0;
                                for (int kTmp = 1; kTmp < k; kTmp++) {
                                    for (int i = s; i < e; i++) {

//                                        if (i-s+1 < kTmp)
//                                            continue;
//
//                                        if (e-i+1 < k-kTmp)
//                                            continue;

                                        double tmp;
                                        if (kTmp == 1) {
//                                            tmp = _table[start][mid][indexK(binaryK)][height - 1][mid] +
//                                                  _table[mid + 1][end][indexK(cluster - binaryK)][height - 1][end];
                                            tmp = _table[s][i][kTmp - 1][h - 1][i] +
                                                  _table[i + 1][e][k - kTmp - 1][h - 1][e];
                                            tmp += _data->getSE(s, i, parentVol);
                                        }
                                        else {
//                                            tmp = _table[start][mid][indexK(binaryK)][height][parentEnd] +
//                                                  _table[mid + 1][end][indexK(cluster - binaryK)][height - 1][end];
                                            tmp = _table[s][i][kTmp - 1][h][parentEnd] +
                                                  _table[i + 1][e][k - kTmp - 1][h - 1][e];
                                        }

                                        tmp += _data->getSE(i + 1, e, parentVol);
                                        if (tmp <= minTmp) {
                                            minTmp = tmp;
                                            minIdx = i;
                                            leftK = kTmp;
                                        }
                                    }
                                }
                            }
                            else {
                                minTmp = std::numeric_limits<double>::infinity();
                                minIdx = 0;
                                leftK = 0;
                            }
//                            _minIndexArray[start][end][indexK(cluster)][height][parentEnd] = minIdx;
//                            _table[start][end][indexK(cluster)][height][parentEnd] = minTmp;
//                            _leftKArray[start][end][indexK(cluster)][height][parentEnd] = leftK;
                            _minIndexArray[s][e][k - 1][h][parentEnd] = minIdx;
                            _table[s][e][k - 1][h][parentEnd] = minTmp;
                            _leftKArray[s][e][k - 1][h][parentEnd] = leftK;
                        }
                    }
                }
                if ( _table[0][_N_ - 1][k-1][_H_ - 1][_N_ - 1] <  _table[0][_N_ - 1][k-2][_H_ - 1][_N_ - 1]){
                    _optimalK_ = k;
                    printf("--------\noptimalK=%d, table=%f\n", _optimalK_, _table[0][_N_ - 1][k-1][_H_ - 1][_N_ - 1]);
                }
                else{
                    if (_VERBOSE_)
                        printf("finish filling db table\n");

                    if (_DEBUG_)
                        printf("filling db table consumes %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);
                    return;
                }
            }
        }

        if (_VERBOSE_)
            printf("finish filling db table\n");

        if (_DEBUG_)
            printf("filling db table consumes %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);
    }


    void Detector::initH(int h)
    {
        for (int s = 0; s < _N_; s++) {
            for (int e = s; e < _N_; e++) {
                for (int k = 0; k < _N_; k++) {
                    _table[s][e][0][h][k] = std::numeric_limits<double>::infinity();
                }
            }
        }
    }


    void Detector::backTrace(int k, int h, bool add)
    {
        initK();

        multiSplit(0, _N_ - 1, k, h - 1, _N_ - 1, add);

        _boundaries.emplace_back(0, 0);

        sort(_boundaries.begin(), _boundaries.end(), utils::cmpBoundary);

        for (int i = 0; i < _boundaries.size(); i++) {
            if (i == _boundaries.size() - 1) {
                _boundaries[i].second = _N_ - 1;
            }
            else {
                _boundaries[i].second = _boundaries[i + 1].first - 1;
            }
        }

        if (_VERBOSE_) {
            printf("boundaries:");
            for (int i = 0; i < _boundaries.size(); i++) {
                printf("(%d, %d)", _boundaries[i].first, _boundaries[i].second);
                if (i < _boundaries.size()-1)
                    printf(", ");
            }
            printf("\n");
        }

        fflush(stdout);
    }


    void Detector::initK()
    {
        _boundaries.clear();
    }


    void Detector::multiSplit (int start, int end, int k, int h, int parentEnd, bool add)
    {
        if (k == 1) {
            if (add)
                _multiTree.insert(start, end);
            if (_DEBUG_) {
                printf("multisplit-------------%d %d %d %d k=1, parentEnd: %d %f\n",
                       start, end, k, h, parentEnd, _table[start][end][k-1][h][parentEnd]);
            }
            return;
        }
        else {
            if (h != 0) {
//                int leftK = _leftKArray[start][end][indexK(k)][h][parentEnd];
                int leftK = _leftKArray[start][end][k-1][h][parentEnd];
                if (leftK == 0) {
                    if (add)
                        _multiTree.insert(start, end);
                    if (_DEBUG_) {
                        printf("multisplit-------------%d %d %d %d leftK=1, parentEnd: %d %f\n",
                               start, end, k, h, parentEnd, _table[start][end][k - 1][h][parentEnd]);
                    }
                    multiSplit(start, end, k, h-1, end, add);
                }
                else {
//                    int midPos = _minIndexArray[start][end][indexK(k)][h][parentEnd];
                    int midPos = _minIndexArray[start][end][k-1][h][parentEnd];
                    _boundaries.emplace_back(midPos + 1, 0);
                    if (add)
                        _multiTree.insert(midPos+1, end);
                    if (_DEBUG_) {
                        printf("multisplit-------------%d %d %d %d h!=1, parentEnd: %d %f\n",
                               start, end, k, h, parentEnd, _table[start][end][k - 1][h][parentEnd]);
                    }
                    multiSplit(start, midPos, leftK, h, parentEnd, add);
                    multiSplit(midPos + 1, end, k-leftK, h-1, end, add);
                }
            }
            else {
//                int midPos = _minIndexArray[start][end][indexK(k)][h][parentEnd];
                int midPos = _minIndexArray[start][end][k-1][h][parentEnd];
                _boundaries.emplace_back(midPos + 1, 0);
                if (add)
                    _multiTree.insert(midPos+1, end);
                if (_DEBUG_) {
                    printf("multisplitmultisplit-------------%d %d %d %d h=1, parentEnd: %d %f\n",
                           start, end, k, h, parentEnd, _table[start][end][k - 1][h][parentEnd]);
                }
                multiSplit(start, midPos, k-1, h, parentEnd, add);
            }
        }
    }

}