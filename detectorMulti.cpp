//
// Created by wang mengbo on 2019-09-03.
//

#include "detectorMulti.h"
#include "data.h"
#include "params.h"

namespace multi {

    Detector::Detector(SuperTAD::Data &data)
    {
        _data = &data;
        _boundaries.reserve(SuperTAD::_K_);
        _table = new double ****[SuperTAD::_N_];
        _minIndexArray = new int ****[SuperTAD::_N_];
        _leftKArray = new int ****[SuperTAD::_N_];
        for (int s = 0; s < SuperTAD::_N_; s++) {
            _table[s] = new double ***[SuperTAD::_N_];
            _minIndexArray[s] = new int ***[SuperTAD::_N_];
            _leftKArray[s] = new int ***[SuperTAD::_N_];
            for (int e = s; e < SuperTAD::_N_; e++) {
                _table[s][e] = new double **[SuperTAD::_K_];
                _minIndexArray[s][e] = new int **[SuperTAD::_K_];
                _leftKArray[s][e] = new int **[SuperTAD::_K_];
                for (int k = 0; k < SuperTAD::_K_; k++) {
                    _table[s][e][k] = new double *[SuperTAD::_H_];
                    _minIndexArray[s][e][k] = new int *[SuperTAD::_H_];
                    _leftKArray[s][e][k] = new int *[SuperTAD::_H_];
                    for (int h = 0; h < SuperTAD::_H_; h++) {
                        _table[s][e][k][h] = new double[SuperTAD::_N_]{};
                        _minIndexArray[s][e][k][h] = new int[SuperTAD::_N_]{};
                        _leftKArray[s][e][k][h] = new int[SuperTAD::_N_]{};
                    }
                }
            }
        }
    }


    Detector::~Detector ()
    {
        for (int s = 0; s < SuperTAD::_N_; s++) {
            for (int e = s; e < SuperTAD::_N_; e++) {
                for (int k = 0; k < SuperTAD::_K_; k++) {
                    for (int h = 0; h < SuperTAD::_H_; h++) {
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

    std::vector<multi::TreeNode *>& Detector::execute()
    {
        std::clock_t tTmp;

        fillTable();

        std::vector<IntDoublePair> sumOfEntropy;
        std::vector<double> sumOfLeaves;
        std::vector<IntDoublePair> normLeaves;

        int kOpt = -1;
        double entropy;

        if (SuperTAD::_DETERMINE_K_) {
            if (SuperTAD::_VERBOSE_) {
                printf("start determine optimal K\n");
                tTmp = std::clock();
            }
            else
                printf("determine optimal K\n");

            for (int k=2; k <= SuperTAD::_optimalK_ + 1; k++) {

                if (SuperTAD::_VERBOSE_)
                    printf("--------\nK=%d\n", k);


                entropy = _table[0][SuperTAD::_N_ - 1][k - 1][SuperTAD::_H_ - 1][SuperTAD::_N_ - 1];
                if (SuperTAD::_VERBOSE_)
                    printf("min structure entropy=%f\n", entropy);

//                sumOfEntropy.emplace_back(num, _table[0][_N_ - 1][indexK(num)][_H_ - 1][_N_ - 1]);
                sumOfEntropy.emplace_back(k, _table[0][SuperTAD::_N_ - 1][k - 1][SuperTAD::_H_ - 1][SuperTAD::_N_ - 1]);

                backTrace(k, SuperTAD::_H_);

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
                    log2(SuperTAD::_N_ / (double) k) + (SuperTAD::_N_ * (k - 1) / (double) (k * (SuperTAD::_N_ - 1))) * log2((double) k);
                normLeaves.emplace_back(k, leafSum / divisor);
                fflush(stdout);

            }
            if (SuperTAD::_VERBOSE_){
                printf("--------\n");
                printf("sumOfEntropy:\n");
                for (int i=0; i<sumOfEntropy.size(); i++)
                    printf("(%d, %f)\n", sumOfEntropy[i].first, sumOfEntropy[i].second);
                printf("normLeaves:\n");
                for (int i=0; i<normLeaves.size(); i++)
                    printf("(%d, %f)\n", normLeaves[i].first, normLeaves[i].second);
            }

//            if (_H_ == 1 || _H_ == 2) {
//                if (_optimalK_ < _K_)
//                    kOpt = _optimalK_;
//                else {
//                    sort(sumOfEntropy.begin(), sumOfEntropy.end(), utils::cmpIntDoublePairBySecond);
//                    kOpt = sumOfEntropy[0].first;
//                    sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePairBySecond);
//                    kOpt = normLeaves[0].first;
//                }
//            }
//            else {
//                sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePairBySecond);
//                kOpt = normLeaves[0].first;
//            }
            kOpt = SuperTAD::_optimalK_;
            printf("optimal K is %d\n", kOpt);

            if (SuperTAD::_VERBOSE_)
                printf("finish determine optimal K\n");

            if (SuperTAD::_DEBUG_)
                printf("determining optimal K consumes %fs\n", (float)(std::clock()-tTmp)/CLOCKS_PER_SEC);

            backTrace(kOpt, SuperTAD::_H_, true);

        }
        else {
            printf("K=%d\n", SuperTAD::_K_);
//            sumOfEntropy.emplace_back(num, _table[0][_N_ - 1][indexK(num)][_H_ - 1][_N_ - 1]);
            sumOfEntropy.emplace_back(SuperTAD::_K_, _table[0][SuperTAD::_N_ - 1][SuperTAD::_K_ - 1][SuperTAD::_H_ - 1][
                SuperTAD::_N_ - 1]);

            backTrace(SuperTAD::_K_, SuperTAD::_H_, true);
        }

        _nodeList = &_multiTree.nodeList();
        if (SuperTAD::_VERBOSE_) {
            printf("nodes:");
            for (int i = 0; i < _nodeList->size(); i++) {
                printf("(%d, %d)", (*_nodeList)[i]->_val[0], (*_nodeList)[i]->_val[1]);
                if (i < _nodeList->size()-1)
                    printf(", ");
            }
            printf("\n");
        }

        if (!SuperTAD::_NO_OUTPUT_)
            _writer.writeTree(SuperTAD::_OUTPUT_ + ".multi", *_nodeList);

        return *_nodeList;
    }


    void Detector::fillTable ()
    {
        std::clock_t t;
        if (SuperTAD::_VERBOSE_) {
            if (SuperTAD::_DEBUG_)
                t = std::clock();

            printf("start filling db table\n");
        }
        else
            printf("fill db table\n");

        for (int s = 0; s < SuperTAD::_N_; s++) {
            for (int e = s; e < SuperTAD::_N_; e++) {
                double currentVolume = _data->getVol(s, e);
                double binSum;
                if (s == 0) {
                    binSum = _data->getGtimesLogG(currentVolume) - _data->_sumOfGtimesLogG[e];
                } else {
                    binSum = _data->getGtimesLogG(currentVolume) - (_data->_sumOfGtimesLogG[e] - _data->_sumOfGtimesLogG[s - 1]);
                }
                _table[s][e][0][0][e] = binSum / _data->_doubleEdgeSum;
            }
        }


        for (int k = 1; k < SuperTAD::_K_; k++) {
            for (int s = 0; s < SuperTAD::_N_; s++) {
                for (int parentEnd = s; parentEnd < SuperTAD::_N_; parentEnd++) {
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

        for (int h = 1; h < SuperTAD::_H_; h++) {
            initH(h);
            for (int k = 2; k < SuperTAD::_K_ + 1; k++) {
                for (int s = 0; s < SuperTAD::_N_; s++) {
                    for (int parentEnd = s; parentEnd < SuperTAD::_N_; parentEnd++) {
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
                            _minIndexArray[s][e][k - 1][h][parentEnd] = minIdx;
                            _table[s][e][k - 1][h][parentEnd] = minTmp;
                            _leftKArray[s][e][k - 1][h][parentEnd] = leftK;
                        }
                    }
                }
//                printf("h=%d, k=%d, table=%f\n", h, k, _table[0][_N_-1][k-1][h][_N_-1]);

                if (h == SuperTAD::_H_ - 1 and SuperTAD::_DETERMINE_K_ ) {
                    if (_table[0][SuperTAD::_N_ - 1][k - 1][h][SuperTAD::_N_ - 1] < _table[0][SuperTAD::_N_ - 1][k - 2][h][
                        SuperTAD::_N_ - 1]) {
                        SuperTAD::_optimalK_ = k;
                        printf("--------\noptimalK=%d, table=%f\n", SuperTAD::_optimalK_, _table[0][SuperTAD::_N_ - 1][k - 1][
                            SuperTAD::_H_ - 1][SuperTAD::_N_ - 1]);
                    } else {
                        printf("--------\nno more optimalK=%d, table=%f\n", SuperTAD::_optimalK_, _table[0][
                            SuperTAD::_N_ - 1][k - 1][SuperTAD::_H_ - 1][SuperTAD::_N_ - 1]);
                        if (SuperTAD::_VERBOSE_)
                            printf("finish filling db table\n");

                        if (SuperTAD::_DEBUG_)
                            printf("filling db table consumes %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);
                        return;
                    }
                }
            }
        }

        if (SuperTAD::_VERBOSE_)
            printf("finish filling db table\n");

        if (SuperTAD::_DEBUG_)
            printf("filling db table consumes %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);
    }


    void Detector::initH(int h)
    {
        for (int s = 0; s < SuperTAD::_N_; s++) {
            for (int e = s; e < SuperTAD::_N_; e++) {
                for (int k = 0; k < SuperTAD::_N_; k++) {
                    _table[s][e][0][h][k] = std::numeric_limits<double>::infinity();
                }
            }
        }
    }


    void Detector::backTrace(int k, int h, bool add)
    {
        initK();

        multiSplit(0, SuperTAD::_N_ - 1, k, h - 1, SuperTAD::_N_ - 1, add);

        _boundaries.emplace_back(0, 0);

        sort(_boundaries.begin(), _boundaries.end(), utils::cmpBoundary);

        for (int i = 0; i < _boundaries.size(); i++) {
            if (i == _boundaries.size() - 1) {
                _boundaries[i].second = SuperTAD::_N_ - 1;
            }
            else {
                _boundaries[i].second = _boundaries[i + 1].first - 1;
            }
        }

        if (SuperTAD::_VERBOSE_) {
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
            if (SuperTAD::_DEBUG_) {
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
                    if (SuperTAD::_DEBUG_) {
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
                    if (SuperTAD::_DEBUG_) {
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
                if (SuperTAD::_DEBUG_) {
                    printf("multisplitmultisplit-------------%d %d %d %d h=1, parentEnd: %d %f\n",
                           start, end, k, h, parentEnd, _table[start][end][k - 1][h][parentEnd]);
                }
                multiSplit(start, midPos, k-1, h, parentEnd, add);
            }
        }
    }

}