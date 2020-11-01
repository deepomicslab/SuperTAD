//
// Created by mengbowang on 4/29/2020.
//

#include "detectorH.h"
#include "data.h"
#include "params.h"
#include "inputAndOutput.h"

namespace multi {

    DetectorH1::DetectorH1(SuperTAD::Data &data) {
        _data = &data;
//        _edgeCount = &data.edgeCount();
//        for (int i=0; i < SuperTAD::_K_; i++) {
//            _kToIdx.emplace(k++, i);
//        }
        SuperTAD::_K_ = SuperTAD::_N_/SuperTAD::_MinSize_;
        printf("set max K to be %d\n", SuperTAD::_K_);
        _table = new double *[SuperTAD::_N_];
        _minIndexArray = new int *[SuperTAD::_N_];
        for (int i=0; i < SuperTAD::_N_; i++) {
            _table[i] = new double [SuperTAD::_N_]{};
            _minIndexArray[i] = new int [SuperTAD::_N_]{};
        }
    }


    DetectorH1::~DetectorH1() {
        for (int i=0; i < SuperTAD::_N_; i++) {
            delete _table[i];
            delete _minIndexArray[i];
        }
        delete _table;
        delete _minIndexArray;
    }


    std::vector<Boundary> DetectorH1::execute(int h) {
        _table[0][0] = _data->getSE(0, 0, _data->_doubleEdgeSum);

        double currentVol, parentVol, binSum;

        if (SuperTAD::_VERBOSE_)
            printf("start k=0\n");

        for (int i=1; i < SuperTAD::_N_; i++) {
            parentVol = _data->_doubleEdgeSum;
            currentVol = _data->getVol(0, i);
            binSum = _data->getGtimesLogG(currentVol) - _data->_sumOfGtimesLogG[i];
            _table[i][0] = _data->getSE(0, i, parentVol, currentVol) + binSum/_data->_doubleEdgeSum;
//            printf("base case: i=%d, parentVol=%f, currentVol=%f,table=%f\n", i, parentVol, currentVol, _table[i][0]);
        }

        if (SuperTAD::_VERBOSE_)
            printf("finish k=0, se=%f\n", _table[SuperTAD::_N_ - 1][0]);

        if (SuperTAD::_VERBOSE_)
            printf("start calculating h=1\n");

        double minSE, tmpSE;
        int minIdx;
        double sumOfLeavesTmp = _table[SuperTAD::_N_ - 1][0];
        for (int a=1; a < SuperTAD::_K_; a++) {
            if (SuperTAD::_VERBOSE_)
                printf("start k=%d\n", a);
            for (int b=a; b < SuperTAD::_N_; b++) {
                minSE = std::numeric_limits<double>::infinity();
                minIdx = 0;
                for (int i=a-1; i<b; i++) {
                    tmpSE = _table[i][a - 1];
                    if (i+1==b)
                        tmpSE += _data->getSE(b, b, _data->_doubleEdgeSum);
                    else {
                        parentVol = _data->_doubleEdgeSum;
                        currentVol = _data->getVol(i+1, b);
                        tmpSE += _data->getSE(i + 1, b, parentVol, currentVol);
                        binSum = _data->getGtimesLogG(currentVol) - (_data->_sumOfGtimesLogG[b] - _data->_sumOfGtimesLogG[i]);
                        tmpSE += binSum / _data->_doubleEdgeSum;
                    }
                    if (tmpSE < minSE) {
                        minSE = tmpSE;
                        minIdx = i;
                    }
                }
                _minIndexArray[b][a] = minIdx;
                _table[b][a] = minSE;
            }
            if (SuperTAD::_VERBOSE_)
                printf("finish k=%d, structure entropy=%f\n", a, minSE);
            if (SuperTAD::_DETERMINE_K_) {
                if (_table[SuperTAD::_N_ - 1][a] < sumOfLeavesTmp){
                    SuperTAD::_optimalK_ = a;
                    sumOfLeavesTmp = minSE;
                }
                else
                    break;
            }
        }

        if (SuperTAD::_VERBOSE_)
            printf("finish calculating h=1\n");

        if (SuperTAD::_DETERMINE_K_) {
//            printf("optimalK=%d, K=%d, _k=%d\n", SuperTAD::_optimalK_ , SuperTAD::_K_, _k);
            if (SuperTAD::_optimalK_ < SuperTAD::_K_)
                _k = SuperTAD::_optimalK_ + 1;
            printf("determine k to be %d\n", _k);
        }

        DetectorH1::backTrace();
        if (h==-1)
            _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".multi2D", _boundaries);
        return _boundaries;
    }


    /*
        Back tracing with max_index_array to find the optimal route.
        param k: num of clustering
        return the positions in the optimal route.
     */
    void DetectorH1::backTrace() {
        printf("#data points: %d \n", SuperTAD::_N_);
        int *boundaries = new int [_k];
        boundaries[_k - 1] = SuperTAD::_N_ - 1;
        for (int i=_k-2; i>-1; i--) {
            int t1 = boundaries[i + 1];
            boundaries[i] = _minIndexArray[t1][i + 1];
        }

        for (int i=0; i<_k; i++) {
            if (i==0)
                _boundaries.emplace_back(1, boundaries[i]+1);
            else
                _boundaries.emplace_back(boundaries[i - 1] + 2, boundaries[i]+1);
        }
//        delete &boundaries;
//        for (int i = 0; i < _boundaries.size(); i++) {
//            std::cout << "boundary[" << i << "]=(" << _boundaries[i].first << ", " << _boundaries[i].second << ")\n";
//        }
    }

    Merge::Merge(SuperTAD::Data &data, std::vector<Boundary> &_preBoundList) {
        _data = &data;
        _preBoundaries = _preBoundList;
        N = (int) _preBoundaries.size();
        _table = new double *[N];
        _minIndexArray = new int *[N];
        for (int i=0; i < N; i++) {
            _table[i] = new double [N]{};
            _minIndexArray[i] = new int [N]{};
        }
    }

    Merge::~Merge() {
        for (int i=0; i < N; i++) {
            delete _table[i];
            delete _minIndexArray[i];
        }
        delete _table;
        delete _minIndexArray;
    }

    std::vector<Boundary> Merge::execute(int h) {
        double currentVol, binSum;
        int start, end;
        // record se of each Pre-node
        for (std::vector<Boundary>::iterator it = _preBoundaries.begin(); it != _preBoundaries.end(); ++it) {
            start = it->first-1;
            end = it->second-1;
            currentVol = _data->getVol(start, end);
            binSum = _data->getGtimesLogG(currentVol) - (_data->_sumOfGtimesLogG[end] - _data->_sumOfGtimesLogG[start-1]);
            _prenodeSE.emplace_back(binSum/_data->_doubleEdgeSum);
//            printf("start=%d, end=%d, currentVol=%f, prenodeSE=%f\n", start, end, currentVol, binSum/_data->_doubleEdgeSum);
        }
        if (SuperTAD::_VERBOSE_) printf("Start calculating base case. #node=%d\n", N);
        // base case
        int nodeStart, nodeEnd;
        double nodeVol;
        for (int i=0; i < N; i++) {
            end = _preBoundaries[i].second- 1;
            currentVol = _data->getVol(0, end);
//            printf("end=%d, currentVol=%f\n", end, currentVol);
            _table[i][0] = 0;
            for (int j=0; j < i+1; j++) {
                _table[i][0] += _prenodeSE[j];
                nodeStart = _preBoundaries[j].first - 1;
                nodeEnd = _preBoundaries[j].second - 1;
                nodeVol = _data->getVol(nodeStart, nodeEnd);
                _table[i][0] += _data->getSE(nodeStart, nodeEnd, currentVol, nodeVol);
//                printf("nodestart=%d, nodeend=%d, nodevol=%f, _table=%f\n", nodeStart, nodeEnd, nodeVol, _table[i][0]);
            }
            _table[i][0] += _data->getSE(0, end, _data->_doubleEdgeSum, currentVol);
        }
        // upper case
        if (SuperTAD::_VERBOSE_) printf("finish k = 0, se=%f\nStart caluclating upper case.\n", _table[N - 1][0]);
        double minSE, tmpSE;
        int minIdx;
        double sumOfLeavesTmp = _table[N-1][0];
        for (int a=1; a < N; a++) {
            for (int b=a; b < N; b++) {
                minSE = std::numeric_limits<double>::infinity();
                minIdx = 0;
                for (int i=a-1; i<b; i++) {
                    tmpSE = _table[i][a - 1];
                    if (i+1==b) {
                        tmpSE += _prenodeSE[b];
                        tmpSE += _data->getSE(_preBoundaries[b].first - 1, _preBoundaries[b].second - 1,
                                              _data->_doubleEdgeSum);
                    }
                    else {
                        start = _preBoundaries[i+1].first-1;
                        end = _preBoundaries[b].second-1;
                        currentVol = _data->getVol(start, end);
                        for (int node=i+1; node<b+1; node++) {
                            nodeStart = _preBoundaries[node].first-1;
                            nodeEnd = _preBoundaries[node].second-1;
                            tmpSE += _prenodeSE[node];
                            tmpSE += _data->getSE(nodeStart, nodeEnd, currentVol);
                        }
                        tmpSE += _data->getSE(start, end, _data->_doubleEdgeSum, currentVol);
                    }
                    if (tmpSE < minSE) {
                        minSE = tmpSE;
                        minIdx = i;
                    }
                }
//                printf("a=%d, b=%d, minIdx=%d, minSE=%f\n", a, b, minIdx, minSE);
                _minIndexArray[b][a] = minIdx;
                _table[b][a] = minSE;
            }
            if (SuperTAD::_VERBOSE_)
                printf("finish k = %d, structure entropy=%f\n", a, minSE);
            if (SuperTAD::_DETERMINE_K_) {
                if (_table[N-1][a] < sumOfLeavesTmp){
                    SuperTAD::_optimalK_ = a;
                    sumOfLeavesTmp = minSE;
                }
                else
                    break;
            }
        }

        if (SuperTAD::_VERBOSE_)
            printf("finish calculating h=1\n");

        if (SuperTAD::_DETERMINE_K_) {
            if (SuperTAD::_optimalK_ < SuperTAD::_K_)
                _k = SuperTAD::_optimalK_ + 1;
            printf("determine k to be %d\n", _k);
        }

        Merge::backTrace();
        if (h==-1)
            _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".multi2D_Merge", _boundaries);
        return _boundaries;
    }

    void Merge::backTrace() {
        printf("#data points: %d \n", N);
        int *boundaries = new int [_k];
        boundaries[_k - 1] = N - 1;
        for (int i=_k-2; i>-1; i--) {
            int t1 = boundaries[i + 1];
            boundaries[i] = _minIndexArray[t1][i + 1];
        }

        for (int i=0; i<_k; i++) {
            if (i==0)
                _boundaries.emplace_back(1, _preBoundaries[boundaries[i]].second);
            else
                _boundaries.emplace_back(_preBoundaries[boundaries[i-1]].second+1 , _preBoundaries[boundaries[i]].second);
        }
    }

    detectorH::detectorH(SuperTAD::Data &data)
    {
        _data = &data;
        std::vector<Boundary> _boundary;
    }

    void detectorH::pipeline(std::string preResult) {
        // acquire the first layer of TAD from pre-detected file or detectorH1
        std::vector<Boundary> _preBoundaries;
        if (preResult=="") {
            SuperTAD::_HD_ -= 1;
            printf("No pre-detected TAD result input, so applying for the first-time dividing.\n");
            multi::DetectorH1 dm(*_data);
            _preBoundaries = dm.execute(-1);
        } else
            SuperTAD::Reader::parseBoundariesIn8ColsFormat(_preBoundaries, preResult);
        _boundary.insert(_boundary.end(), _preBoundaries.begin(), _preBoundaries.end()); // record the first layer of TAD
        _clusters.insert(_clusters.end(),  _preBoundaries.begin(), _preBoundaries.end());

        // merge up
        std::vector<Boundary> _preboundForMerge = _preBoundaries;
        for (int i = SuperTAD::_HU_; i > 0; i--) {
            printf("Start to merge for the time: %d\n", SuperTAD::_HU_ - i + 1);
            multi::Merge dM(*_data, _preboundForMerge);
            _preboundForMerge = dM.execute(SuperTAD::_HU_ - i + 1);
            _boundary.insert(_boundary.end(), _preboundForMerge.begin(), _preboundForMerge.end());
//            _clusters.insert(_clusters.end(), _preboundForMerge.begin(), _preboundForMerge.end());// record the first layer of TAD
        }

        // go down
        std::vector<Boundary> _preboundForDivi = _preBoundaries;
        int start, end;
        double **_subMatrix;
        std::vector<Boundary> bounTmp;
        for (int i = SuperTAD::_HD_; i > 0; i--) {
            printf("Start to divide for the time: %d\n", SuperTAD::_HD_ - i + 1);
            std::vector<Boundary> _bounDiviResult;
            for (int node=0; node<_preboundForDivi.size(); node++) {
                start = _preboundForDivi[node].first;
                end = _preboundForDivi[node].second;
                printf("start=%d, end=%d, node=%d\n", start, end, node);
                if (end-start+1 <= SuperTAD::_MinSize_)
                    continue;
                else {
                    SuperTAD::Data::parseSubMatrix(_subMatrix, *_data, start, end);
                    SuperTAD::Data subdata(_subMatrix, end - start + 1);
                    subdata.init();
                    multi::DetectorH1 dD(subdata);
                    bounTmp = dD.execute(SuperTAD::_HD_ - i + 1);
                    for (int b = 0; b<bounTmp.size(); b++) {
                        _bounDiviResult.emplace_back(bounTmp[b].first+start-1, bounTmp[b].second+start-1);
                    }
                }
            }
            _boundary.insert(_boundary.end(), _bounDiviResult.begin(), _bounDiviResult.end());   // record the layers during dividing
            _preboundForDivi = _bounDiviResult;
        }

        if (!SuperTAD::_NO_OUTPUT_)
            _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".multi2D_All", _boundary);
    }

}
