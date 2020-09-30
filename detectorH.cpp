//
// Created by mengbowang on 4/29/2020.
//

#include "detectorH.h"

namespace multi {

    DetectorH1::DetectorH1(Data &data) {
        _data = &data;
//        _edgeCount = &data.edgeCount();
        int k=1;
        for (int i=0; i < _K_; i++) {
            _kToIdx.emplace(k++, i);
        }
        _table = new double *[_N_];
        _minIndexArray = new int *[_N_];
        for (int i=0; i < _N_; i++) {
            _table[i] = new double [_N_]{};
            _minIndexArray[i] = new int [_N_]{};
        }
    }


    DetectorH1::~DetectorH1() {
        for (int i=0; i < _N_; i++) {
            delete _table[i];
            delete _minIndexArray[i];
        }
        delete _table;
        delete _minIndexArray;
    }


    void DetectorH1::execute() {
        _table[0][0] = _data->getSE(0, 0, _data->_doubleEdgeSum);

        double currentVol, parentVol, binSum;

        if (_VERBOSE_)
            printf("start k=0\n");

        for (int i=1; i < _N_; i++) {
            parentVol = _data->_doubleEdgeSum;
            currentVol = _data->getVol(0, i);
            binSum = _data->getGtimesLogG(currentVol) - _data->_sumOfGtimesLogG[i];
            _table[i][0] = _data->getSE(0, i, parentVol, currentVol) + binSum/_data->_doubleEdgeSum;
//            printf("base case: i=%d, parentVol=%f, currentVol=%f,table=%f\n", i, parentVol, currentVol, _table[i][0]);
        }

        if (_VERBOSE_)
            printf("finish k=0, se=%f\n", _table[_N_ - 1][0]);

        if (_VERBOSE_)
            printf("start calculating h=1\n");

        double minSE, tmpSE;
        int minIdx;
        double sumOfLeavesTmp = std::numeric_limits<double>::infinity();
        for (int a=1; a < _K_; a++) {
            if (_VERBOSE_)
                printf("start k=%d\n", a);
            for (int b=a; b < _N_; b++) {
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
            if (_VERBOSE_)
                printf("finish k=%d, structure entropy=%f\n", a, minSE);
            if (_DETERMINE_K_) {
                if (_table[_N_-1][a] < sumOfLeavesTmp){
                    _optimalK_ = a;
                    sumOfLeavesTmp = minSE;
                }
                else
                    break;
            }
        }

        if (_VERBOSE_)
            printf("finish calculating h=1\n");

        if (_DETERMINE_K_) {
            if (_optimalK_ < _K_)
                _k = _optimalK_ + 1;
            printf("determine k to be %d\n", _k);
        }

        DetectorH1::backTrace();

        _writer.writeBoundaries(_OUTPUT_ + ".multi2D.txt", _boundaries);
    }


    /*
        Back tracing with max_index_array to find the optimal route.
        param k: num of clustering
        return the positions in the optimal route.
     */
    void DetectorH1::backTrace() {
        printf("#data points: %d \n", _N_);
        int *boundaries = new int [_k];
        boundaries[_k - 1] = _N_ - 1;
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

    Merge::Merge(Data &data, std::vector<Boundary> &_preBoundaries) {
        _data = &data;
        _N_ = static_cast<int>(_preBoundaries.size());
        _table = new double *[_N_];
        _minIndexArray = new int *[_N_];
        for (int i=0; i < _N_; i++) {
            _table[i] = new double [_K_]{};
            _minIndexArray[i] = new int [_K_]{};
        }
    }


    Merge::~Merge() {
        for (int i=0; i < _N_; i++) {
            delete _table[i];
            delete _minIndexArray[i];
        }
        delete _table;
        delete _minIndexArray;
    }

    void Merge::execute() {

        double currentVol, binSum;
        int start, end;
        // record se of each Pre-node
        for (std::vector<Boundary>::iterator it=_preBoundaries.begin(); it != _preBoundaries.end(); it++) {
            start = it->first-1;
            end = it->second-1;
            currentVol = _data->getVol(start, end);
            binSum = _data->getGtimesLogG(currentVol) - (_data->_sumOfGtimesLogG[end] - _data->_sumOfGtimesLogG[start-1]);
            _prenodeSE.emplace_back(binSum/_data->_doubleEdgeSum);
        }
        if (_VERBOSE_) printf("Start calculating base case.\n");
        // base case
        int nodeStart, nodeEnd;
        double nodeVol;
        for (int i=0; i < _k; i++) {
            end = _preBoundaries[i].second - 1;
            currentVol = _data->getVol(0, end);
            _table[i][0] = 0;
            for (int j=0; j < i+1; j++) {
                _table[i][0] += _prenodeSE[j];
                nodeStart = _preBoundaries[j].first - 1;
                nodeEnd = _preBoundaries[j].second - 1;
                nodeVol = _data->getVol(nodeStart, nodeEnd);
                _table[i][0] += _data->getSE(nodeStart, nodeEnd, currentVol, nodeVol);
            }
            _table[i][0] += _data->getSE(0, end, _data->_doubleEdgeSum, currentVol);
        }
        // upper case
        if (_VERBOSE_) printf("finish k = 0, se=%f\nStart caluclating upper case.\n", _table[_k-1][0]);
        double minSE, tmpSE;
        int minIdx;
        double sumOfLeavesTmp = std::numeric_limits<double>::infinity();
        for (int a=1; a < _k; a++) {
            for (int b=a; b < _k; b++) {
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
                _minIndexArray[b][a] = minIdx;
                _table[b][a] = minSE;
            }
            if (_VERBOSE_)
                printf("finish k=%d, structure entropy=%f\n", a, minSE);
            if (_DETERMINE_K_) {
                if (_table[_N_-1][a] < sumOfLeavesTmp){
                    _optimalK_ = a;
                    sumOfLeavesTmp = minSE;
                }
                else
                    break;
            }
        }

        if (_VERBOSE_)
            printf("finish calculating h=1\n");

        if (_DETERMINE_K_) {
            if (_optimalK_ < _K_)
                _k = _optimalK_ + 1;
            printf("determine k to be %d\n", _k);
        }

        Merge::backTrace();

        _writer.writeBoundaries(_OUTPUT_ + ".multi2D_Merge.txt", _boundaries);
    }

    void Merge::backTrace() {
        printf("#data points: %d \n", _N_);
        int *boundaries = new int [_k];
        boundaries[_k - 1] = _N_ - 1;
        for (int i=_k-2; i>-1; i--) {
            int t1 = boundaries[i + 1];
            boundaries[i] = _minIndexArray[t1][i + 1];
        }

        for (int i=0; i<_k; i++) {
            if (i==0)
                _boundaries.emplace_back(1, _preBoundaries[boundaries[i]].second);
            else
                _boundaries.emplace_back(_preBoundaries[boundaries[i-1]].first , _preBoundaries[boundaries[i]].second);
        }
        delete &boundaries;
    }

    detectorH::detectorH(Data &data)
    {
        _data = &data;
        std::vector<Boundary> _boundary;
    }

    void detectorH::pipeline(std::string preResult) {
        // acquire the first layer of TAD from pre-detected file or detectorH1
        std::vector<Boundary> _preBoundaries;
        if (preResult=="") {
            printf("No pre-detected TAD result input.\n");
            multi::DetectorH1 dm(*_data);
            dm.execute();
            preResult = _OUTPUT_ + ".multi2D.txt";
        }
        if (_VERBOSE_) printf("Obtain the preResult as %s \n", preResult.c_str());
        Reader::parseBoundariesIn8ColsFormat(_preBoundaries, preResult);
        _boundary.insert(_boundary.end(), _preBoundaries.begin(), _preBoundaries.end()); // record the first layer of TAD

        // merge up
        std::vector<Boundary> _preboundForMerge = _preBoundaries;
        for (int i = _HU_; i>0; i--) {
            multi::Merge dM(*_data, _preboundForMerge);
            dM.execute();
            preResult = _OUTPUT_ + ".multi2D_Merge.txt";
            Reader::parseBoundariesIn8ColsFormat(_preboundForMerge, preResult);
            _boundary.insert(_boundary.end(), _preboundForMerge.begin(), _preboundForMerge.end());   // record the layers during merging
        }

        // go down
        std::vector<Boundary> _preboundForDivi = _preBoundaries;
        int start, end;
        double **_subMatrix;
        std::vector<Boundary> bounTmp;
        for (int i = _HD_-1; i>0; i--) {
            std::vector<Boundary> _bounDiviResult;
            for (int node=0; node<_preboundForDivi.size(); node++) {
                start = _preboundForDivi[node].first;
                end = _preboundForDivi[node].second;

                Data::parsesubMatrix(_subMatrix, *_data, start, end);
                Data subdata(_subMatrix);
                subdata.init();
                multi::DetectorH1 dD(subdata);
                dD.execute();
                Reader::parseBoundariesIn8ColsFormat(bounTmp, _OUTPUT_ + ".multi2D.txt");
                for (int b = 0; b<bounTmp.size(); b++) {
                    _bounDiviResult.emplace_back(bounTmp[b].first+start-1, bounTmp[b].second+start-1);
                }
            }
            _boundary.insert(_boundary.end(), _bounDiviResult.begin(), _bounDiviResult.end());   // record the layers during dividing
            _preboundForDivi = _bounDiviResult;
        }

        _writer.writeBoundaries(_OUTPUT_ + ".multi2D_All.txt", _boundary);
    }

}
