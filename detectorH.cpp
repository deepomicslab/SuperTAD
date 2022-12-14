//
// Created by mengbowang on 4/29/2020.
//

#include "detectorH.h"


namespace SuperTAD::multi {

    Partition::Partition(SuperTAD::Data &data) {
        _data = &data;

        if (SuperTAD::_DETERMINE_K_)
        {
            SuperTAD::_K_ = SuperTAD::_N_;
            printf("Set max K to be %d\n", SuperTAD::_K_);
        } else
        {
            printf("The optimal K is %d\n", SuperTAD::_K_);
            _k = SuperTAD::_K_;
        }
        _table = new double *[SuperTAD::_N_];
        _minIndexArray = new int *[SuperTAD::_N_];
        for (int i=0; i < SuperTAD::_N_; i++) {
            _table[i] = new double [SuperTAD::_N_]{};
            _minIndexArray[i] = new int [SuperTAD::_N_]{};
        }
    }


    Partition::~Partition() {
        for (int i=0; i < SuperTAD::_N_; i++) {
            delete _table[i];
            delete _minIndexArray[i];
        }
        delete _table;
        delete _minIndexArray;
    }


    std::vector<Boundary> Partition::execute(int h) {
        _table[0][0] = _data->getSE(0, 0, 0, SuperTAD::_N_-1);

        double binSum;

        for (int i=1; i < SuperTAD::_N_; i++) {
            _table[i][0] = _data->getSE(0, i, 0, SuperTAD::_N_-1);
            binSum = _data->getVol(0, i) * _data->_logVolTable[0][i] - _data->_sumOfGtimesLogG[i];
            _table[i][0] += binSum / _data->_doubleEdgeSum;
//            printf("base case: i=%d, parentVol=%f, currentVol=%f,table=%f\n", i, parentVol, currentVol, _table[i][0]);
        }

        if (SuperTAD::_VERBOSE_)
            printf("finish k=0, se=%f\n", _table[SuperTAD::_N_ - 1][0]);

        double minSE, tmpSE;
        int minIdx;
        double sumOfLeavesTmp = _table[SuperTAD::_N_ - 1][0];

        for (int a=1; a < SuperTAD::_K_; a++) { // with set K or the max_K

            for (int b=a; b < SuperTAD::_N_; b++) {
                minSE = std::numeric_limits<double>::infinity();
                minIdx = 0;
                for (int i=a-1; i<b; i++) {
                    tmpSE = _table[i][a - 1];
                    if (i+1==b)
                        tmpSE += _data->getSE(b, b, 0, SuperTAD::_N_-1);
                    else {
                        tmpSE += _data->getSE(i+1, b, 0, SuperTAD::_N_-1);
                        binSum = _data->getVol(i+1, b) * _data->_logVolTable[i+1][b-i-1] - (_data->_sumOfGtimesLogG[b] - _data->_sumOfGtimesLogG[i]);
                        tmpSE += binSum / _data->_doubleEdgeSum;
                    }
                    if (tmpSE - minSE < - SuperTAD::_THRESHOLD_) {
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
                if (_table[SuperTAD::_N_ - 1][a] - sumOfLeavesTmp < - SuperTAD::_THRESHOLD_){
                    SuperTAD::_OPTIMAL_K_ = a;
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
            if (SuperTAD::_OPTIMAL_K_ < SuperTAD::_K_)
                _k = SuperTAD::_OPTIMAL_K_ + 1;
            printf("determine k to be %d\n", _k);
        }

        Partition::backTrace();

        if (h==-1)
            if (_BIN_LIST_)
                _writer.writeBoundaries(SuperTAD::_OUTPUT_ + ".multi2D.txt", _boundaries);
            else
                _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".multi2D", _boundaries);
        return _boundaries;
    }


    void Partition::backTrace() {
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
        N = (int) _preBoundaries.size();    // the number of boundaries
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
        double binSum;
        int start, end;
        // record se of each Pre-node
        for (std::vector<Boundary>::iterator it = _preBoundaries.begin(); it != _preBoundaries.end(); ++it) {
            start = it->first-1;
            end = it->second-1;
            binSum = _data->getVol(start, end) * _data->_logVolTable[start][end-start];
            if (start == 0)
                binSum -= _data->_sumOfGtimesLogG[end];
            else
                binSum -= (_data->_sumOfGtimesLogG[end] - _data->_sumOfGtimesLogG[start-1]);
            _prenodeSE.emplace_back(binSum / _data->_doubleEdgeSum);
//            printf("start=%d, end=%d, prenodeSE=%f\n", start, end, binSum/_data->_doubleEdgeSum);
        }
        if (SuperTAD::_VERBOSE_) printf("Start calculating base case. #node=%d\n", N);
        // base case
        int nodeStart, nodeEnd;

        for (int i=0; i < N; i++) {
            end = _preBoundaries[i].second - 1;
            for (int j=0; j < i+1; j++) {
                _table[i][0] += _prenodeSE[j];
                nodeStart = _preBoundaries[j].first - 1;
                nodeEnd = _preBoundaries[j].second - 1;
                _table[i][0] += _data->getSE(nodeStart, nodeEnd, 0, end);
//                printf("nodestart=%d, nodeend=%d, nodevol=%f, _table=%f\n", nodeStart, nodeEnd, nodeVol, _table[i][0]);
            }
            _table[i][0] += _data->getSE(0, end, 0, SuperTAD::_N_ - 1);
        }

        // upper case
        if (SuperTAD::_VERBOSE_) printf("finish k = 0, se=%f\nStart caluclating upper case.\n", _table[N - 1][0]);
        double minSE, tmpSE;
        int minIdx;
        double sumOfLeavesTmp = _table[N-1][0];
        for (int a=1; a < N; a++) {
            for (int b=a; b < N; b++) {
                minSE = std::numeric_limits<double>::infinity();
                for (int i=a-1; i<b; i++) {
                    tmpSE = _table[i][a - 1];
                    if (i+1==b) {
                        tmpSE += _prenodeSE[b];
                        tmpSE += _data->getSE(_preBoundaries[b].first - 1, _preBoundaries[b].second - 1,0, SuperTAD::_N_-1);
                    }
                    else {
                        start = _preBoundaries[i+1].first-1;
                        end = _preBoundaries[b].second-1;
                        for (int node=i+1; node<b+1; node++) {
                            nodeStart = _preBoundaries[node].first-1;
                            nodeEnd = _preBoundaries[node].second-1;
                            tmpSE += _prenodeSE[node];
                            tmpSE += _data->getSE(nodeStart, nodeEnd, start, end);
                        }
                        tmpSE += _data->getSE(start, end, 0, SuperTAD::_N_-1);
                    }
                    if (tmpSE - minSE < -SuperTAD::_THRESHOLD_) {
                        minSE = tmpSE;
                        minIdx = i;
                    }
                }
                if (b == N-1)
                    printf("a=%d, b=%d, minIdx=%d, minSE=%f\n", a, b, minIdx, minSE);
                _minIndexArray[b][a] = minIdx;
                _table[b][a] = minSE;
            }
            if (SuperTAD::_VERBOSE_)
                printf("finish k = %d, structure entropy=%f, mid = %d\n", a, minSE, minIdx);
            if (SuperTAD::_DETERMINE_K_) {
                if (_table[N-1][a] - sumOfLeavesTmp < - SuperTAD::_THRESHOLD_){
                    SuperTAD::_OPTIMAL_K_ = a;
                    sumOfLeavesTmp = minSE;
                }
                else
                    break;
            }
        }

        if (SuperTAD::_VERBOSE_)
            printf("finish calculating h=1\n");

        if (SuperTAD::_DETERMINE_K_) {
            if (SuperTAD::_OPTIMAL_K_ < SuperTAD::_K_)
                _k = SuperTAD::_OPTIMAL_K_ + 1;
            printf("determine k to be %d\n", _k);
        }

        Merge::backTrace();
        if (h==-1)
            _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".multi2D_Merge", _boundaries);
        return _boundaries;
    }


    void Merge::backTrace() {
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
            if (SuperTAD::_V1_)
            {
                multi::Partition dm(*_data);
                _preBoundaries = dm.execute(-1);
            } else
            {
                multi::PartitionV2 dm(*_data);
                _preBoundaries = dm.execute(-1);
            }

        } else
            SuperTAD::Reader::parseBoundariesIn8ColsFormat(_preBoundaries, preResult);
        _boundary.insert(_boundary.end(), _preBoundaries.begin(), _preBoundaries.end()); // record the first layer of TAD
        _clusters.insert(_clusters.end(),  _preBoundaries.begin(), _preBoundaries.end());

        // merge up
        std::vector<Boundary> _preboundForMerge = _preBoundaries;
        std::vector<Boundary> _tmpbound;
        for (int i = SuperTAD::_HU_; i > 0; i--) {
            printf("Start to merge for the time: %d\n", SuperTAD::_HU_ - i + 1);
            if (SuperTAD::_V1_)
            {
                multi::Merge dM(*_data, _preboundForMerge);
                _tmpbound = dM.execute(SuperTAD::_HU_ - i + 1);
            } else
            {
                multi::MergeV2 dM(*_data, _preboundForMerge);
                _tmpbound = dM.execute(SuperTAD::_HU_ - i + 1);
            }

            if (_tmpbound != _preboundForMerge and _tmpbound.size() > 1)
            {
                _preboundForMerge = _tmpbound;
                _boundary.insert(_boundary.end(), _preboundForMerge.begin(),
                                 _preboundForMerge.end());  // record the layers during merging
            } else
                break;
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
//                printf("start=%d, end=%d, node=%d\n", start, end, node);
                if (end-start+1 <= 3) // set the min TAD size as 3
                    continue;
                else {
                    SuperTAD::Data::parseSubMatrix(_subMatrix, *_data, start, end);
                    SuperTAD::Data subdata(_subMatrix, end - start + 1);
                    subdata.init();
                    if (SuperTAD::_V1_)
                    {
                        multi::Partition dD(subdata);
                        bounTmp = dD.execute(SuperTAD::_HD_ - i + 1);
                    } else
                    {
                        multi::PartitionV2 dD(subdata);
                        bounTmp = dD.execute(SuperTAD::_HD_ - i + 1);
                    }

                    for (int b = 0; b<bounTmp.size(); b++) {
                        _bounDiviResult.emplace_back(bounTmp[b].first+start-1, bounTmp[b].second+start-1);
                    }
                }
            }
            _boundary.insert(_boundary.end(), _bounDiviResult.begin(), _bounDiviResult.end());   // record the layers during dividing
            _preboundForDivi = _bounDiviResult;
        }

        int height = SuperTAD::_HD_ + SuperTAD::_HU_ + 1;
        _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".multi2D_AllH" + std::to_string(height), _boundary);

    }

}
