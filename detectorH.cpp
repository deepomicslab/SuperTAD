//
// Created by mengbowang on 4/29/2020.
//

#include "detectorH.h"

namespace multi {

    DetectorH1::DetectorH1(Data &data) {
        _data = &data;
        _edgeCount = &data.edgeCount();
        int k=1;
        for (int i=0; i < _K_; i++) {
            _kToIdx.emplace(k++, i);
        }
        _table = new double *[_N_];
        _minIndexArray = new int *[_N_];
        _leftKArray = new int *[_N_];
        for (int i=0; i < _N_; i++) {
            _table[i] = new double [_N_]{};
            _minIndexArray[i] = new int [_N_]{};
            _leftKArray[i] = new int [_N_]{};
        }
    }


    DetectorH1::~DetectorH1() {
        for (int i=0; i < _N_; i++) {
            delete _table[i];
            delete _minIndexArray[i];
            delete _leftKArray[i];
        }
        delete _table;
        delete _minIndexArray;
        delete _leftKArray;
    }


    void DetectorH1::execute() {
        _table[0][0] = _data->getSE(0, 0, 2. *_data->_edgeSum);

        double currentVol, parentVol, binSum;

        if (_VERBOSE_)
            printf("start k=0\n");

        for (int i=1; i < _N_; i++) {
            parentVol = 2. * _data->_edgeSum;
            currentVol = _data->getVol(0, i);
            binSum = _data->getGtimesLogG(currentVol) - _data->_sumOfGtimesLogG[i];
            _table[i][0] = _data->getSE(0, i, parentVol, currentVol) + binSum/(2. *_data->_edgeSum);
        }

        if (_VERBOSE_)
            printf("finish k=0, se=%f\n", _table[_N_ - 1][0]);

        if (_VERBOSE_)
            printf("start calculating h=1\n");

        double minSE, tmpSE;
        int minIdx;
        for (int a=1; a < _K_; a++) {
            if (_VERBOSE_)
                printf("start k=%d\n", a);
            for (int b=0; b < _N_; b++) {
                minSE = std::numeric_limits<double>::infinity();
                minIdx = 0;
                for (int i=0; i<b; i++) {
                    tmpSE = _table[i][a - 1];
                    if (i+1==b)
                        tmpSE += _data->getSE(b, b, 2. * _data->_edgeSum);
                    else {
                        parentVol = 2. * _data->_edgeSum;
                        currentVol = _data->getVol(i+1, b);
                        tmpSE += _data->getSE(i + 1, b, parentVol, currentVol);
                        binSum = _data->getGtimesLogG(currentVol) - (_data->_sumOfGtimesLogG[b] - _data->_sumOfGtimesLogG[i]);
                        tmpSE += binSum / (2. * _data->_edgeSum);
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
        }

        if (_VERBOSE_)
            printf("finish calculating h=1\n");

        if (_DETERMINE_K_) {
            if (_VERBOSE_)
                std::cout << "start to determine k\n";
            else
                printf("determine k\n");

            double *y = _table[_N_ - 1];
            int tmpIdx = 0;
            double tmpValue = y[0];
            for (int i=1; i < _K_; i++) {
                if (y[i] < tmpValue) {
                    tmpValue = y[i];
                    tmpIdx = i;
                }
            }
            _k = tmpIdx + 1;
            if (_VERBOSE_)
                printf("finish determine k\n");
            printf("optimal k=%d \n", _k);
        }

        backTrace();

        _writer.writeBoundaries(_OUTPUT_ + ".h1.txt", _boundaries);
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
//        for (int i = 0; i < _boundaries.size(); i++) {
//            std::cout << "boundary[" << i << "]=(" << _boundaries[i].first << ", " << _boundaries[i].second << ")\n";
//        }
    }

}
