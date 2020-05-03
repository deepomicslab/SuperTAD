//
// Created by mengbowang on 4/29/2020.
//

#include "detectorH.h"

namespace multi {

    DetectorH1::DetectorH1(Data &data) {
        _data = &data;
        _edgeCount = &data.edgeCount();
        int k=1;
        for (int i=0; i<_K; i++) {
            _kToIdx.emplace(k++, i);
        }
        _table = new double *[_N];
        _minIndexArray = new int *[_N];
        _leftKArray = new int *[_N];
        for (int i=0; i<_N; i++) {
            _table[i] = new double [_N]{};
            _minIndexArray[i] = new int [_N]{};
            _leftKArray[i] = new int [_N]{};
        }
    }


    DetectorH1::~DetectorH1() {
        for (int i=0; i<_N; i++) {
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
        for (int i=1; i<_N; i++) {
            double currentVol = _data->getVol(0, i);
            double binSum = _data->getGtimesLogG(currentVol) - _data->_sumOfGtimesLogG[i];
            _table[i][0] = _data->getSE(0, i, 2. * _data->_edgeSum, currentVol) + binSum/(2. *_data->_edgeSum);
        }
        printf("finish k=0, SE=%f \n", _table[_N-1][0]);
        std::cout << "start to calculate upper case\n";
        for (int a=1; a<_K; a++) {
            double minTmp;
            int minIdx;
            for (int b=0; b<_N; b++) {
                minTmp= std::numeric_limits<double>::infinity();
                minIdx = 0;
                for (int i=0; i<b; i++) {
                    double tmp = _table[i][a-1];
                    if (i+1==b)
                        tmp += _data->getSE(b, b, 2. * _data->_edgeSum);
                    else {
                        double currentVol = _data->getVol(i+1, b);
                        tmp += _data->getSE(i+1, b, 2. * _data->_edgeSum, currentVol);
                        double binSum = _data->getGtimesLogG(currentVol) - (_data->_sumOfGtimesLogG[b] - _data->_sumOfGtimesLogG[i]);
                        tmp += binSum / (2. * _data->_edgeSum);
                    }
                    if (tmp<minTmp) {
                        minTmp = tmp;
                        minIdx = i;
                    }
                }
                _minIndexArray[b][a] = minIdx;
                _table[b][a] = minTmp;
            }
            printf("finish k=%d, structure entropy=%f \n", a, minTmp);
        }
        if (_DETERMINE_K) {
            std::cout << "start to determine k\n";
            double *y = _table[_N-1];
            int tmpIdx = 0;
            double tmpValue = y[0];
            for (int i=1; i<_K; i++) {
                if (y[i] < tmpValue) {
                    tmpValue = y[i];
                    tmpIdx = i;
                }
            }
            _k = tmpIdx + 1;
            printf("the optimal k is %d \n", _k);
        }

        backTrace();

        _writer.writeBoundaries(_INPUT + ".h1.txt", _boundaries);
    }


    /*
        Back tracing with max_index_array to find the optimal route.
        param k: num of clustering
        return the positions in the optimal route.
     */
    void DetectorH1::backTrace() {
        printf("#data points: %d \n", _N);
        int *boundaries = new int [_k];
        boundaries[_k - 1] = _N - 1;
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
        for (int i = 0; i < _boundaries.size(); i++) {
            std::cout << "boundary[" << i << "]=(" << _boundaries[i].first << ", " << _boundaries[i].second << ")\n";
        }
    }

}
