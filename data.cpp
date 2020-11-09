//
// Created by wang mengbo on 2019-09-01.
//

#include "data.h"


namespace SuperTAD
{

    Data::Data(std::string input)
    {
        // Reader::parseMatrix(_contactMat, _INPUT_);
        Reader::parseMatrix2Table(_contactArray, input);
        printf("number of bins is %d\n", _N_);

        if (_K_ <= 0)
            _K_ = _N_ / _MinSize_;

        _logVolTable = new double *[_N_];
        _volTable = new double *[_N_];
        _edgeCountArray = new double *[_N_];
        for (int s = 0; s < _N_; s++) {
            _logVolTable[s] = new double[_N_ - s]{};
            _volTable[s] = new double[_N_ - s]{};
            _edgeCountArray[s] = new double[_N_]{};
        }
        if (_BINARY_ && _FAST_) {
            if (_PENALTY_ < 0) {
                _PENALTY_ = ceil(_N_ / 10);
            }
            printf("set fast mode penalty to %d\n", _PENALTY_);
        }

        if (_DETERMINE_K_ && _K_ < 0) {
            _K_ = _N_ / _MinSize_;
            printf("set max K to %d\n", _K_);
        }
    }


    SuperTAD::Data::Data(double **array, int n)
    {
        _initByPointer = true;
        _N_ = n;
        printf("initing the data class via 2d array\n");
        printf("number of bins is %d\n", _N_);
        _contactArray = array;
        // printf("countactArray=\n");
        // utils::print2Darray(_contactArray, _N_, _N_);

        if (_N_ < _K_) {
            _K_ = _N_ / _MinSize_;
            // printf("reset max K to %d\n", _N_);
        }

        _logVolTable = new double *[_N_];
        _volTable = new double *[_N_];
        _edgeCountArray = new double *[_N_];
        for (int s = 0; s < _N_; s++) {
            _logVolTable[s] = new double[_N_ - s]{};
            _volTable[s] = new double[_N_ - s]{};
            _edgeCountArray[s] = new double[_N_]{};
        }

        if (_BINARY_ && _FAST_) {
            if (_PENALTY_ < 0) {
                _PENALTY_ = ceil(_N_ / 10);
            }
            printf("set fast mode penalty to %d\n", _PENALTY_);
        }

        if (_DETERMINE_K_ && _K_ < 0) {
            _K_ = _N_ / _MinSize_;
            printf("set max K to %d\n", _K_);
        }
    }



    Data::~Data()
    {
        for (int s=0; s<_N_; s++) {
            delete [] _logVolTable[s];
            delete [] _volTable[s];
            delete [] _edgeCountArray[s];
            if (!_initByPointer)
                delete [] _contactArray[s];
        }
        delete [] _logVolTable;
        delete [] _volTable;
        delete [] _edgeCountArray;
        if (!_initByPointer)
            delete [] _contactArray;
    }


    void Data::init()
    {
        std::clock_t tTmp;

        if (_VERBOSE_) {
            std::cout << "start initialization\n";
            tTmp = std::clock();
        }

//    printf("_N_=%d\n", _N_);
//    printf("before init: contactArray=\n");
//    utils::print2Darray(_contactArray, _N_, _N_);
//    printf("before init: edgeCountArray=\n");
//    utils::print2Darray(_edgeCountArray, _N_, _N_);

        // edge sum
        int i, j, k;
        for (k = 1; k < _N_; k++) {
            for (i = 0; i < _N_ - k; i++) {
                j = i + k;
                double intra = _contactArray[i][j];
//            printf("intra=contactArray[%d][%d]=%f\n", i, j, intra);
                if (j - 1 > 0) {
                    intra += _edgeCountArray[i][j - 1];
//                printf("intra+edgeCountArray[%d][%d]=%f\n", i, j-1, intra);
                }
                if (i + 1 < _N_) {
                    intra += _edgeCountArray[i + 1][j];
//                printf("intra+edgeCountArray[%d][%d]=%f\n", i+1, j, intra);
                }
                if (j - 1 > 0 && i + 1 < _N_) {
                    intra -= _edgeCountArray[i + 1][j - 1];
//                printf("intra-dgeCountArray[%d][%d]=%f\n", i+1, j-1, intra);
                }
//            printf("intra=%f, abs(intra)=%f, _THRESHOLD_=%f\n", intra, std::abs(intra), _THRESHOLD_);

                if (std::abs(intra) < _THRESHOLD_ || intra < 0) {
//                fprintf(stderr, "[ERROR] intra=%f, abs(intra)=%f, _THRESHOLD_=%f\n", intra, std::abs(intra), _THRESHOLD_);
//                if (std::abs(intra) < _THRESHOLD_)
//                    printf("intra[%d][%d]<THRESHOLD(%f)\n", i, j, _THRESHOLD_);
//                else if (intra < 0) {
//                    printf("intra[%d][%d]<0\n", i, j);
//                }
                    _edgeCountArray[i][j] = 0;
                } else
                    _edgeCountArray[i][j] = intra;
//            printf("edgeCountArray[%d][%d]=%f\n", i, j, _edgeCountArray[i][j]);
            }
        }
        for (int i = 0; i < _N_; i++) {
            for (int j = i; j < _N_; j++) {
                double inter = _edgeCountArray[0][j] + _edgeCountArray[i][_N_ - 1] - 2 * _edgeCountArray[i][j];
                if (i - 1 > 0)
                    inter -= _edgeCountArray[0][i - 1];
                if (j + 1 < _N_)
                    inter -= _edgeCountArray[j + 1][_N_ - 1];

                if (std::abs(inter) < _THRESHOLD_ || inter < 0) {
//                if (std::abs(inter) < _THRESHOLD_)
//                    printf("inter[%d][%d]<THRESHOLD(%f)\n", i, j, _THRESHOLD_);
//                else if (inter < 0) {
//                    printf("inter[%d][%d]<0\n", i, j);
//                }
                    _edgeCountArray[j][i] = 0;
                } else
                    _edgeCountArray[j][i] = inter;
//            printf("edgeCountArray[%d][%d]=%f\n", i, j, _edgeCountArray[i][j]);
            }
        }
//    printf("edgeCountArray=\n");
//    utils::print2Darray(_edgeCountArray, _N_, _N_);

        _edgeSum = _edgeCountArray[0][_N_ - 1];
        _doubleEdgeSum = 2. * _edgeSum;
//    printf("edgesum=%f, doubleEdgesum=%f\n", _edgeSum, _doubleEdgeSum);

        // calculate volTble and logTable
        for (int s = 0; s < _N_; s++) {
            for (int e = s; e < _N_; e++) {
                _volTable[s][e - s] = getVol(s, e);
                if (_volTable[s][e - s] > 0)
                    _logVolTable[s][e - s] = log2(_volTable[s][e - s]);
                else
                    _logVolTable[s][e - s] = -1;
            }
        }

        // calculate sum of g*log(g)
//    _sumOfGtimesLogG.emplace_back(getGtimesLogG(_edgeCountMat.coeff(0, 0)) );
        _sumOfGtimesLogG.emplace_back(getGtimesLogG(_edgeCountArray[0][0]));
        for (int i = 1; i < _N_; i++) {
//        _sumOfGtimesLogG.emplace_back(_sumOfGtimesLogG[i-1] + getGtimesLogG(_edgeCountMat.coeff(i, i)));
            _sumOfGtimesLogG.emplace_back(_sumOfGtimesLogG[i - 1] + getGtimesLogG(_edgeCountArray[i][i]));
        }

        if (_VERBOSE_)
            printf("finish initialization\n");

        if (_DEBUG_)
            printf("initialization consumes %fs\n", (float) (std::clock() - tTmp) / CLOCKS_PER_SEC);
    }


    void Data::parseSubMatrix(double **&subMatrix, Data &Matrix, int start, int end)
    {
        int N = end - start + 1;
        subMatrix = new double *[N];
        for (int i = 0; i < N; i++) {
            subMatrix[i] = new double[N]{};
            for (int j = 0; j < N; j++) {
                subMatrix[i][j] = Matrix._contactArray[i + start - 1][j + start - 1];
            }
        }
    }


    double Data::getG(int s, int e)
    {
        return _edgeCountArray[e][s];
    }


    double Data::getVol(int s, int e)
    {
        if (e < s)
            fprintf(stderr, "s=%d, e=%d, e<s\n", s, e);

//    if (s!=e)
//        return 2. * _edgeCountMat.coeff(s, e) + _edgeCountMat.coeff(e, s);
//    else
//        return _edgeCountMat.coeff(e, s);
        if (s != e)
            return 2. * _edgeCountArray[s][e] + _edgeCountArray[e][s];
        else
            return _edgeCountArray[e][s];
    }


    double Data::getSE(int s, int e, double parentVol)
    {
        // g / edge_sum * log2(V_p / V)
        double currentVol = getVol(s, e);
//    if(currentVol > 0 && parentVol >= currentVol)
//        return _edgeCountMat(e, s) / _doubleEdgeSum * log2(parentVol / currentVol);
        if (currentVol > 0 && parentVol >= currentVol)
            return _edgeCountArray[e][s] / _doubleEdgeSum * log2(parentVol / currentVol);
        return 0;
    }


    double Data::getSEwithLogPV(int s, int e, double logPV)
    {
//    if (_LOG_VOL_TABLE_)
//        return _logVolTable[s][e-s] > 0 ? _edgeCountMat(e, s) / _doubleEdgeSum * (logPV - _logVolTable[s][e-s]) : 0;
        if (_PRE_LOG_)
            return _logVolTable[s][e - s] > 0 ? _edgeCountArray[e][s] / _doubleEdgeSum *
                                                (logPV - _logVolTable[s][e - s]) : 0;
        else {
            double currentVol = getVol(s, e);
            if (currentVol > 0) {
                double logCV = log2(currentVol);
//            if (logPV >= logCV)
//                return _edgeCountMat(e, s) / _doubleEdgeSum * (logPV - logCV);
                if (logPV >= logCV)
                    return _edgeCountArray[e][s] / _doubleEdgeSum * (logPV - logCV);
            }
            return 0;
        }
    }


    double Data::getSE(int s, int e, double parentVol, double currentVol)
    {
        // g / edge_sum * log2(V_p / V)
        if (currentVol > 0 && parentVol >= currentVol)
            return _edgeCountArray[e][s] / _doubleEdgeSum * log2(parentVol / currentVol);
        return 0;
    }


    double Data::getSEwithLogs(int s, int e, double logPV, double logCV)
    {
//    return logPV>=logCV ? _edgeCountMat(e, s) / _doubleEdgeSum * (logPV - logCV) : 0;
        return logPV >= logCV ? _edgeCountArray[e][s] / _doubleEdgeSum * (logPV - logCV) : 0;
    }


    double Data::getSEwithLogDiff(int s, int e, double logDiff)
    {
        // g / edge_sum * log2(V_p / V)
//    return logDiff > 0 ? _edgeCountMat(e, s) / _doubleEdgeSum * logDiff : 0;
        return logDiff > 0 ? _edgeCountArray[e][s] / _doubleEdgeSum * logDiff : 0;
    }


    double Data::getGtimesLogG(double binG)
    {
        if (binG == 0)
            return 0;
        else
            return binG * log2(binG);
    }

}
