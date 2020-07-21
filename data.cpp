//
// Created by wang mengbo on 2019-09-01.
//

#include "data.h"


Data::Data(std::string fileName)
{
//    Reader::parseMatrix(_contactMat, _INPUT_);
    Reader::parseMatrix2Table(_contactArray, _INPUT_);
    printf("number of bins is %d\n", _N_);

    _logVolTable = new double *[_N_];
    _volTable = new double *[_N_];
    _edgeCountArray = new double *[_N_];
    for (int s=0; s<_N_; s++) {
        _logVolTable[s] = new double [_N_-s];
        _volTable[s] = new double [_N_-s];
        _edgeCountArray[s] = new double [_N_];
    }

    if (_BINARY_ && _FAST_) {
        if (_PENALTY_<0) {
            _PENALTY_ = ceil(_N_/10);
        }
        printf("set fast mode penalty to %d\n", _PENALTY_);
    }

    if (_DETERMINE_K_ && _K_ < 0) {
        _K_ = sqrt(_N_) + 5;
        printf("set max K to %d\n", _K_);
    }
}


Data::~Data()
{
    for (int s=0; s<_N_; s++) {
        delete _logVolTable[s];
        delete _volTable[s];
    }
    delete _logVolTable;
    delete _volTable;
}


void Data::init()
{
    std::clock_t tTmp;

    if (_VERBOSE_) {
        std::cout << "start initialization\n";
        tTmp = std::clock();
    }
    else
        printf("initialize\n");

//    // calculate edge count(sum)
//    _edgeCountMat.resize(_N_, _N_);
//    // intra
//    for (int k=1; k < _N_; k++) {
//        for (int i=0; i < _N_ - k; i++) {
//            int j = i+k;
//            double intra = _contactMat.coeff(i, j);
//            if (j-1 > 0)
//                intra += _edgeCountMat.coeff(i, j - 1);
//            if (i+1 < _N_)
//                intra += _edgeCountMat.coeff(i + 1, j);
//            if (j-1>0 && i+1 < _N_)
//                intra -= _edgeCountMat.coeff(i + 1, j - 1);
//            if (abs(intra) < _THRESHOLD_ || intra < 0)
//                _edgeCountMat(i, j) = 0;
//            else
//                _edgeCountMat(i, j) = intra;
//        }
//    }
//    // inter
//    for (int i=0; i < _N_; i++) {
//        for (int j=i; j < _N_; j++) {
//            double inter = _edgeCountMat.coeff(0, j) + _edgeCountMat.coeff(i, _N_ - 1) - 2 * _edgeCountMat.coeff(i, j);
//            if (i-1 > 0)
//                inter -= _edgeCountMat.coeff(0, i - 1);
//            if (j+1 < _N_)
//                inter -= _edgeCountMat.coeff(j + 1, _N_ - 1);
//            if (abs(inter) < _THRESHOLD_ || inter < 0) {
//                _edgeCountMat(j, i) = 0;
//            } else
//                _edgeCountMat(j, i) = inter;
//        }
//    }
//    setEdgeSum();

    // edge sum
    int i, j, k;
    for (k=1; k<_N_; k++) {
        for (i=0; i<_N_-k; i++) {
            j = i+k;
            double intra = _contactArray[i][j];
            if (j-1 > 0)
                intra += _edgeCountArray[i][j-1];
            if (i+1 < _N_)
                intra += _edgeCountArray[i+1][j];
            if (j-1>0 && i+1 < _N_)
                intra -= _edgeCountArray[i+1][j-1];
            if (abs(intra) < _THRESHOLD_ || intra < 0)
                _edgeCountArray[i][j] = 0;
            else
                _edgeCountArray[i][j] = intra;
        }
    }
    for (int i=0; i < _N_; i++) {
        for (int j=i; j < _N_; j++) {
            double inter = _edgeCountArray[0][j] + _edgeCountArray[i][_N_-1] - 2 * _edgeCountArray[i][j];
            if (i-1 > 0)
                inter -= _edgeCountArray[0][i-1];
            if (j+1 < _N_)
                inter -= _edgeCountArray[j+1][_N_-1];
            if (abs(inter) < _THRESHOLD_ || inter < 0) {
                _edgeCountArray[j][i] = 0;
            } else
                _edgeCountArray[j][i] = inter;
        }
    }
    _edgeSum = _edgeCountArray[0][_N_-1];
    _doubleEdgeSum = 2. * _edgeSum;

    // calculate volTble and logTable
    for (int s=0; s<_N_; s++) {
        for (int e=s; e<_N_; e++) {
            _volTable[s][e-s] = getVol(s, e);
            if (_volTable[s][e-s]>0)
                _logVolTable[s][e-s] = log2(_volTable[s][e-s]);
            else
                _logVolTable[s][e-s] = -1;
        }
    }

    // calculate sum of g*log(g)
//    _sumOfGtimesLogG.emplace_back(getGtimesLogG(_edgeCountMat.coeff(0, 0)) );
    _sumOfGtimesLogG.emplace_back(getGtimesLogG(_edgeCountArray[0][0]));
    for (int i=1; i < _N_; i++) {
//        _sumOfGtimesLogG.emplace_back(_sumOfGtimesLogG[i-1] + getGtimesLogG(_edgeCountMat.coeff(i, i)));
        _sumOfGtimesLogG.emplace_back(_sumOfGtimesLogG[i-1] + getGtimesLogG(_edgeCountArray[i][i]));
    }

    if (_VERBOSE_)
        printf("finish initialization\n");

    if (_DEBUG_)
        printf("initialization consumes %fs\n", (float)(std::clock()-tTmp)/CLOCKS_PER_SEC);
}


//void Data::setEdgeSum() {
//    _edgeSum = _edgeCountMat.coeff(0, _N_ - 1);
//    _doubleEdgeSum = 2. * _edgeSum;
//}


double Data::getVol(int s, int e)
{
    if (e<s)
        fprintf(stderr, "s=%d, e=%d, e<s\n", s, e);

//    if (s!=e)
//        return 2. * _edgeCountMat.coeff(s, e) + _edgeCountMat.coeff(e, s);
//    else
//        return _edgeCountMat.coeff(e, s);
    if (s!=e)
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
    if(currentVol > 0 && parentVol >= currentVol)
        return _edgeCountArray[e][s] / _doubleEdgeSum * log2(parentVol / currentVol);
    return 0;
}


double Data::getSEwithLogPV(int s, int e, double logPV)
{
//    if (_LOG_VOL_TABLE_)
//        return _logVolTable[s][e-s] > 0 ? _edgeCountMat(e, s) / _doubleEdgeSum * (logPV - _logVolTable[s][e-s]) : 0;
    if (_PRE_LOG_)
        return _logVolTable[s][e-s] > 0 ? _edgeCountArray[e][s] / _doubleEdgeSum * (logPV - _logVolTable[s][e-s]) : 0;
    else {
        double currentVol = getVol(s, e);
        if(currentVol > 0) {
            double logCV = log2(currentVol);
//            if (logPV >= logCV)
//                return _edgeCountMat(e, s) / _doubleEdgeSum * (logPV - logCV);
            if (logPV >= logCV)
                return _edgeCountArray[e][s] / _doubleEdgeSum * (logPV - logCV);
        }
        return 0;
    }
}


double Data::getSE(int s, int e, double parentVol, double currentVol) {
    // g / edge_sum * log2(V_p / V)
//    if(currentVol > 0 && parentVol >= currentVol)
//        return _edgeCountMat(e, s) / _doubleEdgeSum * log2(parentVol / currentVol);
    if(currentVol > 0 && parentVol >= currentVol)
        return _edgeCountArray[e][s] / _doubleEdgeSum * log2(parentVol / currentVol);
    return 0;
}


double Data::getSEwithLogs(int s, int e, double logPV, double logCV)
{
//    return logPV>=logCV ? _edgeCountMat(e, s) / _doubleEdgeSum * (logPV - logCV) : 0;
    return logPV>=logCV ? _edgeCountArray[e][s] / _doubleEdgeSum * (logPV - logCV) : 0;
}


double Data::getSEwithLogDiff(int s, int e, double logDiff)
{
    // g / edge_sum * log2(V_p / V)
//    return logDiff > 0 ? _edgeCountMat(e, s) / _doubleEdgeSum * logDiff : 0;
    return logDiff > 0 ? _edgeCountArray[e][s] / _doubleEdgeSum * logDiff : 0;
}


double Data::getGtimesLogG(double binG) {
    if (binG==0)
        return 0;
    else
        return binG * log2(binG);
}
