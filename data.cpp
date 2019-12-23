//
// Created by wang mengbo on 2019-09-01.
//

#include "data.h"


Data::Data(std::string fileName)
{
    _reader = new Reader(fileName);
    _N = _reader->parse(_contactMat);
    std::cout << "#bins: " << _N << std::endl;
    if (_K < 0)
        _K = _N / 5;
}


Data::~Data()
{
    delete _reader;
}


void Data::init()
{
    std::clock_t t = std::clock();
    if (_VERBOSE)
        std::cout << "start data initialization\n";
    _edgeCount.resize(_N, _N);
    for (int i = 0; i < _N; i++) {
//    if (_VERBOSE)
//      std::cout << "i=" << i << std::endl;
        for (int j = i; j < _N; j++) {
            int n = j - i + 1;
//      if (_VERBOSE)
//        std::cout << "i=" << i << ", j=" << j << ", n=" << n << std::endl;
            Eigen::MatrixXd currentMat = _contactMat.block(i, i, n, n);
            double countIntra = (currentMat.sum() - currentMat.diagonal().sum()) * .5;

            double sum1 = _contactMat.block(0, i, i, j-i+1).sum();
//      if (_VERBOSE)
//        std::cout << "sum1=" << sum1 << std::endl;

            double sum2;
            if (j + 1 < _N)
                sum2 = _contactMat.block(i, j+1, j-i+1, _N-1-j).sum();
            else
                sum2 = 0;
//      if (_VERBOSE)
//        std::cout << "sum2=" << sum2 << std::endl;
            double countInter = sum1 + sum2;
            _edgeCount(i, j) = countIntra;
            _edgeCount(j, i) = countInter;
        }
    }
    if (_VERBOSE) {
//    std::cout << "_edgeCount=" << _edgeCount << std::endl;
        std::cout << "finish data initialization; running time=" << (float)(std::clock()-t) / CLOCKS_PER_SEC << "s\n";;
        Writer::dumpMatrix(_edgeCount, _INPUT+".init1.txt");
    }
}
void Data::init2()
{
    std::clock_t t = std::clock();
    if (_VERBOSE)
        std::cout << "start data initialization\n";
    _edgeCount.resize(_N, _N);
//  _edgeCount.setZero();
    for (int k=1; k<_N; k++) {
        for (int i=0; i<_N-k; i++) {
            int j = i+k;
            double intra = _contactMat.coeff(i, j);
            if (j-1 > 0)
                intra += _edgeCount.coeff(i, j-1);
            if (i+1 < _N)
                intra += _edgeCount.coeff(i+1, j);
            if (j-1>0 && i+1 < _N)
                intra -= _edgeCount.coeff(i+1, j-1);
            if (abs(intra) < 1e-6) {
                _edgeCount(i, j) = 0;
                if (_VERBOSE)
                    printf("_edgeCount.coeff(%d, %d)=%f\n",i, j, _edgeCount.coeff(0, i-1)); fflush(stdout);
            }
            else
                _edgeCount(i, j) = intra;
        }
    }
    for (int i=0; i<_N; i++) {
        for (int j=i; j<_N; j++) {
            double inter = _edgeCount.coeff(0, j) + _edgeCount.coeff(i, _N-1) - 2*_edgeCount.coeff(i, j);
            if (i-1 > 0)
                inter -= _edgeCount.coeff(0, i-1);
            if (j+1 < _N)
                inter -= _edgeCount.coeff(j+1, _N-1);
            if (abs(inter) < 1e-6 || inter < 0) {
                _edgeCount(j, i) = 0;
                if (_VERBOSE) {
                    printf("i=%d, j=%d", i, j);
                    if (i-1 > 0)
                        printf(", _edgeCount.coeff(%d, %d)=%f", 0, i - 1, _edgeCount.coeff(0, i - 1));
                    if (j+1 < _N)
                        printf(", _edgeCount.coeff(%d, %d)=%f", j + 1, _N - 1, _edgeCount.coeff(j + 1, _N - 1));
                    std::cout << "\n";
                    fflush(stdout);
                }
            } else
                _edgeCount(j, i) = inter;
        }
    }
    if (_VERBOSE) {
//    std::cout << "_edgeCount=" << _edgeCount << std::endl;
        std::cout << "finish data initialization; running time=" << (float)(std::clock()-t) / CLOCKS_PER_SEC << "s\n";;
        Writer::dumpMatrix(_edgeCount, _INPUT+".init2.txt");
    }
}

double Data::getVol (int s, int e)
{
    if (s != e)
        return 2. * _edgeCount.coeff(s, e) + _edgeCount.coeff(e, s);
    else
        return _edgeCount.coeff(e, s);
}
