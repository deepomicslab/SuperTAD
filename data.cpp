//
// Created by wang mengbo on 2019-09-01.
//

#include "data.h"


Data::Data(std::string fileName)
{
    _reader = new Reader(fileName);
    _N_ = _reader->parse(_contactMat);
    std::cout << "#bins=" << _N_ << std::endl;
    if (_K_ < 0) {
        _K_ = sqrt(_N_) + 5;
        std::cout << "k=" << _K_ << "\n";
    }
}


Data::~Data()
{
    delete _reader;
}

void Data::init()
{
    if (_VERBOSE_)
        std::cout << "start data initialization\n";

    std::clock_t t = std::clock();

    // calculate edge count(sum)
    _edgeCount.resize(_N_, _N_);
//    _edgeCount.setZero();
    for (int k=1; k < _N_; k++) {
        for (int i=0; i < _N_ - k; i++) {
            int j = i+k;
            double intra = _contactMat.coeff(i, j);
            if (j-1 > 0)
                intra += _edgeCount.coeff(i, j-1);
            if (i+1 < _N_)
                intra += _edgeCount.coeff(i+1, j);
            if (j-1>0 && i+1 < _N_)
                intra -= _edgeCount.coeff(i+1, j-1);
            if (abs(intra) < _THRESHOLD_ || intra < 0)
                _edgeCount(i, j) = 0;
            else
                _edgeCount(i, j) = intra;
        }
    }
    for (int i=0; i < _N_; i++) {
        for (int j=i; j < _N_; j++) {
            double inter = _edgeCount.coeff(0, j) + _edgeCount.coeff(i, _N_ - 1) - 2 * _edgeCount.coeff(i, j);
            if (i-1 > 0)
                inter -= _edgeCount.coeff(0, i-1);
            if (j+1 < _N_)
                inter -= _edgeCount.coeff(j+1, _N_ - 1);
            if (abs(inter) < _THRESHOLD_ || inter < 0) {
                _edgeCount(j, i) = 0;
            } else
                _edgeCount(j, i) = inter;
        }
    }
    setEdgeSum();

    if (_DEBUG_) {
//    std::cout << "_edgeCount=" << _edgeCount << std::endl;
        std::cout << "finish calculating edge count(sum); running time=" << (float)(std::clock()-t) / CLOCKS_PER_SEC << "s\n";
//        Writer::dumpMatrix(_edgeCount, _INPUT+".init.txt");
    }

    // calculate sum of g*log(g)
    t = std::clock();
    _sumOfGtimesLogG.emplace_back( getGtimesLogG(_edgeCount.coeff(0,0)) );
    for (int i=1; i < _N_; i++) {
        _sumOfGtimesLogG.emplace_back( _sumOfGtimesLogG[i-1] + getGtimesLogG(_edgeCount.coeff(i,i)));
    }
    if (_DEBUG_)
        std::cout << "finish calculating sum of g*log(g); running time=" << (float)(std::clock()-t) / CLOCKS_PER_SEC << "s\n";

    if (_VERBOSE_)
        printf("finish initialization\n");
}


void Data::setEdgeSum() {
    _edgeSum = _edgeCount.coeff(0, _N_ - 1);
}


double Data::getVol(int s, int e)
{
    if (s!=e)
        return 2. * _edgeCount.coeff(s, e) + _edgeCount.coeff(e, s);
    else
        return _edgeCount.coeff(e, s);
}


double Data::getSE(int start, int end, double parentVol)
{
    // g / edge_sum * log2(V_p / V)
    double currentVol = getVol(start, end);
    if(currentVol > 0 && parentVol >= currentVol)
        return _edgeCount(end, start) / (2. * _edgeSum) * log2(parentVol / currentVol);
    return 0;
}


double Data::getSE(int start, int end, double parentVol, double currentVol) {
    // g / edge_sum * log2(V_p / V)
    if(currentVol > 0 && parentVol >= currentVol)
        return _edgeCount(end, start) / (2. * _edgeSum) * log2(parentVol / currentVol);
    return 0;
}


double Data::getGtimesLogG(double binG) {
    if (binG==0)
        return 0;
    else
        return binG * log2(binG);
}
