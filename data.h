//
// Created by wang mengbo on 2019-09-01.
//

#ifndef PROGRAM_DATA_H
#define PROGRAM_DATA_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <map>
#include <ctime>
#include "params.h"
#include "inputAndOutput.h"


class Data {
private:
    Eigen::MatrixXd _contactMat;
    bool _sym;

    // upper tri is intra; lower tri is inter
    Eigen::MatrixXd _edgeCount;

    double _edgeSum;
    double **** _asymEdgeCount;
    Reader *_reader;

public:
    Data(std::string fileName);
    ~Data();
//    void init0();
    void init();
    Eigen::MatrixXd & edgeCount() { return _edgeCount; }
    double getVol(int s, int e);
    double getSE(int start, int end, double parentVol);
    void setEdgeSum() { _edgeSum = _edgeCount.coeff(0, _N-1); }
    double getEdgeSum() { return _edgeSum; }
};

#endif //PROGRAM_DATA_H
