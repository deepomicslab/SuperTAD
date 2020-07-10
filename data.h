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
#include "utils.h"


class Data {
private:
    std::string _chrom1;

    std::string _chrom2;

    std::map<int, std::pair<int64_t, int64_t>> _chrom1Idx2Interval;

    std::map<int, std::pair<int64_t, int64_t>> _chrom2Idx2Interval;

//    Eigen::MatrixXd _contactMat;
    double **_contactArray;

    // upper tri is intra; lower tri is inter
//    Eigen::MatrixXd _edgeCountMat;

public:
    double **_edgeCountArray;

    double _edgeSum;

    double _doubleEdgeSum;

    double **_logVolTable;

    double **_volTable;

    std::vector<double> _sumOfGtimesLogG;

    Data(std::string fileName);

    ~Data();

    void init();

//    Eigen::MatrixXd & edgeCount() { return _edgeCountMat; }

    double getVol(int s, int e);

    double getSE(int s, int e, double parentVol);

    double getSEwithLogPV(int s, int e, double logPV);

    double getSE(int s, int e, double parentVol, double currentVol);

    double getSEwithLogs(int s, int e, double logPV, double logCV);

    double getSEwithLogDiff(int s, int e, double logDiff);

//    void setEdgeSum();

    double getGtimesLogG(double binG);
};

#endif //PROGRAM_DATA_H
