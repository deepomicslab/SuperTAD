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
#include <string>
#include <map>
#include <ctime>
#include <cmath>
#include "params.h"
#include "inputAndOutput.h"
#include "utils.h"


namespace SuperTAD
{
    class Data {
    private:
        std::string _chrom1;

        std::string _chrom2;

        std::map<int, std::pair<int64_t, int64_t>> _chrom1Idx2Interval;

        std::map<int, std::pair<int64_t, int64_t>> _chrom2Idx2Interval;

        double **_contactArray=NULL;

        bool _initByPointer=false;

    public:
        // upper tri is intra; lower tri is inter
        double **_edgeCountArray=NULL, _edgeSum, _doubleEdgeSum, **_logVolTable=NULL, **_volTable=NULL;
        std::vector<double> _sumOfGtimesLogG;

        Data(std::string input);

        Data(double **array, int n);

        ~Data();

        void init();

        static void parseSubMatrix(double **&subMatrix, Data &Matrix, int start, int end);

        // debug use
        double getG(int s, int e);

        double getVol(int s, int e);

        double getSE(int s, int e, double parentVol);

        double getSEwithLogPV(int s, int e, double logPV);

        double getSE(int s, int e, double parentVol, double currentVol);

        double getSEwithLogs(int s, int e, double logPV, double logCV);

        double getSEwithLogDiff(int s, int e, double logDiff);

        double getGtimesLogG(double binG);
    };
}

#endif //PROGRAM_DATA_H
