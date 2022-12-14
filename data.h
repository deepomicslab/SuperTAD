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
        double **_contactTable=NULL;
        bool _initByPointer=false;

    public:
        double **_edgeCountTable=NULL;  // upper triangular is intra cluster conductance; lower triangular is inter clusters conductance
        double _edgeSum, _doubleEdgeSum, **_logVolTable=NULL;
        std::vector<double> _sumOfGtimesLogG;   // the sum of glog(g) of bins

        Data(std::string input);

        Data(double **array, int n);

        ~Data();

        void init();

        static void parseSubMatrix(double **&subMatrix, Data &Matrix, int start, int end);

        // return intra conductance(g) of node(s, e)
        double getG(int s, int e);

        // calculate and return volume of node(s, e)
        double getVol(int s, int e);

        // calculate and return  structure entropy given s and e; se=g/edge_sum*log2(V_p/V)
        double getSE(int s, int e, double parentVol);

        // calculate and return  structure entropy given s and e and parent volume; se=g/edge_sum*log2(V_p/V)
        double getSEwithLogPV(int s, int e, double logPV);

        // calculate and return  structure entropy given s, e, volumes of parent and current nodes; se=g/edge_sum*log2(V_p/V)
        double getSE(int s, int e, double parentVol, double currentVol);

        // calculate and return  structure entropy given s, e, ps, pe
        double getSE(int s, int e, int ps, int pe);

        // calculate and return  g*log(g)
        double getGtimesLogG(double binG);
    };
}

#endif //PROGRAM_DATA_H
