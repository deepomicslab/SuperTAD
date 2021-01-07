//
// Created by mengbowang on 4/29/2020.
//

#ifndef PROGRAM_DETECTORH_H
#define PROGRAM_DETECTORH_H

#include <limits>
#include "params.h"
#include "data.h"
#include "inputAndOutput.h"
#include "multiTree.h"
#include <map>
#include "utils.h"
#include <ctime>
#include "inputAndOutput.h"

namespace SuperTAD::multi {

    class DetectorH1 {
    private:
        SuperTAD::Data * _data;
//        Eigen::MatrixXd * _edgeCount;
        SuperTAD::Writer _writer;
        double **_table;
        int **_minIndexArray;
//        std::map<int, int> _kToIdx;
        std::vector<Boundary> _boundaries;
        int _k = 1;

    public:
        DetectorH1(SuperTAD::Data &data);

        ~DetectorH1();

        std::vector<Boundary> execute(int h=-1);

        void backTrace();
    };

    
    class Merge {
    private:
        SuperTAD::Data * _data;
        SuperTAD::Writer _writer;
        std::vector<Boundary> _preBoundaries;
        std::vector<double> _prenodeSE;
        double **_table;
        int **_minIndexArray;
        int N;
        std::vector<Boundary> _boundaries;
        int _k;

    public:
        Merge(SuperTAD::Data &data, std::vector<Boundary> &_preBoundList);

        ~Merge();

        std::vector<Boundary> execute(int h=-1);

        void backTrace();
    };

    class detectorH {
    private:
        SuperTAD::Data * _data;
        SuperTAD::Writer _writer;
    public:
        std::vector<Boundary> _boundary;
        std::vector<Boundary> _clusters;
        detectorH(SuperTAD::Data &data);
        void pipeline(std::string preResult);
    };

};

#endif //PROGRAM_DETECTORH_H
