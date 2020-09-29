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

namespace multi {

    class DetectorH1 {
    private:
        Data * _data;
//        Eigen::MatrixXd * _edgeCount;
        Writer _writer;
        multi::Tree _multiTree;
        std::vector<multi::TreeNode *> * _nodeList;
        double **_table;
        int **_minIndexArray;
        std::map<int, int> _kToIdx;
        std::vector<Boundary> _boundaries;
        int _k;

    public:
        DetectorH1(Data &data);
        ~DetectorH1();
        void execute();
        void backTrace();
    };

    class Merge {
    private:
        Data * _data;
        std::vector<Boundary> _preBoundaries;
        std::vector<double> _prenodeSE;
        double **table;
        int **minIndexArray;
        int _k;
        std::vector<Boundary> _boundaries;

    public:
        Merge(Data &data, std::vector<Boundary> &_preBoundaries);
        ~Merge();
        void execute();
        void backTrace();
    };

    class detectorH {
    private:
        Data * _data;
        std::vector<Boundary> _boundary;
    public:
        detectorH(Data &data);
        ~detectorH();
        void pipeline(std::string preResult);
    };

};

#endif //PROGRAM_DETECTORH_H
