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
        Eigen::MatrixXd * _edgeCount;
        Writer _writer;
        multi::Tree _multiTree;
        std::vector<multi::TreeNode *> * _nodeList;
        double **_table;
        int **_minIndexArray;
        int **_leftKArray;
        std::map<int, int> _kToIdx;
        std::vector<boundary> _boundaries;
        int _k;

    public:
        DetectorH1(Data &data);
        ~DetectorH1();
        void execute();
        void backTrace();
    };

};

#endif //PROGRAM_DETECTORH_H
