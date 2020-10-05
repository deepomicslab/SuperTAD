//
// Created by wang mengbo on 2019-09-03.
//

#ifndef PROGRAM_DETECTORMULTI_H
#define PROGRAM_DETECTORMULTI_H

#include <limits>
#include "params.h"
#include "data.h"
#include "inputAndOutput.h"
#include "multiTree.h"
#include <map>
#include "utils.h"
#include <ctime>
#include <algorithm>


namespace multi {

    class Detector {
    private:
        SuperTAD::Data * _data;
//        Eigen::MatrixXd * _edgeCount;
        SuperTAD::Writer _writer;
        multi::Tree _multiTree;
        std::vector<multi::TreeNode *> * _nodeList;
        double *****_table;
        int *****_minIndexArray;
        int *****_leftKArray;
//        std::map<int, int> _kToIdx;
        std::vector<Boundary> _boundaries;

    public:
        Detector(SuperTAD::Data &data);
        ~Detector ();
//        int indexK(int k) { return _kToIdx.find(k)->second; }
        void execute();
        void initK();
        void fillTable();
        void initH(int h);
        void backTrace(int k, int h, bool add = false);
        void multiSplit(int start, int end, int k, int h, int parentEnd, bool add = false);
    };

}

#endif //PROGRAM_DETECTORMULTI_H
