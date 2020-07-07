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


namespace multi {

//    class DetectorBase {
//    protected:
//        Data * _data;
//        Eigen::MatrixXd * _edgeCount;
//        multi::Tree _multiTree;
//        std::vector<multi::TreeNode *> * _nodeList;
//        std::map<int, int> _kToIdx;
//        std::vector<utils::boundary> _boundaryList;
//
//    public:
//        DetectorBase(Data & data);
//        ~DetectorBase() {};
//        int indexK(int k) { return _kToIdx.find(k)->second; }
//
//    };

    class Detector {
    private:
        Data * _data;
        Eigen::MatrixXd * _edgeCount;
        Writer _writer;
        multi::Tree _multiTree;
        std::vector<multi::TreeNode *> * _nodeList;
        double *****_table;
        int *****_minIndexArray;
        int *****_leftKArray;
        std::map<int, int> _kToIdx;
        std::vector<boundary> _boundaries;

    public:
        Detector(Data &data);
        ~Detector ();
        int indexK(int k) { return _kToIdx.find (k)->second; }
        void execute();
        void initK();
        void fillTable();
        void initH(int h);
        void backTrace(int k, int h, bool add = false);
        void multiSplit(int start, int end, int k, int h, int parentEnd, bool add = false);
    };


}

#endif //PROGRAM_DETECTORMULTI_H
