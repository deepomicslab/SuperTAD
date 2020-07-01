//
// Created by wang mengbo on 2019-09-01.
//

#ifndef PROGRAM_DETECTORBINARY_H
#define PROGRAM_DETECTORBINARY_H

#include <limits>
#include <algorithm>
#include <set>
#include "params.h"
#include "data.h"
#include "binaryTree.h"
#include "inputAndOutput.h"
#include "utils.h"


namespace binary {

    class Detector {
    private:
        Data *_data;
        Eigen::MatrixXd *_edgeCount;
        Writer _writer;
        binary::Tree *_binaryTree;
        std::vector<binary::TreeNode *> *_nodeList;
        double ***_table;
//        double **_baseTable;
//        double _gLogSum{};
        int ***_minIndexArray;
        int ***_leftKArray;
        std::vector<std::pair<int, int>> _boundary;
        std::set<binary::TreeNode *> _trueNodeList;
//        std::map<int, int> _kToIdx;
        int *_numBins;
        int *_kTmpIdx;
        int *_kMinusTmpIdx;

    public:
        Detector(Data &data);

        ~Detector();

        void execute();

        void init();

        void fillTable();

//        int indexK(int k) { return _kToIdx.find(k)->second; };
        void indexKtmp(int k) { *_kTmpIdx = k - 1; }
        void indexK(int k, int &kIdx) { kIdx = k - 1; }
//        int indexK(int k) { return k-1; }

        void numBins(int s, int e) { *_numBins = e - s + 1; }
//        int numBins(int &s, int &e) { return e-s+1; }

        void backTrace(int k, bool add = false);

        void binarySplit(int s, int e, int k, bool add=false, int lv=0);

        void calculateD(binary::TreeNode &node);

        void calculateDensity(binary::TreeNode &node);

        double minusParent(double d, binary::TreeNode &node);

        void filterNodes();

        double getX(binary::TreeNode &node);

        double getY(binary::TreeNode &node);

        bool simpleLinearRegression(std::vector<std::pair<int, binary::TreeNode *>> &nodeList, double ab[]);

//        void preSumGlog();
//
//        double getGlog(double binG);

        bool meaningfulComb(int s, int e, int k);

        int iForLastK(int s, int e, int k, int kTmp);

        int indexEnd(int s, int e) { return e-s; }
    };
}


#endif //PROGRAM_DETECTORBINARY_H
