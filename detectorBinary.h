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
#include <string.h>


namespace SuperTAD::binary {

    class Detector {
    private:
        Data *_data;
        Tree *_binaryTree;
        std::vector<TreeNode *> *_nodeList;
        double ***_table;
        int ***_minIndexTableForBold;
        int ***_minIndexArray;
        int ***_leftKArray;
        std::vector<Boundary> _boundaries;
        std::vector<TreeNode*> _trueNodeList;
        int *_numBins;
        int *_kTmpIdx;
        int *_kMinusKtmpIdx;
        float *_scoreTable;

    public:
        Detector(SuperTAD::Data &data);

        ~Detector();

        void execute();

        void executeFILTER(std::string result);

        void init();

        void filter();

        void fillTable();

        static bool sortStart (Boundary a, Boundary b);

        void indexKtmp(int k) { *_kTmpIdx = k - 1; }

        void indexK(int k, int &kIdx) { kIdx = k - 1; }

        void numBins(int s, int e) { *_numBins = e - s + 1; }

        void backTrace(int k, bool add = false);

        void binarySplit(int s, int e, int k, bool add=false, int lv=0);

        void calculateD(TreeNode &node);

        void calculateDensity(TreeNode &node);

        double minusParent(double d, TreeNode &node);

        int * filterNodes();

        double getX(TreeNode &node);

        double getY(TreeNode &node);

        bool simpleLinearRegression(std::vector<std::pair<int, TreeNode*>> &nodeList, double ab[]);
    };
}


#endif //PROGRAM_DETECTORBINARY_H
