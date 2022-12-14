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
        Data *_data=NULL;
        Tree *_binaryTree=NULL;
        std::vector<TreeNode *> *_nodeList=NULL;
        double ***_table=NULL;
        int ***_minIndexTableForBold=NULL, ***_minIndexTable=NULL, ***_leftKtable=NULL;
        std::vector<Boundary> _boundaries;
        std::vector<TreeNode*> _trueNodeList;
        int _numBins=0, _kTmpIdx=0, _kMinusKtmpIdx=0;   // split into kTmp nodes and k-kTmp nodes
        float *_scoreTable=NULL;  // for tree pruning
        bool _breakFlag;

    public:
        Detector(SuperTAD::Data &data);

        ~Detector();

        void execute();

        void executeFilter(std::string result);

        void init();

        void filter();

        void fillTable();

        static bool sortByStart(Boundary a, Boundary b);

        void setIndex(int k, int &kIdx) { kIdx = k - 1; };

        void setNumBins(int s, int e) { _numBins = e - s + 1;};

        void backTrace(int k, bool add = false);

        void binarySplit(int s, int e, int k, bool add=false, int lv=0);

        void calculateD(TreeNode &node);

        void calculateDensity(TreeNode &node);

        double minusParent(double d, TreeNode &node);

        // filter nodes
        int *filterNodes();

        double getX(TreeNode &node);

        double getY(TreeNode &node);

        bool simpleLinearRegression(std::vector<std::pair<int, TreeNode*>> &nodeList, double ab[]);
    };
}


#endif //PROGRAM_DETECTORBINARY_H
