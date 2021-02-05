//
// Created by yuwzhang7 on 2020/11/2.
//

#ifndef PROGRAM_DETECTORDEEPBINARY_H
#define PROGRAM_DETECTORDEEPBINARY_H

#include <limits>
#include <algorithm>
#include <set>
#include <string.h>
#include "params.h"
#include "data.h"
#include "binaryTree.h"
#include "inputAndOutput.h"
#include "utils.h"


namespace SuperTAD::deepBinary {

    class Detector {
    private:
    public:
        SuperTAD::Data *_data;
        binary::Tree *_binaryTree;
        std::vector<binary::TreeNode *> *_nodeList;
        double **_table;
        int **_minIndexArray;
        std::vector<Boundary> _boundaries;
        binary::BasePruner *_pruner=NULL;
        
        Detector(SuperTAD::Data &data);

        ~Detector();

        void execute();

        void fillTable();

        void backTrace(bool add = false);

        void binarySplit(int s, int e, bool add=false);
    };

}

#endif //PROGRAM_DETECTORDEEPBINARY_H



