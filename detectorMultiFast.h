//
// Created by yuwzhang7 on 2021/10/18.
//

#ifndef PROGRAM_DETECTORMULTIFAST_H
#define PROGRAM_DETECTORMULTIFAST_H

#include <limits>
#include "params.h"
#include "data.h"
#include "inputAndOutput.h"
#include "multiTree.h"
#include "utils.h"
#include <algorithm>
#include <vector>


namespace SuperTAD::multifast {
    class Discretization{
    private:
        SuperTAD::Data *_data = NULL;
        double ***_tableF = NULL;
        double ****_tableT = NULL;
        int ***_minStepF = NULL;    // the searching result for each value in F
        int ****_minIndexT = NULL;  // the searching result for each value in T
        int ****_tableStepIndexT = NULL;    // the table storing the step index for each value in T
        double **_tablestep = NULL; // the table storing the step size given start and end
        double ****_sumofgcT = NULL;    // the table storing the sum of gc

        std::vector<Boundary> _boundaries;

    public:
        multi::Tree _multiTree;
        Discretization(SuperTAD::Data &data);
        ~Discretization();
        multi::Tree execute();
        void fillTable();
        void multiSplit(int start, int end, int l, int h);
    };

    class NeighborSearch{
    private:
        SuperTAD::Data *_data;
        multi::Tree _multiTree;
        std::vector<multi::TreeNode*> _nodeList;
        int _Nnodes = 0;
        SuperTAD::Writer _writer;
        double *****_table = NULL;
        int *****_minIndex = NULL;

    public:
        NeighborSearch(SuperTAD::Data &data, multi::Tree multiTree);
        ~NeighborSearch();
        void execute();
        void fillTable();
    };
}

#endif //PROGRAM_DETECTORMULTIFAST_H
