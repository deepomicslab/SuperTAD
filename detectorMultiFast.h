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


namespace SuperTAD { namespace multi {

    class Discretization{
    private:
        SuperTAD::Data *_data = NULL;
        double ****_tableT = NULL;

        int ***_minStepF = NULL;    // the searching result for each value in F
        int ****_minIndexT = NULL;  // the searching result for each value in T
        int ****_tableStepIndexT = NULL;    // the table storing the step index for each value in T
        double **_tablestep = NULL; // the table storing the step size given start and end
        double ****_sumofgcT = NULL;    // the table storing the sum of gc (accurate value)

        std::vector<Boundary> _boundaries;

    public:
        double ***_tableF = NULL;   // Will be used in NeighborSearching
        multi::Tree _multiTree;
        Discretization(SuperTAD::Data &data);
        ~Discretization();

        void execute();
        void fillTable();
        void multiSplit(int start, int end, int l, int h);
    };

    class NeighborSearch{
    private:
        SuperTAD::Data *_data;
        multi::Tree _multiTree;
        double *** _tableF;
        std::vector<multi::TreeNode*> _nodeList;    // include root
        int _Nnodes = 0;
        SuperTAD::Writer _writer;
        double *****_table = NULL;
        int *****_minIndex = NULL;

    public:
        NeighborSearch(SuperTAD::Data &data, multi::Tree *multiTree, double ***tableF);
        ~NeighborSearch();
        void execute();
        void fillTable();
        int mapWin(int index) { return index + SuperTAD::_WINDOW_;};
        void multiSplit(multi::TreeNode &Node, int s_win, int e_win);
    };
} }

#endif //PROGRAM_DETECTORMULTIFAST_H
