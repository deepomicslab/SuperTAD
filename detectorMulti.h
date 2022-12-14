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
#include <vector>


namespace SuperTAD::multi {

    class Detector {
    private:
        SuperTAD::Data *_data = NULL;
        SuperTAD::Writer _writer;
        double *****_table = NULL;
        int *****_minIndexArray = NULL;
        int *****_leftKArray = NULL;
        std::vector<Boundary> _boundaries;

    public:
        multi::Tree _multiTree;
        explicit Detector(SuperTAD::Data &data);
        ~Detector ();

        void execute();
        void initK();
        void fillTable();
        void initH(int h);
        void backTrace(int k, int h, bool add = false);
        void multiSplit(int start, int end, int k, int h, int parentEnd, bool add = false);
    };

}

#endif //PROGRAM_DETECTORMULTI_H
