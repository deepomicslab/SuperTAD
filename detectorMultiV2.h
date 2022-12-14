//
// Created by yuwzhang7 on 2022/10/4.
//

#ifndef PROGRAM_DETECTORMULTIV2_H
#define PROGRAM_DETECTORMULTIV2_H

#include <limits>
#include "data.h"
#include "multiTree.h"

namespace SuperTAD::multi {

    class MultiV2 {
    private:
        SuperTAD::Data *_data = NULL;
        SuperTAD::Writer _writer;
        double ****_table = NULL;
        int ****_minIndexArray = NULL;
        std::vector<Boundary> _boundaries;
        multi::Tree _multiTree;

    public:
        MultiV2(SuperTAD::Data &data);

        ~MultiV2 ();

        void execute();

        void fillTable();

        void backTrace(int h, bool add = true);

        void multiSplit(int start, int end, int h, int parentEnd, bool add = false);
    };

}

#endif //PROGRAM_DETECTORMULTIV2_H
