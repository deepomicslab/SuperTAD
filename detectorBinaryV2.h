//
// Created by yuwzhang7 on 2022/10/4.
//

#ifndef PROGRAM_DETECTORBINARYV2_H
#define PROGRAM_DETECTORBINARYV2_H

#include "data.h"
#include "detectorBinary.h"
#include "binaryTree.h"

namespace SuperTAD::binary {

    class BinaryV2 {
    private:
        Data *_data=NULL;
        SuperTAD::Writer _writer;
        double **_table=NULL;
        int **_minIndexTable=NULL;
        std::vector<Boundary> _boundaries;

    public:
        BinaryV2(SuperTAD::Data &data);

        ~BinaryV2();

        void execute();

        void init() {_boundaries.clear(); };

        void fillTable();

        void backTrace();

        void binarySplit(int s, int e);

    };
}

#endif //PROGRAM_DETECTORBINARYV2_H
