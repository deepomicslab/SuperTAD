//
// Created by yuwzhang7 on 2022/10/3.
//

#ifndef PROGRAM_PARTITIONANDMERGE_H
#define PROGRAM_PARTITIONANDMERGE_H

#include "data.h"

namespace SuperTAD::multi
{

    class PartitionV2
    {
    private:
        SuperTAD::Data *_data;
        SuperTAD::Writer _writer;
        double *_table;
        int *_minIndexArray;
        std::vector <Boundary> _boundaries;
        int _k = 0;

    public:
        PartitionV2(SuperTAD::Data &data);

        ~PartitionV2();

        std::vector <Boundary> execute(int h = -1);

        void backTrace();
    };


    class MergeV2
    {
    private:
        SuperTAD::Data *_data;
        SuperTAD::Writer _writer;
        std::vector <Boundary> _preBoundaries;
        std::vector<double> _prenodeSE; // the sum of bin se for each leave
        double *_table;
        int *_minIndexArray;
        int N;
        std::vector <Boundary> _boundaries;
        int _k;

    public:
        MergeV2(SuperTAD::Data &data, std::vector <Boundary> &_preBoundList);

        ~MergeV2();

        std::vector <Boundary> execute(int h = -1);

        void backTrace();
    };
};

#endif //PROGRAM_PARTITIONANDMERGE_H
