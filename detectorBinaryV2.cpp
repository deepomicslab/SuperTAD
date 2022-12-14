//
// Created by yuwzhang7 on 2022/10/4.
//

#include "detectorBinaryV2.h"

namespace SuperTAD::binary
{

    BinaryV2::BinaryV2(SuperTAD::Data &data)
    {
        _data = &data;
        _table = new double *[SuperTAD::_N_];
        _minIndexTable = new int *[SuperTAD::_N_];
        for (int s = 0; s < SuperTAD::_N_; s++)
        {
            _table[s] = new double[SuperTAD::_N_]{};
            _minIndexTable[s] = new int[SuperTAD::_N_]{};
        }
        _boundaries.reserve(SuperTAD::_K_);
    }


    BinaryV2::~BinaryV2()
    {
        for (int s = 0; s < SuperTAD::_N_; s++)
        {
            delete _table[s];
            delete _minIndexTable[s];
        }
        delete _table;
        delete _minIndexTable;
    }

    void BinaryV2::execute()
    {
        fillTable();

        backTrace();

        _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".binary.original", _boundaries);

        if (SuperTAD::_PRUNE_)
        {
            binary::Detector db(*_data);
            db.executeFilter(SuperTAD::_OUTPUT_ + ".binary.original");
        }
    }

    void BinaryV2::fillTable()
    {
        for (int l = 1; l < SuperTAD::_N_; l++)
        {
            for (int s = 0; s + l < SuperTAD::_N_ ; s++)
            {
                // basic case with only one domain
                _table[s][s+l] = _data->getVol(s, s+l) * _data->_logVolTable[s][l];
                if (s == 0)
                    _table[s][s+l] -= _data->_sumOfGtimesLogG[s+l];
                else
                    _table[s][s+l] -= _data->_sumOfGtimesLogG[s+l] - _data->_sumOfGtimesLogG[s-1];

                _table[s][s+l] /= _data->_doubleEdgeSum;
                _minIndexTable[s][s+l] = -1;    // record only one domain

                // upper case
                double tmpSE;

                for (int i = s; i < s+l; i++)
                {
                    tmpSE = _table[s][i] + _table[i+1][s+l];
                    tmpSE += _data->getSE(s, i, s, s+l);
                    tmpSE += _data->getSE(i+1, s+l, s, s+l);
                    if (tmpSE - _table[s][s+l] < - SuperTAD::_THRESHOLD_)
                    {
                        _table[s][s+l] = tmpSE;
                        _minIndexTable[s][s+l] = i;
                    }
                }
            }

        }

        if (SuperTAD::_VERBOSE_)
            printf("finish filling dp table, the minimum se is %f\n", _table[0][SuperTAD::_N_-1]);

    }

    void BinaryV2::backTrace()
    {
        init();

        int i = _minIndexTable[0][SuperTAD::_N_-1];

        binarySplit(0, i);
        binarySplit(i+1, SuperTAD::_N_-1);

        for (int i = 0; i < _boundaries.size(); i++)
        {
            printf("%d, %d\n", _boundaries[i].first, _boundaries[i].second);
        }

    }

    void BinaryV2::binarySplit(int s, int e)
    {
        int i = _minIndexTable[s][e];

        _boundaries.emplace_back(s+1, e+1);

        if (i == -1 or s == e)
        {
            return;
        }
        else
        {

            binarySplit(s, i);

            binarySplit(i + 1, e);
        }
    }

}
