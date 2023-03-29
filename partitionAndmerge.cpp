//
// Created by yuwzhang7 on 2022/10/3.
//

#include "partitionAndmerge.h"

namespace SuperTAD::multi
{

    PartitionV2::PartitionV2(SuperTAD::Data &data)
    {
        _data = &data;

        _table = new double [SuperTAD::_N_]{};
        _minIndexArray = new int [SuperTAD::_N_]{};
    }


    PartitionV2::~PartitionV2()
    {
        delete _table;
        delete _minIndexArray;
    }


    std::vector<Boundary> PartitionV2::execute(int h)
    {
        _table[0] = _data->getSE(0, 0, 0, SuperTAD::_N_ - 1);
        _minIndexArray[0] = -1;

        double binSum, tmpSE;

        for (int i = 1; i < SuperTAD::_N_; i++)
        {
            // base case, where 0-i belong to a domain
            _table[i] = _data->getSE(0, i, 0, SuperTAD::_N_ - 1);
            binSum = _data->getVol(0, i) * _data->_logVolTable[0][i] - _data->_sumOfGtimesLogG[i];
            _table[i] += binSum / _data->_doubleEdgeSum;
            _minIndexArray[i] = -1;

//            if (SuperTAD::_VERBOSE_) printf("finish calculating base case. se of %d is %f\n", i, _table[i]);

            // upper case
            for (int j = std::max(i-SuperTAD::_MAX_SIZE_, 0); j < i; j++)
            {
                tmpSE = _table[j] + _data->getSE(j+1, i, 0, SuperTAD::_N_-1);
                if (j+1 != i)
                {
                    binSum = _data->getVol(j+1, i) * _data->_logVolTable[j+1][i-j-1] - (_data->_sumOfGtimesLogG[i] - _data->_sumOfGtimesLogG[j]);
                    tmpSE += binSum / _data->_doubleEdgeSum;
                }

                if (tmpSE - _table[i] < -SuperTAD::_THRESHOLD_)
                {
                    _table[i] = tmpSE;
                    _minIndexArray[i] = j;
                }
            }
//            if (SuperTAD::_VERBOSE_) printf("finish i = %d, se is %f with %d\n", i, _table[i], _minIndexArray[i]);
        }

        if (SuperTAD::_VERBOSE_)
            printf("finish calculating h=1, the optimal se is %f\n", _table[SuperTAD::_N_-1]);

        PartitionV2::backTrace();

        if (h == -1)
            if (_BIN_LIST_)
                _writer.writeBoundaries(SuperTAD::_OUTPUT_ + ".multi2D.txt", _boundaries);
            else
                _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".multi2D", _boundaries);
        return _boundaries;
    }


    void PartitionV2::backTrace()
    {
        std::vector<int> bounds;
        int right = SuperTAD::_N_ - 1;
        while (right != -1)
        {
            bounds.emplace_back(right);
            right = _minIndexArray[right];
        }
        _boundaries.emplace_back(1, bounds[bounds.size()-1]+1);
        for (int i = bounds.size()-1; i > 0; i--)
            _boundaries.emplace_back(bounds[i]+2, bounds[i-1]+1);

//        for (int i = 0; i < _boundaries.size(); i++) {
//            std::cout << "boundary[" << i << "]=(" << _boundaries[i].first << ", " << _boundaries[i].second << ")\n";
//        }
    }

    MergeV2::MergeV2(SuperTAD::Data &data, std::vector<Boundary> &_preBoundList)
    {
        _data = &data;
        _preBoundaries = _preBoundList;
        N = (int) _preBoundaries.size();    // the number of boundaries
        _table = new double [N]{};
        _minIndexArray = new int [N]{};
    }

    MergeV2::~MergeV2()
    {
        delete _table;
        delete _minIndexArray;
    }

    std::vector<Boundary> MergeV2::execute(int h)
    {
        double binSum;
        int start, end;
        // record se of each Pre-node
        for (std::vector<Boundary>::iterator it = _preBoundaries.begin(); it != _preBoundaries.end(); ++it)
        {
            start = it->first - 1;
            end = it->second - 1;
            if (start == end)
                _prenodeSE.emplace_back(0);
            else
            {
                binSum = _data->getVol(start, end) * _data->_logVolTable[start][end - start];
                if (start == 0)
                    binSum -= _data->_sumOfGtimesLogG[end];
                else
                    binSum -= (_data->_sumOfGtimesLogG[end] - _data->_sumOfGtimesLogG[start - 1]);
                _prenodeSE.emplace_back(binSum / _data->_doubleEdgeSum);
//            printf("start=%d, end=%d, prenodeSE=%f\n", start, end, binSum/_data->_doubleEdgeSum);
            }
        }

        int nodeStart, nodeEnd;
        double tmpSE;

        for (int i = 0; i < N; i++)
        {
            // base case
            end = _preBoundaries[i].second - 1;
            _table[i] = _data->getSE(0, end, 0, SuperTAD::_N_-1);
            for (int j = 0; j < i + 1; j++)
            {
                _table[i] += _prenodeSE[j];
                nodeStart = _preBoundaries[j].first - 1;
                nodeEnd = _preBoundaries[j].second - 1;
                _table[i] += _data->getSE(nodeStart, nodeEnd, 0, end);
//                printf("nodestart=%d, nodeend=%d, nodevol=%f, _table=%f\n", nodeStart, nodeEnd, nodeVol, _table[i][0]);
            }
            _minIndexArray[i] = -1;

//            if (SuperTAD::_VERBOSE_) printf("finish calculating base case. #node=%d, se of %d is %f\n", N, i, _table[i]);
            // upper case
            for (int j = 0; j < i; j++)
            {
                tmpSE = _table[j];
                start = _preBoundaries[j+1].first - 1;
                tmpSE += _data->getSE(start, end, 0, SuperTAD::_N_-1);
                for (int l = j+1; l < i+1; l++)
                {
                    tmpSE += _prenodeSE[l];
                    nodeStart = _preBoundaries[l].first - 1;
                    nodeEnd = _preBoundaries[l].second - 1;
                    tmpSE += _data->getSE(nodeStart, nodeEnd, start, end);
                }

                if (tmpSE - _table[i] < - SuperTAD::_THRESHOLD_)
                {
                    _table[i] = tmpSE;
                    _minIndexArray[i] = j;
                }

            }
//            if (SuperTAD::_VERBOSE_) printf("finish calculating upper case. se of %d is %f with %d\n", i, _table[i], _minIndexArray[i]);

        }

        if (SuperTAD::_VERBOSE_)
            printf("finish calculating h=1, the optimal se is %f\n", _table[N-1]);

        MergeV2::backTrace();
        if (h == -1)
            _writer.writeBoundIn8Cols(SuperTAD::_OUTPUT_ + ".multi2D_Merge", _boundaries);
        return _boundaries;
    }


    void MergeV2::backTrace()
    {
        std::vector<int> bounds;
        int right = N - 1;
        while (right != -1)
        {
            bounds.emplace_back(right);
            right = _minIndexArray[right];
        }
        _boundaries.emplace_back(_preBoundaries[0].first, _preBoundaries[bounds[bounds.size()-1]].second);
        for (int i = bounds.size()-1; i > 0; i--)
        {
            _boundaries.emplace_back(_preBoundaries[bounds[i]].second+1, _preBoundaries[bounds[i - 1]].second);
        }
    }
}