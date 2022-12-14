//
// Created by yuwzhang7 on 2022/10/4.
//

#include "detectorMultiV2.h"

namespace SuperTAD::multi
{

    MultiV2::MultiV2(SuperTAD::Data &data)
    {
        _data = &data;
        _multiTree.setData(data);

        _table = new double ***[SuperTAD::_N_];
        _minIndexArray = new int ***[SuperTAD::_N_];
        for (int s = 0; s < SuperTAD::_N_; s++)
        {
            _table[s] = new double **[SuperTAD::_N_];
            _minIndexArray[s] = new int **[SuperTAD::_N_];
            for (int e = s; e < SuperTAD::_N_; e++)
            {
                _table[s][e] = new double *[SuperTAD::_H_];
                _minIndexArray[s][e] = new int *[SuperTAD::_H_];
                for (int h = 0; h < SuperTAD::_H_; h++)
                {
                    _table[s][e][h] = new double[SuperTAD::_N_]{};
                    _minIndexArray[s][e][h] = new int[SuperTAD::_N_]{};
                }
            }
        }
    }


    MultiV2::~MultiV2()
    {
        for (int s = 0; s < SuperTAD::_N_; s++)
        {
            for (int e = s; e < SuperTAD::_N_; e++)
            {
                for (int h = 0; h < SuperTAD::_H_; h++)
                {
                    delete[] _table[s][e][h];
                    delete[] _minIndexArray[s][e][h];
                }
                delete[] _table[s][e];
                delete[] _minIndexArray[s][e];
            }
            delete[] _table[s];
            delete[] _minIndexArray[s];
        }
        delete[] _table;
        delete[] _minIndexArray;
    }

    void MultiV2::execute()
    {
        fillTable();
        _boundaries.clear();

        backTrace(SuperTAD::_H_);

        if (SuperTAD::_VERBOSE_) {
            printf("tree nodes:");
            for (int i = 0; i < _multiTree._nodeList.size(); i++) {
                printf("(%d, %d)", _multiTree._nodeList[i]->_val[0], _multiTree._nodeList[i]->_val[1]);
                if (i < _multiTree._nodeList.size()-1)
                    printf(", ");
            }
            printf("\n");
        }
        _writer.writeTree(SuperTAD::_OUTPUT_ + ".multi", _multiTree._nodeList);
    }


    void MultiV2::fillTable()
    {
        std::clock_t t = std::clock();
        if (SuperTAD::_VERBOSE_)
        {
            printf("start filling db table\n");
        }

        // base case, h=0
        double binSum, tmpSE;
        for (int l = 1; l < SuperTAD::_N_; l++)
        {
            for (int s = 0; s + l < SuperTAD::_N_; s++)
            {
                // basic case with only one domain
                binSum = _data->getVol(s, s + l) * _data->_logVolTable[s][l];
                if (s == 0)
                    binSum -= _data->_sumOfGtimesLogG[s + l];
                else
                    binSum -= _data->_sumOfGtimesLogG[s + l] - _data->_sumOfGtimesLogG[s - 1];

                for (int p = s + l; p < SuperTAD::_N_; p++)
                {
                    _minIndexArray[s][s][0][p] = -1;
                    _table[s][s + l][0][p] = binSum / _data->_doubleEdgeSum + _data->getSE(s, s+l, s, p);
                    _minIndexArray[s][s + l][0][p] = -1;

                    // multiple domains
                    for (int i = s; i < s + l; i++)
                    {
                        // left part
                        if (i == s)
                            tmpSE = _data->getSE(s, s, s, p);
                        else
                            tmpSE = _table[s][i][0][p];
                        // rightmost domain
                        if (i + 1 == s + l)
                            tmpSE += _data->getSE(i + 1, i + 1, s, p);
                        else
                        {
                            tmpSE += _data->getSE(i + 1, s + l, s, p);
                            tmpSE += (_data->getVol(i + 1, s + l) * _data->_logVolTable[i + 1][s + l - i - 1]) /
                                     _data->_doubleEdgeSum;
                            tmpSE -= (_data->_sumOfGtimesLogG[s + l] - _data->_sumOfGtimesLogG[i]) /
                                     _data->_doubleEdgeSum;
                        }
                        if (tmpSE - _table[s][s + l][0][p] < -SuperTAD::_THRESHOLD_)
                        {
                            _table[s][s + l][0][p] = tmpSE;
                            _minIndexArray[s][s + l][0][p] = i;
                        }

                    }
                }
            }
        }
        if (SuperTAD::_VERBOSE_)
            printf("finish filling the base case, table[0][N-1][h=0][N-1]=%f\n", _table[0][SuperTAD::_N_-1][0][SuperTAD::_N_-1]);

        for (int h = 1; h < SuperTAD::_H_; h ++)
        {
            // start == end
            for (int s = 0; s < SuperTAD::_N_; s++)
            {
                for (int p = s+1; p < SuperTAD::_N_; p++)
                {
                    _table[s][s][h][p] = _data->getSE(s, s, s, p);
                    _minIndexArray[s][s][h][p] = -1;
                }
            }
            // end > start
            for (int l = 1; l < SuperTAD::_N_; l++)
            {
                for (int s = 0; s+l < SuperTAD::_N_; s++)
                {
                    for (int p = s + l; p < SuperTAD::_N_; p++)
                    {
                        _table[s][s+l][h][p] = _table[s][s+l][h-1][s+l] + _data->getSE(s, s+l, s, p);
                        _minIndexArray[s][s+l][h][p] = -1;  // height minus one

                        for (int i = s; i < s+l; i++)
                        {
                            tmpSE = _table[s][i][h][p] + _table[i+1][s+l][h-1][s+l] + _data->getSE(i+1, s+l, s, p);
                            if (tmpSE - _table[s][s+l][h][p] < - SuperTAD::_THRESHOLD_)
                            {
                                _table[s][s+l][h][p] = tmpSE;
                                _minIndexArray[s][s+l][h][p] = i;
                            }
                        }
                    }
                }
            }
            if (SuperTAD::_VERBOSE_)
                printf("finish filling h=%d, table[0][N-1][%d][N-1]=%f\n", h, h, _table[0][SuperTAD::_N_-1][h][SuperTAD::_N_-1]);
        }
        if (SuperTAD::_VERBOSE_)
            printf("finish filling db table, consumes %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);
    }

    void MultiV2::backTrace(int h, bool add)
    {

        multiSplit(0, SuperTAD::_N_ - 1, SuperTAD::_H_ - 1, SuperTAD::_N_ - 1, add);

        _boundaries.emplace_back(0, 0);

        sort(_boundaries.begin(), _boundaries.end(), utils::cmpBoundary);

        for (int i = 0; i < _boundaries.size(); i++) {
            if (i == _boundaries.size() - 1) {
                _boundaries[i].second = SuperTAD::_N_ - 1;
            }
            else {
                _boundaries[i].second = _boundaries[i + 1].first - 1;
            }
        }

        if (SuperTAD::_VERBOSE_) {
            printf("leaves:");
            for (int i = 0; i < _boundaries.size(); i++) {
                printf("(%d, %d)", _boundaries[i].first, _boundaries[i].second);
                if (i < _boundaries.size()-1)
                    printf(", ");
            }
            printf("\n");
        }

    }

    void MultiV2::multiSplit (int start, int end, int h, int parentEnd, bool add)
    {
        int mid = _minIndexArray[start][end][h][parentEnd];

        if (mid == -1)
        {
            if (add)
                _multiTree.add(start, end);
            if (SuperTAD::_VERBOSE_) {
                printf("multisplit-------------%d %d %d height minus one, parentEnd: %d, se: %f\n",
                       start, end, h-1, parentEnd, _table[start][end][h][parentEnd]);
            }
            if (start != end and h > 0)
                multiSplit(start, end, h-1, end, add);
        } else
        {
            _boundaries.emplace_back(mid + 1, 0);
            if (add)
                _multiTree.add(mid + 1, end);
            if (SuperTAD::_VERBOSE_) {
                printf("multisplit-------------%d %d %d parentEnd: %d partition, mid: %d, se: %f\n",
                       start, end, h, parentEnd, mid, _table[start][end][h][parentEnd]);
            }
            if (h > 0)
            {
                multiSplit(start, mid, h, parentEnd, add);
                multiSplit(mid + 1, end, h-1, end, add);
            } else
                multiSplit(start, mid, h, parentEnd, add);
        }
    }

}

