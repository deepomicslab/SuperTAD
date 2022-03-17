//
// Created by yuwzhang7 on 2021/10/18.
//

#include "detectorMultiFast.h"

namespace SuperTAD::multifast{

    Discretization::Discretization(SuperTAD::Data &data)
    {
        _data = &data;
        _multiTree.setData(data);
        if (SuperTAD::_STEP_ < 0)
            SuperTAD::_STEP_ = SuperTAD::_N_;
        _tableF = new double **[SuperTAD::_N_];
        _minStepF = new int **[SuperTAD::_N_];
        _tableT = new double ***[SuperTAD::_N_];
        _minIndexT = new int ***[SuperTAD::_N_];
        _tableStepIndexT = new int ***[SuperTAD::_N_];
        _sumofgcT = new double ***[SuperTAD::_N_];
        _tablestep = new double *[SuperTAD::_N_];
        for (int s = 0; s < SuperTAD::_N_; s++) {
            _tablestep[s] = new double[SuperTAD::_N_ - s]{};    // efficiently stored
            _tableF[s] = new double *[SuperTAD::_N_];
            _minStepF[s] = new int *[SuperTAD::_N_];
            _tableT[s] = new double **[SuperTAD::_N_];
            _minIndexT[s] = new int **[SuperTAD::_N_];
            _tableStepIndexT[s] = new int **[SuperTAD::_N_];
            _sumofgcT[s] = new double **[SuperTAD::_N_];
            for (int e = s; e < SuperTAD::_N_; e++) {
                _tableF[s][e] = new double [SuperTAD::_H_]{};
                _minStepF[s][e] = new int [SuperTAD::_H_]{};
                _tableT[s][e] = new double *[SuperTAD::_STEP_ + 2];
                _minIndexT[s][e] = new int *[SuperTAD::_STEP_ + 2];
                _tableStepIndexT[s][e] = new int *[SuperTAD::_STEP_+2];
                _sumofgcT[s][e] = new double *[SuperTAD::_STEP_ + 2];
                for (int i = 0; i < SuperTAD::_STEP_ + 2; i++) {
                    _tableT[s][e][i] = new double [SuperTAD::_H_]{};
                    _minIndexT[s][e][i] = new int [SuperTAD::_H_]{};
                    _tableStepIndexT[s][e][i] = new int [SuperTAD::_H_]{};
                    _sumofgcT[s][e][i] = new double [SuperTAD::_H_]{};
                }
            }
        }
    }

    Discretization::~Discretization()
    {
        for (int s = 0; s < SuperTAD::_N_; s++){
            for (int e = s; e < SuperTAD::_N_; e++){
                for (int i = 0; i < SuperTAD::_STEP_ + 2; i++){
                    delete []_tableT[s][e][i];
                    delete []_minIndexT[s][e][i];
                    delete []_tableStepIndexT[s][e][i];
                    delete []_sumofgcT[s][e][i];
                }
                delete []_tableF[s][e];
                delete []_minStepF[s][e];
                delete []_tableT[s][e];
                delete []_minIndexT[s][e];
                delete []_tableStepIndexT[s][e];
                delete []_sumofgcT[s][e];
            }
            delete []_tableF[s];
            delete []_minStepF[s];
            delete []_tableT[s];
            delete []_minIndexT[s];
            delete []_tableStepIndexT[s];
            delete []_sumofgcT[s];
        }
        delete []_tableF;
        delete []_minStepF;
        delete []_tableT;
        delete []_minIndexT;
        delete []_tableStepIndexT;
        delete []_sumofgcT;
    }

    multi::Tree Discretization::execute()
    {
        // fill table T and F
        fillTable();

        multiSplit(0, SuperTAD::_N_-1, _minStepF[0][SuperTAD::_N_-1][SuperTAD::_H_-1], SuperTAD::_H_-1);

        return _multiTree;
    }

    void Discretization::fillTable()
    {
        if (SuperTAD::_VERBOSE_)
            printf("Start filling db table. The step of discretization varies\n");
        // init the leaves of F (H=0)
        for (int s = 0; s < SuperTAD::_N_; s++ ){
            for (int e = s; e < SuperTAD::_N_; e++){
                double currentVolumn = _data->getVol(s, e);
                _tablestep[s][e-s] = currentVolumn / SuperTAD::_STEP_;
                _tableF[s][e][0] = currentVolumn * _data->_logVolTable[s][e-s];
                if (s == 0)
                    _tableF[s][e][0] -= _data->_sumOfGtimesLogG[e];
                else
                    _tableF[s][e][0] -= _data->_sumOfGtimesLogG[e] - _data->_sumOfGtimesLogG[s-1];
                _tableF[s][e][0] /= _data->_doubleEdgeSum;
            }
        }

        // start to fill table T
        for (int h = 1; h < SuperTAD::_H_; h++) {
            for (int s = 0; s < SuperTAD::_N_; s++) {
                for (int e = s; e < SuperTAD::_N_; e++) {
                    _tableF[s][e][h] = std::numeric_limits<double>::infinity();
                    for (int k = 0; k < SuperTAD::_STEP_ + 2; k++) {
                        _tableT[s][e][k][h] = std::numeric_limits<double>::infinity();
                    }
                }
                for (int i = s; i < SuperTAD::_N_; i++)
                {
                    for (int e = i; e < SuperTAD::_N_; e++)
                    {
                        if (std::isnan(_data->_logVolTable[s][e - s]))
                        {
                            _tableF[s][e][h] = 0;
                            _tableT[s][e][0][h] = 0;
                        } else
                        {
                            double minTmp, sumofg;
                            int stepIndex = 0;
                            if (e == i)
                            {
                                if (_tablestep[s][e-s] > SuperTAD::_THRESHOLD_)
                                    stepIndex = int(_data->getG(s, e) / _tablestep[s][e - s]);

                                _tableT[s][e][stepIndex][h] = _tableF[s][e][h - 1] -
                                                              (_data->getG(s, e) * _data->_logVolTable[s][e - s] /
                                                               _data->_doubleEdgeSum);
                                _minIndexT[s][e][stepIndex][h] = i;
                                _tableStepIndexT[s][e][stepIndex][h] = stepIndex;
                                _minStepF[s][e][h] = stepIndex;
                                _sumofgcT[s][e][stepIndex][h] = _data->getG(s, e);
                            } else
                            {
                                for (int l = 0; l < SuperTAD::_STEP_ + 1; l++)
                                {
                                    if (not std::isinf(_tableT[s][i][l][h]))
                                    {
                                        minTmp = _tableT[s][i][l][h] + _tableF[i + 1][e][h - 1] -
                                                (_data->getG(i + 1, e) * _data->_logVolTable[i + 1][e-i-1]) / _data->_doubleEdgeSum;
                                        sumofg = _sumofgcT[s][i][l][h] + _data->getG(i + 1, e);
                                        if (_tablestep[s][e-s] > SuperTAD::_THRESHOLD_)
                                            stepIndex = int(sumofg / _tablestep[s][e - s]);

                                        if (minTmp < _tableT[s][e][stepIndex][h])
                                        {
                                            _tableT[s][e][stepIndex][h] = minTmp;
                                            _minIndexT[s][e][stepIndex][h] = i;
                                            _tableStepIndexT[s][e][stepIndex][h] = l;
                                            _sumofgcT[s][e][stepIndex][h] = sumofg;
                                            minTmp = _tableT[s][e][stepIndex][h] +
                                                     (_sumofgcT[s][e][stepIndex][h] * _data->_logVolTable[s][e-s] /
                                                      _data->_doubleEdgeSum);
                                            if (minTmp < _tableF[s][e][h])
                                            {
                                                _tableF[s][e][h] = minTmp;
                                                _minStepF[s][e][h] = stepIndex;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }

       }

    }

    void Discretization::multiSplit(int start, int end, int l, int h)
    {
        int mid_pos = _minIndexT[start][end][l][h];
        printf("start = %d, end= %d, l=%d, h=%d, mid=%d\n", start, end, l, h, mid_pos);
        if (mid_pos == end) {
            _multiTree.add(start, end);
            if (h > 1) {
                multiSplit(start, end, _minStepF[start][end][h-1], h-1);
            }
        }
        else {
            _multiTree.add(mid_pos+1, end);
            multiSplit(start, mid_pos, _tableStepIndexT[start][end][l][h], h);
            if (h > 1) {
                multiSplit(mid_pos+1, end, _minStepF[mid_pos+1][end][h-1], h-1);
            }
        }
    }

    NeighborSearch::NeighborSearch(SuperTAD::Data &data, SuperTAD::multi::Tree multiTree)
    {
        _data = &data;
        _multiTree = multiTree;
        _Nnodes = _multiTree._nodeList.size();  // containing root
        printf("#node=%d\n", _Nnodes);
        for (int i = 0; i < _Nnodes; i++) {
                printf("(%d, %d)", _multiTree._nodeList[i]->_val[0], _multiTree._nodeList[i]->_val[1]);
                if (i < _Nnodes-1)
                    printf(", ");
            }
        printf("%d\n", _Nnodes);
        _table = new double ****[_Nnodes];
        _minIndex = new int ****[_Nnodes];
        for (int u = 0; u < _Nnodes; u++) {
            _table[u] = new double ***[SuperTAD::_WINDOW_];
            _minIndex[u] = new int ***[SuperTAD::_WINDOW_];
            for (int s = 0; s < SuperTAD::_WINDOW_; s++) {
                _table[u][s] = new double **[SuperTAD::_WINDOW_];
                _minIndex[u][s] = new int **[SuperTAD::_WINDOW_];
                for (int e = 0; e < SuperTAD::_WINDOW_; e++) {
                    _table[u][s][e] = new double *[_Nnodes];
                    _minIndex[u][s][e] = new int *[_Nnodes];
                    for (int c = 0; c < _Nnodes; c++) {
                        _table[u][s][e][c] = new double [SuperTAD::_WINDOW_]{};
                        _minIndex[u][s][e][c] = new int [SuperTAD::_WINDOW_]{};
                    }
                }
            }
        }
    }

    NeighborSearch::~NeighborSearch()
    {
        for (int u = 0; u < _Nnodes; u++) {
            for (int s = 0; s < SuperTAD::_WINDOW_; s++) {
                for (int e = 0; e < SuperTAD::_WINDOW_; e++) {
                    for (int c = 0; c < _Nnodes; c++) {
                        delete []_table[u][s][e][c];
                        delete []_minIndex[u][s][e][c];
                    }
                    delete []_table[u][s][e];
                    delete []_minIndex[u][s][e];
                }
                delete []_table[u][s];
                delete []_minIndex[u][s];
            }
            delete []_table[u];
            delete []_minIndex[u];
        }
        delete []_table;
        delete []_minIndex;
    }

    void NeighborSearch::execute()
    {
        _multiTree.getNodeList(_nodeList);  // preorder traversal, including root.
        _Nnodes = _nodeList.size();
        for (int i = 0; i < _Nnodes; i++) {
            printf("(%d, %d)", _nodeList[i]->_val[0], _nodeList[i]->_val[1]);
            if (i < _Nnodes-1)
                printf(", ");
        }
        printf("%d\n", _Nnodes);



    }

    void NeighborSearch::fillTable()
    {
        if (SuperTAD::_VERBOSE_)
            printf("Start filling db table in the neighbour searching step.\n");
    }
}

