//
// Created by yuwzhang7 on 2021/10/18.
//

#include "detectorMultiFast.h"

namespace SuperTAD { namespace multi{

    Discretization::Discretization(SuperTAD::Data &data)
    {
        _data = &data;
        SuperTAD::_H_ += 1;
        _multiTree.setData(data);
        if (SuperTAD::_STEP_ < 0)
            SuperTAD::_STEP_ = SuperTAD::_N_;   // default setting

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

    void Discretization::execute()
    {
        // fill table T and F
        fillTable();

        multiSplit(0, SuperTAD::_N_-1, _minStepF[0][SuperTAD::_N_-1][SuperTAD::_H_-1], SuperTAD::_H_-1);
    }

    void Discretization::fillTable()
    {
        if (SuperTAD::_VERBOSE_)
            printf("Start filling db table. The step of discretization varies.\n");
        // init the base case of F (H = 0)
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
        if (SuperTAD::_VERBOSE_)
            printf("base case with height = 1, se = %f\n", _tableF[0][SuperTAD::_N_-1][0]);

        // start to fill table T
        for (int h = 1; h < SuperTAD::_H_; h++) {
            // initialization
            for (int s = 0; s < SuperTAD::_N_; s++) {
                for (int e = s; e < SuperTAD::_N_; e++) {
                    _tableF[s][e][h] = std::numeric_limits<double>::infinity();
                    for (int k = 0; k < SuperTAD::_STEP_ + 2; k++) {
                        _tableT[s][e][k][h] = std::numeric_limits<double>::infinity();
                    }
                }
            }
            // upper case, H > 0
            for (int s = 0; s < SuperTAD::_N_; s++) {
                for (int i = s; i < SuperTAD::_N_; i++)
                {
                    for (int e = i; e < SuperTAD::_N_; e++)
                    {
                        if (_data->_logVolTable[s][e - s] <= 1e-15 or _tablestep[s][e-s] <= 1e-15)
                        {
                            _tableF[s][e][h] = 0;
                            _tableT[s][e][0][h] = 0;
                            _tableStepIndexT[s][e][0][h] = 0;
                            _minIndexT[s][e][0][h] = s;
                            _sumofgcT[s][e][0][h] = 0;
                            _minStepF[s][e][h] = 0;
                        }
                        else
                        {
                            double minTmp, sumofg, minTmpF;
                            int stepIndex = 0;
                            if (e == i)
                            {
                                stepIndex = int(_data->getG(s, e) / _tablestep[s][e - s]);

                                minTmp = _tableF[s][e][h - 1] - ((_data->getG(s, e) * _data->_logVolTable[s][e - s]) /_data->_doubleEdgeSum);
                                if (minTmp < _tableT[s][e][stepIndex][h]) {
                                    _tableT[s][e][stepIndex][h] = minTmp;
                                    _minIndexT[s][e][stepIndex][h] = i;
                                    _tableStepIndexT[s][e][stepIndex][h] = stepIndex;
                                    _sumofgcT[s][e][stepIndex][h] = _data->getG(s, e);  // store the accurate value
//                                    _tableF[s][e][h] = _tableT[s][e][stepIndex][h] + ((_data->getG(s, e)*
//                                            _data->_logVolTable[s][e-s]) / _data->_doubleEdgeSum);    // ???????? check.
//                                    _minStepF[s][e][h] = stepIndex;   // ???????? check.
                                }
                            }
                            else
                            {
                                for (int l = 0; l < SuperTAD::_STEP_ + 1; l++)
                                {
                                    if (not std::isinf(_tableT[s][i][l][h]))    // only consider the possible state in T
                                    {
                                        minTmp = _tableT[s][i][l][h] + _tableF[i + 1][e][h - 1] -
                                                ((_data->getG(i + 1, e) * _data->_logVolTable[i + 1][e-i-1]) / _data->_doubleEdgeSum);
                                        sumofg = _sumofgcT[s][i][l][h] + _data->getG(i + 1, e);
//                                        if (_tablestep[s][e - s] > SuperTAD::_THRESHOLD_)
                                        stepIndex = int(sumofg / _tablestep[s][e - s]);

                                        if (minTmp < _tableT[s][e][stepIndex][h])
                                        {
                                            _tableT[s][e][stepIndex][h] = minTmp;
                                            _minIndexT[s][e][stepIndex][h] = i;
                                            _tableStepIndexT[s][e][stepIndex][h] = l;
                                            _sumofgcT[s][e][stepIndex][h] = sumofg;
                                            minTmpF = _tableT[s][e][stepIndex][h] + ((_sumofgcT[s][e][stepIndex][h] *
                                                    _data->_logVolTable[s][e-s]) / _data->_doubleEdgeSum);
//                                            minTmpF = _tableT[s][e][stepIndex][h] + (stepIndex * _tablestep[s][e-s]) *
//                                                    (_data->_logVolTable[s][e-s] / _data->_doubleEdgeSum);  // check!
                                            if (minTmpF < _tableF[s][e][h])
                                            {
                                                _tableF[s][e][h] = minTmpF;
                                                _minStepF[s][e][h] = stepIndex;
//                                                if (SuperTAD::_VERBOSE_ && s == 0 && e == SuperTAD::_N_ - 1)
//                                                    printf("upper case, start = %d end = %d se = %f step = %d i = %d\n", s, e, minTmpF, stepIndex, i);
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
//        printf("start = %d, end= %d, l=%d, h=%d, mid=%d\n", start, end, l, h, mid_pos);
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

    NeighborSearch::NeighborSearch(SuperTAD::Data &data, SuperTAD::multi::Tree *multiTree, double ***tableF)
    {
        _data = &data;
        _multiTree = *multiTree;
        _tableF = tableF;

        _nodeList = _multiTree._nodeList;

        _multiTree.getPreOrder();
        _nodeList = _multiTree._nodeList;

        std::reverse(_nodeList.begin(), _nodeList.end());
        _Nnodes = _nodeList.size();  // containing root

        for (int i = 0; i < _Nnodes; i++)
            _nodeList[i]->_index = i;

        _table = new double ****[_Nnodes];
        _minIndex = new int ****[_Nnodes];

        for (int u = 0; u < _Nnodes; u++) {
            _table[u] = new double ***[2 * SuperTAD::_WINDOW_ + 1];
            _minIndex[u] = new int ***[2 * SuperTAD::_WINDOW_ + 1];

            for (int s = 0; s < 2 * SuperTAD::_WINDOW_ + 1; s++) {
                _table[u][s] = new double **[2 * SuperTAD::_WINDOW_ + 1];
                _minIndex[u][s] = new int **[2 * SuperTAD::_WINDOW_ + 1];

                for (int e = 0; e < 2 * SuperTAD::_WINDOW_ + 1; e++) {
                    _table[u][s][e] = new double *[_Nnodes];
                    _minIndex[u][s][e] = new int *[_Nnodes];

                    for (int c = 0; c < _Nnodes; c++) {
                        _table[u][s][e][c] = new double [2 * SuperTAD::_WINDOW_ + 1]{};
                        _minIndex[u][s][e][c] = new int [2 * SuperTAD::_WINDOW_ + 1]{};
                    }
                }
            }
        }
    }

    NeighborSearch::~NeighborSearch()
    {
        for (int u = 0; u < _Nnodes; u++) {
            for (int s = 0; s < 2 * SuperTAD::_WINDOW_ + 1; s++) {
                for (int e = 0; e < 2 * SuperTAD::_WINDOW_ + 1; e++) {
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
        for (int i = 0; i < _Nnodes; i++) {
            printf("(%d, %d)", _nodeList[i]->_val[0], _nodeList[i]->_val[1]);
            if (i < _Nnodes-1)
                printf(", ");
        }
        printf("\nThe number of tree nodes is %d\n", _Nnodes);

        NeighborSearch::fillTable();

        multi::TreeNode *root = _nodeList[_Nnodes-1];
        NeighborSearch::multiSplit(*root, 0, 0);

        if (SuperTAD::_VERBOSE_)
        {
            for (int i = 0; i < _Nnodes; i++) {
                printf("(%d, %d)", _nodeList[i]->_val[0], _nodeList[i]->_val[1]);
                if (i < _Nnodes-1)
                    printf(", ");
            }
            printf("\n");
        }
        _nodeList.pop_back();
        _writer.writeTree(SuperTAD::_OUTPUT_ + ".multifast", _nodeList);

    }

    void NeighborSearch::fillTable()
    {
        if (SuperTAD::_VERBOSE_)
            printf("Start filling db table in the neighbour searching step.\n");

        // initialization
        for (int s = 0; s < _Nnodes; s++) {
            for (int i = 0; i < 2 * SuperTAD::_WINDOW_ + 1; i++) {
                for (int j = 0; j < 2 * SuperTAD::_WINDOW_ + 1; j++) {
                    for (int e = 0; e < _Nnodes; e++) {
                        for (int p = 0; p < 2 * SuperTAD::_WINDOW_ + 1; p++) {
                            _table[s][i][j][e][p] = std::numeric_limits<double>::infinity();

                        }
                    }
                }
            }
        }

        int s, e;   //the start and end (bin) of each node
        int start, end; // the indicated bin after embedding window
        int child_s, child_e;   // the start and end (bin) of child node
        int child_e_win, child_end; // win_e_index, the indicated bin of child bound after embedding window
        int child_s_win, child_start;
        double minTmp;
        for (int i = 0; i < _Nnodes; i ++ )
        {
            s = _nodeList[i]->_val[0];
            e = _nodeList[i]->_val[1];
            for (int s_win = -SuperTAD::_WINDOW_; s_win < SuperTAD::_WINDOW_ + 1; s_win ++)
            {
                start = s + s_win;
                for (int e_win = -SuperTAD::_WINDOW_; e_win < SuperTAD::_WINDOW_+1; e_win ++)
                {
                    end = e + e_win;
                    // init for leaf se
                    if (start > -1 && end < SuperTAD::_N_ && start <= end)
                    {
                        if (_nodeList[i]->_children.size() == 0)
                            _table[i][mapWin(s_win)][mapWin(e_win)][0][0] = _tableF[start][end][0];
                        else
                        {
                            for (int m = _nodeList[i]->_children.size()-1; m >= 0; m--)    // child node index, from left to right
                            {
                                child_s = _nodeList[i]->_children[m]->_val[0];
                                child_e = _nodeList[i]->_children[m]->_val[1];
                                for (int c_e = -SuperTAD::_WINDOW_; c_e < SuperTAD::_WINDOW_ + 1; c_e++)
                                {
                                    if (m == 0)    // rightmost child
                                    {
                                        child_end = end;    // same right bound as parent
                                        child_e_win = 0;
                                    } else
                                    {
                                        child_end = child_e + c_e;  // modified right bound
                                        child_e_win = c_e;
                                    }

                                    for (int c_s = -SuperTAD::_WINDOW_; c_s < SuperTAD::_WINDOW_ + 1; c_s++)
                                    {
                                        if (m == _nodeList[i]->_children.size()-1)  // leftmost child
                                        {
                                            child_start = start;    // same left bound as parent
                                            child_s_win = 0;
                                        } else
                                        {
                                            child_start = child_s + c_s;    // modified left bound
                                            child_s_win = c_s;
                                        }
                                        if (child_start <= child_end && child_start > -1 && child_end < SuperTAD::_N_)
                                        {
                                            if (m == _nodeList[i]->_children.size()-1)  // leftmost child
                                            {
                                                _table[i][mapWin(s_win)][mapWin(e_win)][0][mapWin(child_e_win)] =
                                                        _data->getSE(child_start, child_end, start, end);
                                                if (_nodeList[i]->_children[m]->_children.size() == 0)
                                                {
                                                    _table[i][mapWin(s_win)][mapWin(e_win)][0][mapWin(child_e_win)] +=
                                                            _table[_nodeList[i]->_children[m]->_index][mapWin(
                                                                    child_s_win)][mapWin(child_e_win)][0][0];
                                                } else
                                                {
                                                    _table[i][mapWin(s_win)][mapWin(e_win)][0][mapWin(child_e_win)] +=
                                                            _table[_nodeList[i]->_children[m]->_index][mapWin(
                                                                    child_s_win)][mapWin(child_e_win)][
                                                                    _nodeList[i]->_children[m]->_children.size()-1]
                                                                    [mapWin(0)];
                                                }
                                                _minIndex[i][mapWin(s_win)][mapWin(e_win)][0][mapWin(child_e_win)] = child_s_win;

                                            } else
                                            {
                                                minTmp = _table[i][mapWin(s_win)][mapWin(e_win)][_nodeList[i]->_children.size() - m -2]
                                                        [mapWin(child_s_win)] + _data->getSE(child_start, child_end, start, end);
                                                if (_nodeList[i]->_children[m]->_children.size() == 0)
                                                    minTmp += _table[_nodeList[i]->_children[m]->_index][mapWin(
                                                            child_s_win)][mapWin(child_e_win)][0][0];
                                                else
                                                    minTmp += _table[_nodeList[i]->_children[m]->_index][mapWin(
                                                            child_s_win)][mapWin(child_e_win)][
                                                            _nodeList[i]->_children[m]->_children.size() - 1][mapWin(
                                                            0)];

                                                if (minTmp <
                                                    _table[i][mapWin(s_win)][mapWin(e_win)][_nodeList[i]->_children.size() - m -1][mapWin(child_e_win)])
                                                {
                                                    _table[i][mapWin(s_win)][mapWin(e_win)][_nodeList[i]->_children.size() - m -1][mapWin(
                                                            child_e_win)] = minTmp;
                                                    _minIndex[i][mapWin(s_win)][mapWin(e_win)][_nodeList[i]->_children.size() - m -1][mapWin(
                                                            child_e_win)] = child_s_win;
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
        if (SuperTAD::_VERBOSE_)
            printf("final se after neighborSearch is %f\n", _table[_Nnodes-1][mapWin(0)][mapWin(0)][_multiTree._root->_children.size()-1][mapWin(0)]);
    }

    void NeighborSearch::multiSplit(SuperTAD::multi::TreeNode &Node, int s_win, int e_win)
    {
        int right_pos_tmp = 0;
        for (int m = Node._children.size()-1; m >=0; m--)
        {
            multi::TreeNode *childNode = Node._children[Node._children.size()-1 -m];
            if (m == Node._children.size()-1)
            {
                right_pos_tmp = _minIndex[Node._index][mapWin(s_win)][mapWin(e_win)][Node._children.size()-1][mapWin(0)];
                childNode->_val[0] += right_pos_tmp;
                childNode->_val[1] += e_win;
                multiSplit(*childNode, right_pos_tmp, e_win);
            }
            else
            {
                int left_pos = _minIndex[Node._index][mapWin(s_win)][mapWin(e_win)][m][mapWin(right_pos_tmp)];
                if (m == 0)
                {
                    childNode->_val[0] += s_win;
                    childNode->_val[1] += right_pos_tmp;
                    multiSplit(*childNode, s_win, right_pos_tmp);
                    right_pos_tmp = left_pos;
                }
                else
                {
                    childNode->_val[0] += left_pos;
                    childNode->_val[1] += right_pos_tmp;
                    multiSplit(*childNode, left_pos, right_pos_tmp);
                    right_pos_tmp = left_pos;
                }
            }
        }
    }

} }

