//
// Created by yuwzhang7 on 2020/11/2.
//

#include "detectorDeepBinary.h"


namespace SuperTAD::deepBinary
{

    Detector::Detector(Data &data)
    {
        _data = &data;
        _table = new double *[_N_];
        _minIndexArray = new int *[_N_];
        for (int s = 0; s < _N_; s++)
        {
            _table[s] = new double[_N_]{};
            _minIndexArray[s] = new int[_N_]{};
        }
        _boundaries.reserve(2 * _N_);
        _binaryTree = new binary::Tree();
        _binaryTree->setData(*_data);
    }


    Detector::~Detector()
    {
        delete _binaryTree;
        for (int s = 0; s < _N_; s++)
        {
            delete _table[s];
            delete _minIndexArray[s];
        }
        delete _table;
        delete _minIndexArray;
    }


    void Detector::execute()
    {
        std::clock_t tTmp;

        fillTable();

        backTrace(true);

        _nodeList = &_binaryTree->_nodeList;
        if (!_NO_OUTPUT_) {
            Writer::writeTree(_OUTPUT_ + ".deepbinary", *_nodeList);
        }

        if (_VERBOSE_)
            printf("start pruning deep binary tree\n");
        else
            printf("prune deep binary tree\n");

        binary::BasePruner *p;
        if (_PRUNE_) {
            switch (_PRUNE_METHOD_) {
                case binary::PruneMethod1:
                    p = new binary::Pruner1(*_binaryTree);
                    break;
                case binary::PruneMethod2:
                    p = new binary::Pruner2(*_binaryTree);
                    break;
                default:
                    p = new binary::Pruner2(*_binaryTree);

            }
            p->execute();
            multi::treeNodeVerbose(*(p->_prunedTree._root), 0);
            if (_VERBOSE_)
                printf("finish pruning deep binary tree\n");
        }


        if (!_NO_OUTPUT_) {

            if (_PRUNE_)
                Writer::writeTree(_OUTPUT_ + ".deepbinary.pruned", p->_prunedTree._nodeList);
            if (_SE_RESULT_PATH_!="") {
                std::ofstream outfile;
                if (_APPEND_RESULT_)
                    outfile.open(_SE_RESULT_PATH_, std::ios_base::app); // append instead of overwrite
                else
                    outfile.open(_SE_RESULT_PATH_);
                outfile << _table[0][_N_-1];
                if (_PRUNE_)
                    outfile << "\t" << p->_minHtable[p->_optimalK-1][p->_tree->_root->_idx];
                outfile << "\n";
                outfile.close();
            }
        }
        return;
    }


    void Detector::fillTable()
    {
        std::clock_t t;
        if (_VERBOSE_)
        {
            printf("start filling dp table\n");
            t = std::clock();
        }
        else
            printf("fill dp table\n");

//        for (int s = 0; s < _N_; s++)
//        {
//            printf("s: %d, _table[s][s]: %f, zero: %d\n", s, _table[s][s], _table[s][s] == 0);
//            _table[s][s] = 0;
//        }

        double parentVol, tmpSE, minSE;
        int leftE;
        for (int step = 1; step < _N_; step++)
        {
//            if (_VERBOSE_)
//                printf("step = %d\n", step);
            for (int s = 0; s < _N_ - step; s++)
            {
                parentVol = _data->getVol(s, s + step);
                minSE = std::numeric_limits<double>::infinity();
                for (int i = s; i < s + step; i++)
                {
                    tmpSE = _table[s][i] + _table[i + 1][s + step];
                    tmpSE += _data->getSE(s, i, parentVol);
                    tmpSE += _data->getSE(i + 1, s + step, parentVol);
                    if (tmpSE < minSE)
                    {
                        minSE = tmpSE;
                        leftE = i;
                    }
                }
                _table[s][s + step] = minSE;
                _minIndexArray[s][s + step] = leftE;
//                printf("step: %d, s: %d, _table[%d][%d]: %f, i=%d\n", step, s, s, s+step, _table[s][s+step], leftE);
            }
        }
    }


    void Detector::backTrace(bool add)
    {
//        init();
        printf("the SE for optimal deepbinary coding tree is %f\n", _table[0][_N_-1]);
        binarySplit(0, _N_ - 1, add);

//        sort(_boundaries.begin(), _boundaries.end(), utils::cmpBoundary);

        if (_VERBOSE_) {
            printf("boundaries:");
            for (int i = 0; i < _boundaries.size(); i++) {
                printf("(%d, %d)", _boundaries[i].first, _boundaries[i].second);
                if (i < _boundaries.size()-1)
                    printf(", ");
            }
            printf("\n");
        }

    }


    void Detector::binarySplit(int s, int e, bool add)
    {
        if (add)
        {
            int label;
            if (e-s==0){
                label = 0;
                _binaryTree->add(s, e, label);
            }
            else {
                label = 1;
                _binaryTree->add(s, e, label);
            }
        }

        if (e - s == 0)
            return;
        else {
            int i = _minIndexArray[s][e];
            _boundaries.emplace_back(s, i);
            _boundaries.emplace_back(i + 1, e);
            binarySplit(s, i, add);
            binarySplit(i+1, e, add);
        }
    }

}
