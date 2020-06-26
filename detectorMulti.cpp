//
// Created by wang mengbo on 2019-09-03.
//

#include "detectorMulti.h"

namespace multi {

    Detector::Detector(Data &data)
    {
        _data = &data;
        _edgeCount = &data.edgeCount();

        int k = 1;
        for (int i = 0; i < _K; i++) {
            _kToIdx.emplace(k++, i);
        }
        _boundaries.reserve(_K);
        _table = new double ****[_N];
        _minIndexArray = new int ****[_N];
        _leftKArray = new int ****[_N];
        for (int i = 0; i < _N; i++) {
            _table[i] = new double ***[_N];
            _minIndexArray[i] = new int ***[_N];
            _leftKArray[i] = new int ***[_N];
            for (int j = 0; j < _N; j++) {
                _table[i][j] = new double **[_K];
                _minIndexArray[i][j] = new int **[_K];
                _leftKArray[i][j] = new int **[_K];
                for (int k = 0; k < _K; k++) {
                    _table[i][j][k] = new double *[_H];
                    _minIndexArray[i][j][k] = new int *[_H];
                    _leftKArray[i][j][k] = new int *[_H];
                    for (int h = 0; h < _H; h++) {
                        _table[i][j][k][h] = new double[_N]{};
                        _minIndexArray[i][j][k][h] = new int[_N]{};
                        _leftKArray[i][j][k][h] = new int[_N]{};
                    }
                }
            }
        }
    }


    Detector::~Detector ()
    {
        for (int i = 0; i < _N; i++) {
            for (int j = 0; j < _N; j++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        delete _table[i][j][k][h];
                        delete _minIndexArray[i][j][k][h];
                        delete _leftKArray[i][j][k][h];
                    }
                    delete _table[i][j][k];
                    delete _minIndexArray[i][j][k];
                    delete _leftKArray[i][j][k];
                }
                delete _table[i][j];
                delete _minIndexArray[i][j];
                delete _leftKArray[i][j];
            }
            delete _table[i];
            delete _minIndexArray[i];
            delete _leftKArray[i];
        }
        delete _table;
        delete _minIndexArray;
        delete _leftKArray;
    }

    void Detector::execute ()
    {
        std::clock_t t = std::clock();
        fillTable();
        t = std::clock() - t;
        std::cout << "fillTable took: " << (float)t/CLOCKS_PER_SEC << "s\n";

        std::vector<utils::intDoublePair> sumOfEntropy;
        std::vector<double> sumOfLeaves;
        std::vector<utils::intDoublePair> normLeaves;

//        std::cout << "_K=" << _K << std::endl;
        int index = -1;
        if (_DETERMINE_K) {
            // determine K
            for (int num = 2; num < _K + 1; num++) {
                printf("--------\nk=%d\n", num);
                sumOfEntropy.emplace_back(num, _table[0][_N - 1][indexK(num)][_H - 1][_N - 1]);

                backTrace(num, _H);

                double leafSum = 0;
                for (int leaf = 0; leaf < _boundaries.size(); leaf++) {
                    int currentStart = _boundaries[leaf].first;
                    int currentEnd = _boundaries[leaf].second;
                    leafSum += _data->getSE(currentStart, currentEnd, 2 * _data->_edgeSum);
                    std::cout << "currentStart=" << currentStart << ", currentEnd=" << currentEnd << "\n";
                    leafSum += _table[currentStart][currentEnd][0][0][currentEnd];
                }
                sumOfLeaves.emplace_back(leafSum);
                double divisor =
                    log2(_N / (double) num) + (_N * (num - 1) / (double) (num * (_N - 1))) * log2((double) num);
                normLeaves.emplace_back(num, leafSum / divisor);
            }
            printf("sumOfEntropy:\n");
            for (int i=0; i<sumOfEntropy.size(); i++)
                printf("(%d, %f)\n", sumOfEntropy[i].first, sumOfEntropy[i].second);

            printf("normLeaves:\n");
            for (int i=0; i<normLeaves.size(); i++)
                printf("(%d, %f)\n", normLeaves[i].first, normLeaves[i].second);

            if (_H == 1 || _H == 2) {
                sort(sumOfEntropy.begin(), sumOfEntropy.end(), utils::cmpIntDoublePairBySecond);
                index = sumOfEntropy[0].first;
            }
            else {
                sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePairBySecond);
                index = normLeaves[0].first;
            }
            std::cout << "k chosen=" << index << std::endl;
        } else {
            int num = _K;
            printf("--------\nk=%d\n", num);
            sumOfEntropy.emplace_back(num, _table[0][_N - 1][indexK(num)][_H - 1][_N - 1]);
            backTrace(num, _H);
            double leafSum = 0;
            for (int leaf = 0; leaf < _boundaries.size(); leaf++) {
                int currentStart = _boundaries[leaf].first;
                int currentEnd = _boundaries[leaf].second;
                leafSum += _data->getSE(currentStart, currentEnd, 2 * _data->_edgeSum);
                std::cout << "currentStart=" << currentStart << ", currentEnd=" << currentEnd << "\n";
                leafSum += _table[currentStart][currentEnd][0][0][currentEnd];
            }
            sumOfLeaves.emplace_back(leafSum);
            double divisor =
                log2(_N / (double) num) + (_N * (num - 1) / (double) (num * (_N - 1))) * log2((double) num);
            normLeaves.emplace_back(num, leafSum / divisor);
            index = _K;
        }

        backTrace(index, _H, true);

        _nodeList = &_multiTree.nodeList();
        for (int i = 0; i < _nodeList->size(); i++) {
            std::cout << (*_nodeList)[i]->_val[0] << ", " <<  (*_nodeList)[i]->_val[1] << std::endl;
        }
        _writer.writeTree(_INPUT+".original_boundaries.txt", *_nodeList);
    }


    void Detector::fillTable ()
    {
        std::clock_t t = std::clock();
        for (int start = 0; start < _N; start++) {
            for (int end = start; end < _N; end++) {
                double currentVolume = _data->getVol(start, end);
                double binSum;
                if (start == 0) {
                    binSum = _data->getGtimesLogG(currentVolume) - _data->_sumOfGtimesLogG[end];
                } else {
                    binSum = _data->getGtimesLogG(currentVolume) - (_data->_sumOfGtimesLogG[end] - _data->_sumOfGtimesLogG[start-1]);
                }
                _table[start][end][0][0][end] = binSum / (2. * _data->_edgeSum);
            }
        }

        t = std::clock() - t;
        std::cout << "part1 took: " << (float)t/CLOCKS_PER_SEC << "s\n";

        t = std::clock();
//        std::cout << "_K=" << _K << std::endl;
        for (int cluster = 1; cluster < _K; cluster++) {
            for (int start = 0; start < _N; start++) {
                for (int parentEnd = start; parentEnd < _N; parentEnd++) {
                    for (int end = start; end < parentEnd + 1; end++) {
                        double parentVol = _data->getVol(start, parentEnd);
                        double minTmp = std::numeric_limits<double>::infinity();
                        int minIdx = 0;
                        for (int i = start; i < end; i++) {
                            double tmp;
                            if (cluster - 1 == 0) {
                                tmp = _table[start][i][cluster - 1][0][i];
                                tmp += _data->getSE(start, i, parentVol);
                            }
                            else {
                                tmp = _table[start][i][cluster - 1][0][parentEnd];
                            }

                            tmp += _data->getSE(i + 1, end, parentVol);

                            tmp += _table[i + 1][end][0][0][end];
                            if (tmp <= minTmp) {
                                minTmp = tmp;
                                minIdx = i;
                            }
                        }
                        _minIndexArray[start][end][cluster][0][parentEnd] = minIdx;
                        _table[start][end][cluster][0][parentEnd] = minTmp;
                    }
                }
            }
        }
        for (int cluster = 0; cluster < _K; cluster++) {
            std::cout << _table[0][_N-1][cluster][0][_N-1] << " ";
        }
        std::cout << "\n";
        t = std::clock() - t;
        std::cout << "part2 spent: " << (float)t/CLOCKS_PER_SEC << "s\n";

        t = std::clock();
        for (int height = 1; height < _H; height++) {
            initH(height);
            for (int cluster = 2; cluster < _K + 1; cluster++) {
                for (int start = 0; start < _N; start++) {
                    for (int parentEnd = start; parentEnd < _N; parentEnd++) {
                        for (int end = start; end < parentEnd + 1; end++) {
                            double minTmp, currentVol, parentVol;
                            int minIdx, leftK;
                            if (end - start + 1 >= cluster) {
                                minTmp = _table[start][end][indexK(cluster)][height - 1][end];
                                parentVol = _data->getVol(start, parentEnd);
                                minTmp += _data->getSE(start, end, parentVol);
                                minIdx = start;
                                leftK = 0;
                                for (int binaryK = 1; binaryK < cluster; binaryK++) {
                                    for (int mid = start; mid < end; mid++) {
                                        double tmp;
                                        if (binaryK == 1) {
                                            tmp = _table[start][mid][indexK(binaryK)][height - 1][mid] +
                                                  _table[mid + 1][end][indexK(cluster - binaryK)][height - 1][end];
                                            tmp += _data->getSE(start, mid, parentVol);
                                        }
                                        else {
                                            tmp = _table[start][mid][indexK(binaryK)][height][parentEnd] +
                                                  _table[mid + 1][end][indexK(cluster - binaryK)][height - 1][end];
                                        }

                                        tmp += _data->getSE(mid + 1, end, parentVol);
                                        if (tmp <= minTmp) {
                                            minTmp = tmp;
                                            minIdx = mid;
                                            leftK = binaryK;
                                        }
                                    }
                                }
                            }
                            else {
                                minTmp = std::numeric_limits<double>::infinity();
                                minIdx = 0;
                                leftK = 0;
                            }
                            _minIndexArray[start][end][indexK(cluster)][height][parentEnd] = minIdx;
                            _table[start][end][indexK(cluster)][height][parentEnd] = minTmp;
                            _leftKArray[start][end][indexK(cluster)][height][parentEnd] = leftK;
                        }
                    }
                }
            }
        }
        t = std::clock() - t;
        std::cout << "part3 took: " << (float)t/CLOCKS_PER_SEC << "s\n";
    }


    void Detector::initH(int h)
    {
        for (int i = 0; i < _N; i++) {
            for (int j = 0; j < _N; j++) {
                for (int k = 0; k < _N; k++) {
                    _table[i][j][0][h][k] = std::numeric_limits<double>::infinity();
                }
            }
        }
    }


    void Detector::backTrace(int k, int h, bool add)
    {
        initK();

        multiSplit(0, _N-1, k, h-1, _N-1, add);
        _boundaries.emplace_back(0, 0);
        sort(_boundaries.begin(), _boundaries.end(), utils::cmpBoundary);

        for (int i = 0; i < _boundaries.size(); i++) {
            if (i == _boundaries.size() - 1) {
                _boundaries[i].second = _N - 1;
            }
            else {
                _boundaries[i].second = _boundaries[i + 1].first - 1;
            }
        }

        for (int i = 0; i < _boundaries.size(); i++) {
            std::cout << "boundary[" << i << "]=(" << _boundaries[i].first << ", " << _boundaries[i].second << ")\n";
        }
    }


    void Detector::initK()
    {
        _boundaries.clear();
    }


    void Detector::multiSplit (int start, int end, int k, int h, int parentEnd, bool add)
    {
        if (k == 1) {
            if (add)
                _multiTree.insert(start, end);
            printf("filling-------------%d %d %d %d k=1, parentEnd: %d %f\n", start, end, k, h, parentEnd, _table[start][end][indexK(k)][h][parentEnd]);
            return;
        }
        else {
            if (h != 0) {
//        std::cout << "1. start=" << start << ", end=" << end << ", k=" << k << ", indexK(k)=" << indexK (k) << ", h=" << h << ", parentEnd=" << parentEnd << std::endl;
//        std::cout << "_leftKArray[start][end][indexK(k)][h][parentEnd]=" << _leftKArray[start][end][indexK(k)][h][parentEnd] << "\n\n";
                int leftK = _leftKArray[start][end][indexK(k)][h][parentEnd];
                if (leftK == 0) {
                    if (add)
                        _multiTree.insert(start, end);
                    printf("filling-------------%d %d %d %d leftK=1, parentEnd: %d %f\n", start, end, k, h, parentEnd, _table[start][end][indexK(k)][h][parentEnd]);
                    multiSplit(start, end, k, h-1, end, add);
                }
                else {
                    int midPos = _minIndexArray[start][end][indexK(k)][h][parentEnd];
//          std::cout << "2. start=" << start << ", end=" << end << ", k=" << k << ", indexK(k)=" << indexK(k) << ", h=" << h << ", parentEnd=" << parentEnd << ", midPos=" << midPos << std::endl;
//          std::cout << "_minIndexArray[start][end][indexK(k)][h][parentEnd]=" << _minIndexArray[start][end][indexK(k)][h][parentEnd] << "\n\n";
                    _boundaries.emplace_back(midPos + 1, 0);
                    if (add)
                        _multiTree.insert(midPos+1, end);
                    printf("filling-------------%d %d %d %d h!=1, parentEnd: %d %f\n", start, end, k, h, parentEnd, _table[start][end][indexK(k)][h][parentEnd]);
                    multiSplit(start, midPos, leftK, h, parentEnd, add);
                    multiSplit(midPos + 1, end, k-leftK, h-1, end, add);
                }
            }
            else {
//        std::cout << "3. start=" << start << ", end=" << end << ", k=" << k << ", indexK(k)=" << indexK(k) << ", h=" << h << ", parentEnd=" << parentEnd << std::endl;
//        std::cout << "_minIndexArray[start][end][indexK(k)][h][parentEnd]=" << _minIndexArray[start][end][indexK(k)][h][parentEnd] << "\n\n";
                int midPos = _minIndexArray[start][end][indexK(k)][h][parentEnd];
                _boundaries.emplace_back(midPos + 1, 0);
                if (add)
                    _multiTree.insert(midPos+1, end);
                printf("filling-------------%d %d %d %d h=1, parentEnd: %d %f\n", start, end, k, h, parentEnd, _table[start][end][indexK(k)][h][parentEnd]);
                multiSplit(start, midPos, k-1, h, parentEnd, add);
            }
        }
    }

}