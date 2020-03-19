//
// Created by wang mengbo on 2019-09-03.
//

#include "detectorMulti.h"


namespace multi {

    DetectorBase::DetectorBase(Data & data)
    {
        _data = &data;
        _edgeCount = &data.edgeCount();
        int k = 1;
        for (int i = 0; i < _K; i++) {
            _kToIdx.emplace(k++, i);
        }
    }


    DetectorH1::DetectorH1(Data & data) : DetectorBase(data)
    {
        _table = new double * [_N];
        _minIndexArray = new int * [_N];
        _leftKArray = new int * [_N];
//        std::cout << "_K=" << _K << std::endl;
        for (int i = 0; i < _N; i++) {
            _table[i] = new double [_K]{};
            _minIndexArray[i] = new int [_K]{};
            _leftKArray[i] = new int [_K]{};
        }
    }


    void DetectorH1::printTable()
    {
        printf("table:\n");
        for (int i=0; i<_N; i++) {
            for (int j=0; j<_K; j++) {
                printf("%f " ,_table[i][j]);
            }
            printf("\n");
        }
    }


    DetectorH1::~DetectorH1()
    {
        for (int i = 0; i < _N; i++) {
            delete _table[i];
            delete _minIndexArray[i];
            delete _leftKArray[i];
        }
    }


    void DetectorH1::execute()
    {
        printf("execute multi, h=1\n");
        fillTable();
        if (_DETERMINE_K) {
            printf("determine k\n");
            std::vector<utils::intDoublePair> tmp;
            for (int i=0; i<_K; i++)
                tmp.emplace_back(i, _table[_N-1][i]);
            std::sort(tmp.begin(), tmp.end(), utils::cmpIntDoublePair);
            _K = tmp.front().first;
            printf("optimal k=%d\n", _K);
        }
        backTrace();
        Writer::writerBoundaryList(_INPUT + ".original_boundaries.txt", _boundaryList);
    }


    void DetectorH1::initBoundary()
    {
        _boundaryList.clear();
    }


    void DetectorH1::fillTable()
    {
        printf("fill table\n");
        _table[0][0] = _data->getSE(0, 0, 2* _data->getEdgeSum());
        for (int i=1; i<_N; i++) {
            double currentVol = _data->getVol(0, i);
//            printf("currentVol=%f\n", currentVol);
            _table[i][0] = _data->getSE(0, i, 2*_data->getEdgeSum());
//            printf("_table[%d][0]=%f\n", i, _table[i][0]);
            for (int leaf=0; leaf < i + 1; leaf++) {
                _table[i][0] += _data->getSE(leaf, leaf, currentVol);
            }
        }
//        printTable();
        for (int a=1; a<_K; a++) {
            double minTmp;
            int minIdx;
            for (int b=0; b<_N; b++) {
                minTmp = std::numeric_limits<double>::infinity();
                minIdx = 0;
                double tmp;
                for (int i=0; i<b; i++) {
                    tmp = _table[i][a-1];
                    if (i+1==b)
                        tmp += _data->getSE(b, b, 2*_data->getEdgeSum());
                    else {
                        double currentVol = _data->getVol(i+1, b);
                        tmp += _data->getSE(i+1, b, 2*_data->getEdgeSum());
                        for (int leaf=i+1; leaf<b+1; leaf++)
                            tmp += _data->getSE(leaf, leaf, currentVol);
                    }
                    if (tmp<=minTmp) {
                        minTmp = tmp;
                        minIdx = i;
                    }
                }
                _minIndexArray[b][a] = minIdx;
                _table[b][a] = minTmp;
            }
            printf("finish k=%d, structure entropy=%f\n", a, minTmp);
        }
//        printTable();
    }


    void DetectorH1::backTrace()
    {
        initBoundary();
        int lastBoundary = _N;
        int newBoundary;
        for (int i=_K-2; i>=0; i--) {
            newBoundary = _minIndexArray[lastBoundary-1][i+1]+1;
            _boundaryList.emplace_back(newBoundary, lastBoundary);
            lastBoundary = newBoundary;
        }
        _boundaryList.emplace_back(1, lastBoundary);
        std::sort(_boundaryList.begin(), _boundaryList.end(), utils::cmpBoundary);
//        printf("_boundaryList.size()=%d\n", _boundaryList.size());
//        for (std::vector<utils::boundary>::iterator it=_boundaryList.begin(); it!=_boundaryList.end(); it++) {
//            for (int i=it->first; i<=it->second; i++)
//                std::cout << i << " ";
//            std::cout << "\n";
//        }
    }


    Detector::Detector(Data &data)
    {
        _data = &data;
        _edgeCount = &data.edgeCount();

        int k = 1;
        for (int i = 0; i < _K; i++) {
            _kToIdx.emplace(k++, i);
        }

        _table = new double ****[_N];
        _minIndexArray = new int ****[_N];
        _leftKArray = new int ****[_N];
//        std::cout << "_K=" << _K << std::endl;
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
                for (int leaf = 0; leaf < _boundary.size(); leaf++) {
                    int currentStart = _boundary[leaf].first;
                    int currentEnd = _boundary[leaf].second;
                    leafSum += _data->getSE(currentStart, currentEnd, 2 * _data->getEdgeSum());
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
                sort(sumOfEntropy.begin(), sumOfEntropy.end(), utils::cmpIntDoublePair);
                index = sumOfEntropy[0].first;
            }
            else {
                sort(normLeaves.begin(), normLeaves.end(), utils::cmpIntDoublePair);
                index = normLeaves[0].first;
            }
            std::cout << "k chosen=" << index << std::endl;
        } else {
            int num = _K;
            printf("--------\nk=%d\n", num);
            sumOfEntropy.emplace_back(num, _table[0][_N - 1][indexK(num)][_H - 1][_N - 1]);
            backTrace(num, _H);
            double leafSum = 0;
            for (int leaf = 0; leaf < _boundary.size(); leaf++) {
                int currentStart = _boundary[leaf].first;
                int currentEnd = _boundary[leaf].second;
                leafSum += _data->getSE(currentStart, currentEnd, 2 * _data->getEdgeSum());
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
//                  printf("currentVolume=%f\n", currentVolume);
                for (int leaf = start; leaf < end + 1; leaf++) {
//                      printf("leafDegree=%f\n", _data->getVol(leaf, leaf));
                    double SE = _data->getSE(leaf, leaf, currentVolume);
//                      printf("SE=%f\n", SE);
                    _table[start][end][0][0][end] += SE;
                }
            }
        }
//    for (int i = 0; i < _N; i++) {
//      std::cout << _table[0][i][0][0][i] << " ";
//    }
//    std::cout << "\n";
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
                            } else {
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
//          for (int tmpIdx = start; tmpIdx < parentEnd + 1; tmpIdx++) {
//            std::cout << _table[start][tmpIdx][cluster][0][parentEnd] << " ";
//          }
//          std::cout << "\n";
                }
            }
        }
        for (int cluster = 0; cluster < _K; cluster++) {
            std::cout << _table[0][_N-1][cluster][0][_N-1] << " ";
        }
        std::cout << "\n";
        t = std::clock() - t;
        std::cout << "part2 took: " << (float)t/CLOCKS_PER_SEC << "s\n";

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


    void Detector::initH (int h)
    {
        for (int i = 0; i < _N; i++) {
            for (int j = 0; j < _N; j++) {
                for (int k = 0; k < _N; k++) {
                    _table[i][j][0][h][k] = std::numeric_limits<double>::infinity();
                }
            }
        }
    }


    void Detector::backTrace (int k, int h, bool add)
    {
        initK();

        multiSplit(0, _N-1, k, h-1, _N-1, add);
        _boundary.emplace_back(0, 0);
        sort(_boundary.begin(), _boundary.end(), utils::cmpBoundary);

        for (int i = 0; i < _boundary.size(); i++) {
            if (i == _boundary.size() - 1) {
                _boundary[i].second = _N - 1;
            }
            else {
                _boundary[i].second = _boundary[i+1].first - 1;
            }
        }

        for (int i = 0; i < _boundary.size(); i++) {
            std::cout << "boundary[" << i << "]=(" << _boundary[i].first << ", " << _boundary[i].second << ")\n";
        }
    }


    void Detector::initK()
    {
        _boundary.clear();
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
                    _boundary.emplace_back(midPos+1, 0);
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
                _boundary.emplace_back(midPos+1, 0);
                if (add)
                    _multiTree.insert(midPos+1, end);
                printf("filling-------------%d %d %d %d h=1, parentEnd: %d %f\n", start, end, k, h, parentEnd, _table[start][end][indexK(k)][h][parentEnd]);
                multiSplit(start, midPos, k-1, h, parentEnd, add);
            }
        }
    }

}