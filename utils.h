//
// Created by wang mengbo on 2019-09-02.
//

#ifndef PROGRAM_UTILS_H
#define PROGRAM_UTILS_H

#include <cmath>
#include <utility>
#include <stdlib.h>
#include <iostream>
#include <map>
#include <sstream>


typedef std::map<int, double> i2dMap;
typedef std::map<std::string, std::map<int, double>> str_2_i2dMap;
typedef std::pair<int, int> boundary;
typedef std::pair<int, double> intDoublePair;

namespace utils {

    int randInt(int low=0, int high=10);

    double randDouble(double low=0., double high=1.);

//    typedef std::pair<int, int> boundary;

    inline bool cmpBoundary(const boundary &p1, const boundary &p2) { return p1.first < p2.first; }

//    typedef std::pair<int, double> intDoublePair;

    inline bool cmpIntDoublePairBySecond(const intDoublePair &p1, const intDoublePair &p2) { return p1.second < p2.second; }

    bool equalDoubleArrays(double *a1, double *a2, int n);

    void copyDoubleArray(double from[], double to[], int n);

    template<typename T>
    void print3DArray(T ***array, int n, int m, int w);

//    template<typename T, U>
//    void printPairVec(std::vector<std::pair<T, U> ) {
//        for (v=V.begin(); v!=V.end(); v++) {
//            std::cout << "(" << v.first << ", " << v.second << ")"
//        }
//    }

    template<typename T>
    static void print2Dtable(T **table);

    std::string license();

}

#endif //PROGRAM_UTILS_H
