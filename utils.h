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
#include <limits>


typedef std::map<int, double> Int2DoubleMap;
typedef std::map<std::string, std::map<int, double>> Str_2_Int2DoubleMap;
//typedef std::pair<int, int> Boundary;
typedef std::pair<int, double> IntDoublePair;


struct Boundary{
    int first, second, size;
    Boundary() {};
    Boundary(int s, int e) {
        first = s;
        second = e;
        size = e - s + 1;
    };
};

namespace utils {

    int boundariesIntersection(Boundary &b1, Boundary &b2);

    int randInt(int low=0, int high=10);

    double randDouble(double low=0., double high=1.);

    inline bool cmpBoundary(const Boundary &p1, const Boundary &p2) { return p1.first < p2.first; }

    inline bool cmpIntDoublePairBySecond(const IntDoublePair &p1, const IntDoublePair &p2) { return p1.second < p2.second; }

    bool equalDoubleArrays(double *a1, double *a2, int n);

    void copyDoubleArray(double from[], double to[], int n);

    template<typename T>
    void print3Darray(T ***array, int n, int m, int w)
    {
        for (int i = 0; i < n; i++) {
            printf("----i=%d----\n", i);
            for (int j=0; j<m; j++) {
                for (int k=0; k<w; k++) {
                    std::cout << array[i][j][k] << " ";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
    }

//    template<typename T, U>
//    void printPairVec(std::vector<std::pair<T, U> ) {
//        for (v=V.begin(); v!=V.end(); v++) {
//            std::cout << "(" << v.first << ", " << v.second << ")"
//        }
//    }

    template<typename T>
    void print2Darray(T **array, int n, int m)
    {
        for (int i = 0; i < n; i++) {
            for (int j=0; j < m; j++) {
//                printf("%s\t", std::to_string(array[i][j]).c_str());
                if (array[i][j] >= std::numeric_limits<int>::max())
                    std::cout << "inf\t";
                else
                std::cout << array[i][j] << "\t";
            }
            printf("\n");
        }
        printf("\n");
    }

    std::string version();

}

#endif //PROGRAM_UTILS_H
