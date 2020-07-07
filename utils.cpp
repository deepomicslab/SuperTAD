//
// Created by wang mengbo on 2019-09-02.
//

#include "utils.h"

namespace utils {

    // including low but excluding high
    int randInt (int low, int high) {
        double d = (double)rand() / (double)RAND_MAX * (high - low);
        double intpart;
        if (modf(d, &intpart) > 0.5) {
            int r = ceil (d);
            if (r >= high)
                return high - 1;
            return r;
        }
        else
            return floor(d);
    }


    double randDouble(double low, double high) {
        return (double)rand() / (double)RAND_MAX * (high - low);
    }


    bool equalDoubleArrays(double *a1, double *a2, int n)
    {
        for (int i = 0; i < n; i++) {
            if (a1[i] != a2[i])
                return false;
        }
        return true;
    }


    void copyDoubleArray (double from[], double to[], int n)
    {
        for (int i = 0; i < n; i++) {
            to[i] = from[i];
        }
    }

    template<typename T>
    void print3DArray(T ***array, int n, int m, int w)
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


    std::string license()
    {
        return "SuperTAD v1.1";
    }

}