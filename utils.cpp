//
// Created by wang mengbo on 2019-09-02.
//

#include "utils.h"

namespace utils {

    int boundariesIntersection(Boundary &b1, Boundary &b2)
    {
        if (b1.first > b2.second || b2.first > b1.second)
            return 0;
        else {
            int s = std::max(b1.first, b2.first);
            int e = std::min(b1.second, b2.second);
//            printf("b1=(%d, %d), b2=(%d, %d)\n", b1.first, b1.second, b2.first, b2.second);
            return e-s+1;
        }
    }


    // including low but excluding high
    int randInt (int low, int high)
    {
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


    double randDouble(double low, double high)
    {
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


    int LCM(int num1, int num2)
    {
        int min, max, temp;
        if (num1 < num2)
        {
            max = num2;
            min = num1;
        }
        else
        {
            max = num1;
            min = num2;
        }
        while (max%min != 0)
        {
            temp = max%min;
            max = min;
            min = temp;
        }

        temp = num1*num2 / temp;

        return temp;
    }


    std::string version()
    {
        return "SuperTAD v1.1";
    }

}