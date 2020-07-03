//
// Created by wang mengbo on 2019-09-04.
//

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <ctime>
#include "utils.h"


int main() {
    srand(time(0));
    int n = 10;
    int count = 0;
    int r;
//    double ***array = new double **[n];
//    for (int i=0; i<n; i++) {
//        r = rand() % 10 + 1;
//        array[i] = new double *[r];
//        for (int j=0; j<r; j++) {
//            r = (int)rand()%10+1;
//            array[i][j] = new double [r];
//            count += r;
//        }
//    }
//    double **array = new double *[n];
    int sizes[n];
    for (int i=0; i<n; i++) {
        r = rand() % 10 + 1;
        printf("r=%d\n", r);
//        array[i] = new double [r]{};
        sizes[i] = r;
        count += r;
    }
    fflush(stdout);
    double **array = (double**) malloc(count*sizeof(double));
    memset(array, -1, sizeof(array));
    for (int i=0; i<n; i++) {
        for (int j=0; j<sizes[i]; j++) {
            printf("%f ", array[i][j]);
        }
        printf("\n");
    }
}


//int *_n = new int(0);
//
//int indexK1(int k) {
//    return k-1;
//}
//
//void indexK2(int k) {
//    *_n = k-1;
//}
//
//float testIndex1(int n)
//{
//    std::clock_t t = std::clock();
//    for (int i = 0; i < n; i++) {
//        indexK1(i);
//    }
//    t = std::clock() - t;
//    return (float)t/CLOCKS_PER_SEC;
//}
//
//float testIndex2(int n)
//{
//    std::clock_t t = std::clock();
//    for (int i = 0; i < n; i++) {
//        indexK2(i);
//    }
//    return (float) (std::clock() - t)/CLOCKS_PER_SEC;
//}
//
//float testIndex3(int n)
//{
//    std::clock_t t = std::clock();
//    for (int i = 0; i < n; i++) {
//        i-1;
//    }
//    return (float) (std::clock() - t)/CLOCKS_PER_SEC;
//}
//
//void testIndex()
//{
//    int n = 1e+9;
//    float t1, t2, t3;
//    for (int i=0; i<5; i++) {
//        t1 = testIndex1(n);
//        t2 = testIndex2(n);
//        t3 = testIndex3(n);
//        printf("t1=%f, t2=%f, t3=%f\n", t1, t2, t3);
//        std::fflush(stdout);
//    }
//    //    std::clock_t t1 = 0;
////    std::clock_t t2 = 0;
////    std::clock_t t3 = 0;
////    std::clock_t tmp;
////    for (int k=0; k<10; k++) {
////        t1 = 0;
////        t2 = 0;
////        t3 = 0;
////        for (int i = 0; i < 1e+7; i++) {
////            tmp = std::clock();
////            indexK1(i);
////            t1 += std::clock() - tmp;
////
////            tmp = std::clock();
////            indexK2(i);
////            t2 += std::clock() - tmp;
////
////            tmp = std::clock();
////            i - 1;
////            t3 += std::clock() - tmp;
////        }
////        printf("t1=%f, t2=%f, t3=%f\n", (float) t1 / CLOCKS_PER_SEC, (float) t2 / CLOCKS_PER_SEC,
////               (float) t3 / CLOCKS_PER_SEC);
////        std::fflush(stdout);
////    }
//}
//
//
//int main()
//{
//    testIndex();
//}


//int n = 10;
//int k = 3;
//
//void print2DIntArray (int **array, int m, int n) {
//  for (int i = 0; i < m; i++) {
//    for (int j = 0; j < n; j++) {
//      std::cout << array[i][j] << " ";
//    }
//    std::cout << "\n";
//  }
//}
//
//void print2DDoubleArray (double **array, int m, int n) {
//  for (int i = 0; i < m; i++) {
//    for (int j = 0; j < n; j++) {
//      std::cout << array[i][j] << " ";
//    }
//    std::cout << "\n";
//  }
//}
//
//
//bool intArrayEqual (int a1[], int a2[], int n)
//{
//  for (int i = 0; i < n; i++) {
//    if (a1[i] != a2[i])
//      return false;
//  }
//  return true;
//}
//
//
//void modifyIntArray (int a[], int n)
//{
//  a[0] = -999;
//}
//int main()
//{
//    std::vector<utils::intDoublePair> t;
//    for (int i=0; i<10; i++) {
//        t.emplace_back(rand()%10, (double)rand()/RAND_MAX);
//    }
//    std::sort(t.begin(), t.end(), utils::cmpIntDoublePairBySecond);
//    for (int i=0; i<10; i++) {
//        std::cout << t[i].first << ", " << t[i].second << std::endl;
//    }
//}


//int main ()
//{
//    int n = 10;
//
//    Eigen::MatrixXd M1(n,n);
//    double **M2 = new double *[n];
//    for (int i=0; i<n; i++) {
//        M2[i] = new double [n];
//        for (int j=0; j<n; j++) {
//            double tmp = (double)rand() / RAND_MAX;
//            M1(i, j) = tmp;
//            M2[i][j] = tmp;
//        }
//    }
//    std::cout << "M1:\n" << M1 << "\n";
//    std::cout << "M2:\n";
//    print2DDoubleArray(M2, n, n);
//
//    int times = 1e+5;
//    std::clock_t t1=0;
//    std::clock_t t2=0;
//    std::clock_t tmpT;
//    for (int i=0; i<times; i++) {
//        int x = rand() % 10;
//        int y = rand() % 10;
//        tmpT = std::clock();
//        M1.coeff(x, y);
//        t1 += std::clock() - tmpT;
//
//        tmpT = std::clock();
//        M2[x][y];
//        t2 += std::clock() - tmpT;
//    }
//    std::cout << "M1 running time: " << (float)t1/CLOCKS_PER_SEC << "s\n";
//    std::cout << "M2 running time: " << (float)t2/CLOCKS_PER_SEC << "s\n";
//}


//int main () {
//  int t1[2] {};
//  int t2[2] {0, 1};
//  if (intArrayEqual (t1, t2, 2))
//    std::cout << "true" << std::endl;
//  else
//    std::cout << "false" << std::endl;
//
//  modifyIntArray (t1, 2);
//  std::cout << t1[0] << " " << t1[1] << std::endl;
//
//}


//int main () {
//  int _N = 70;
//  int num = 11;
//  std::cout << log2((double)_N / (double)num) << std::endl;
//  std::cout << log2(_N / (double)num) + (_N * (num - 1) / (double)(num * (_N - 1))) * log2(num) << std::endl;
//}

//int main (int argc, char *argv[]) {
//  std::ofstream out;
//  std::string path = "./testFileOutput.txt";
//  std::cout << "path=" << path << "\n";
//  out.open(path);
//  out << "ok\n";
//  out.close();
//  return 0;
//}

//int main (int argc, char *argv[])
//{
//  double **table = new double *[n];
//  for (int i = 0; i < n; i++) {
//    table[i] = new double [k] {};
//  }
//
//  print2DDoubleArray (table, n, k);
//
//  return 0;
//}
