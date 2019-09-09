//
// Created by wang mengbo on 2019-09-04.
//

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

int n = 10;
int k = 3;

void print2DIntArray (int **array, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << array[i][j] << " ";
    }
    std::cout << "\n";
  }
}

void print2DDoubleArray (double **array, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << array[i][j] << " ";
    }
    std::cout << "\n";
  }
}


bool intArrayEqual (int a1[], int a2[], int n)
{
  for (int i = 0; i < n; i++) {
    if (a1[i] != a2[i])
      return false;
  }
  return true;
}


void modifyIntArray (int a[], int n)
{
  a[0] = -999;
}


int main ()
{

}


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
