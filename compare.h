//
// Created by wang mengbo on 2019/12/14.
//

#ifndef PROGRAM_COMPARE_H
#define PROGRAM_COMPARE_H

#include<cstdio>
#include<queue>
#include<cstring>
#include<vector>
#include<iostream>
#include <string.h>
#include "inputAndOutput.h"
#include "Hungarian.h"


class Comparator {
private:
    int _n, _n1, _n2;
    int **_graph;
//    int **_rGraph;
    std::vector<Boundary> _boundaries1, _boundaries2;

public:
    Comparator(){}

    Comparator(std::string path1, std::string path2);

    ~Comparator();

    void init(std::string path1, std::string path2);

//    void execute();

    void execute();

//    bool bfs(int **rGraph, int s, int t, int *parent);
//
//    int fordFulkerson(int **graph, int **rGraph, int s, int t);

};
#endif //PROGRAM_COMPARE_H
