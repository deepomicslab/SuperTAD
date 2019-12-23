//
// Created by MSI-PC on 2019/12/14.
//

#ifndef PROGRAM_COMPARE_H
#define PROGRAM_COMPARE_H

#include<cstdio>
#include<queue>
#include<cstring>
#include<vector>
#include<iostream>

class Compare {
private:
    int _n, _m;
    std::vector<std::vector<int>> capacity;
    std::vector<std::vector<int>> adj;

public:
    int bfs(int s, int t, std::vector<int>& parent);

    int maxflow(int s, int t);
};
#endif //PROGRAM_COMPARE_H
