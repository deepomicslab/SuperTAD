//
// Created by wang mengbo on 2019-09-02.
//

#ifndef PROGRAM_INPUTANDOUTPUT_H
#define PROGRAM_INPUTANDOUTPUT_H

#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <Eigen/Dense>
#include "binaryTree.h"
#include "multiTree.h"


std::string concatePath(std::string path1, std::string path2);

bool isPathExist(const std::string &s);

class Reader {
private:
    std::ifstream _infile;
    std::string _fileName;
    std::vector<std::string> _fileNames;

public:
    Reader(std::string fileName);

    ~Reader() {};

    int parse(Eigen::MatrixXd &contactMat, std::string fileName="");

    int parseTree(Eigen::MatrixXd &contactMat, std::vector<std::string> &fileNames);
};


class Writer {
private:
    std::ofstream _outfile;

public:
    Writer() {};

    ~Writer() {};

    void writeTree(std::string filePath, std::vector<binary::TreeNode *> &nodeList);

    void writeTree(std::string filePath, std::vector<multi::TreeNode *> &nodeList);

    static void dumpMatrix(Eigen::MatrixXd &mat, std::string outpath);
};

#endif //PROGRAM_INPUTANDOUTPUT_H
