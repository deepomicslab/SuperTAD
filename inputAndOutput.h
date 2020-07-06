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
#include "utils.h"
#include <map>
#include "data.h"


std::string concatePath(std::string path1, std::string path2);

bool pathExist(const std::string &s);

class Reader {
private:
//    std::ifstream _infile;
//    std::string _filePath;
//    std::vector<std::string> _fileNames;

public:
//    Reader(std::string fileName);

    Reader() {};

    ~Reader() {};

    static int parseMatrix(Eigen::MatrixXd &contactMat, std::string filePath);

    static int parseShort(Eigen::MatrixXd &contactMat, std::string filePath);

//    int parseTree(Eigen::MatrixXd &contactMat, std::vector<std::string> &fileNames);
};


class Writer {
private:
//    std::ofstream _outfile;

public:
    Writer() {};

    ~Writer() {};

    static void writeTree(std::string filePath, std::vector<binary::TreeNode *> &nodeList);

    static void writeTree(std::string filePath, std::vector<multi::TreeNode *> &nodeList);

    template<class T>
    static void writeTree(std::string filePath, std::vector<T *> &nodeList);

    template<class T>
    static void writeTreeAsBinList(std::string filePath, std::vector<T *> &nodeList);

    template<class T>
    static void writeTreeAsBedpe(std::string filePath, std::vector<T *> &nodeList);

    static void writeBoundaries(std::string filePath, std::vector<utils::boundary> &boundaryList);

    static void dumpMatrix(Eigen::MatrixXd &mat, std::string outPath);

    static void dumpCoordinates(i2dMap &map, std::string outPath, std::ofstream *f=NULL);

    static void writeListOfCoordinates(str_2_i2dMap &map, std::string outPath);

    static void dumpListOfCoordinates(str_2_i2dMap &map, std::string outPath);
};

#endif //PROGRAM_INPUTANDOUTPUT_H
