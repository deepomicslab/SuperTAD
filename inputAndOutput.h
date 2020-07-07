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

//    static void writeTree(std::string filePath, std::vector<binary::TreeNode *> &nodeList);

//    static void writeTree(std::string filePath, std::vector<multi::TreeNode *> &nodeList);

    template<class T>
    static void writeTreeAsBinList(std::string filePath, std::vector<T *> &nodeList)
    {
        std::ofstream outFile;
        outFile.open(filePath);
        if (outFile.is_open()) {
            if (_VERBOSE_)
                printf("start writing tree into %s\n", filePath.c_str());

            for (int i = 0; i < nodeList.size(); i++) {
                for (int j = nodeList[i]->_val[0]; j <= nodeList[i]->_val[1]; j++)
                    outFile << std::to_string(j + 1) << " ";
                outFile << "\n";
            }
            outFile.close();

            if (_VERBOSE_)
                std::cout << "finish writing tree\n";
            else
                printf("write tree into: %s\n", filePath.c_str());
        }
        else {
            std::cerr << "cannot open file: " << filePath << "\n";
        }
    }

    template<class T>
    static void writeTreeAsBedpe(std::string filePath, std::vector<T *> &nodeList)
    {
        FILE *outFile = NULL;
        outFile = std::fopen(filePath.c_str(), "w");
        if (outFile) {
            if (_VERBOSE_)
                printf("start writing tree into %s\n", filePath.c_str());

            fprintf(outFile, "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\n");
            fprintf(outFile, "#%s\n", utils::version().c_str());
            int chr1Start, chr1End, chr2Start, chr2End;
            for (int i = 0; i < nodeList.size(); i++) {
                chr1Start = _CHROM1_START_ + (nodeList[i]->_val[0]-1) * _RESOLUTION_;
                chr1End =   _CHROM1_START_ + (nodeList[i]->_val[1]-1) * _RESOLUTION_;
                chr2Start = _CHROM2_START_ + (nodeList[i]->_val[0]-1) * _RESOLUTION_;
                chr2End =   _CHROM2_START_ + (nodeList[i]->_val[1]-1) * _RESOLUTION_;
                fprintf(outFile, "%s\t%d\t%d\t%s\t%d\t%d\tnode%d\n",
                        _CHROM1_.c_str(), chr1Start, chr1End, _CHROM2_.c_str(), chr2Start, chr2End, i+1);
            }
            fclose(outFile);

            if (_VERBOSE_)
                printf("finish writing tree\n");
            else
                printf("write tree into: %s\n", filePath.c_str());
        }
        else
            fprintf(stderr, "cannot open file: %s\n", filePath);
    }

    template<class T>
    static void writeTree(std::string filePath, std::vector<T *> &nodeList)
    {
        if (_BEDPE_)
            writeTreeAsBedpe(filePath+".bedpe", nodeList);
        else
            writeTreeAsBinList(filePath+".txt", nodeList);
    }

    static void writeBoundaries(std::string filePath, std::vector<boundary> &boundaryList);

    static void dumpMatrix(Eigen::MatrixXd &mat, std::string outPath);

    static void dumpCoordinates(i2dMap &map, std::string outPath, std::ofstream *f=NULL);

    static void writeListOfCoordinates(str_2_i2dMap &map, std::string outPath);

    static void dumpListOfCoordinates(str_2_i2dMap &map, std::string outPath);
};

#endif //PROGRAM_INPUTANDOUTPUT_H
