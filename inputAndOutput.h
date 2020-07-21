//
// Created by wang mengbo on 2019-09-02.
//

#ifndef PROGRAM_INPUTANDOUTPUT_H
#define PROGRAM_INPUTANDOUTPUT_H

#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include "binaryTree.h"
#include "multiTree.h"
#include "utils.h"
#include <map>
#include "data.h"


bool pathExist(const std::string &s);

class Reader {
private:
public:
    Reader() {};

    ~Reader() {};

    static void parseMatrix2Array(double **&table, std::string path);
};


class Writer {
public:
    Writer() {};

    ~Writer() {};

    template<class T>
    static void writeTreeAsBinList(std::string filePath, std::vector<T *> &nodeList)
    {
        filePath += ".txt";
        std::ofstream outFile;
        outFile.open(filePath);
        if (outFile.is_open()) {
            if (_VERBOSE_)
                printf("start writing tree into %s\n", filePath.c_str());
            else
                printf("write tree into %s\n", filePath.c_str());

            for (int i = 0; i < nodeList.size(); i++) {
                for (int j = nodeList[i]->_val[0]; j <= nodeList[i]->_val[1]; j++)
                    outFile << std::to_string(j + 1) << " ";
                outFile << "\n";
            }
            outFile.close();

            if (_VERBOSE_)
                std::cout << "finish writing tree\n";

        }
        else {
            std::cerr << "cannot open file: " << filePath << "\n";
        }
    }

    template<class T>
    static void writeTreeIn7Cols(std::string filePath, std::vector<T *> &nodeList)
    {
        filePath += ".tsv";
        FILE *outFile = NULL;
        outFile = std::fopen(filePath.c_str(), "w");
        if (outFile) {
            if (_VERBOSE_)
                printf("start writing tree into %s\n", filePath.c_str());
            else
                printf("write tree into %s\n", filePath.c_str());

            int bin1Idx, bin1Start, bin1End, bin2Idx, bin2Start, bin2End;
            for (int i = 0; i < nodeList.size(); i++) {
                bin1Idx = nodeList[i]->_val[0];
                bin1Start = _CHROM1_START_ + bin1Idx * _RESOLUTION_;
                bin1End = _CHROM1_START_ + (bin1Idx+1) * _RESOLUTION_;
                bin2Idx = nodeList[i]->_val[1];
                bin2Start = _CHROM1_START_ + bin2Idx * _RESOLUTION_;
                bin2End = _CHROM1_START_ + (bin2Idx+1) * _RESOLUTION_;
                fprintf(outFile, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
                        _CHROM1_.c_str(), bin1Idx+1, bin1Start, bin1End, bin2Idx+1, bin2Start, bin2End);
            }
            fclose(outFile);

            if (_VERBOSE_)
                printf("finish writing tree\n");

        }
        else
            fprintf(stderr, "cannot open file: %s\n", filePath);
    }

    template<class T>
    static void writeTreeAsBedpe(std::string filePath, std::vector<T *> &nodeList)
    {
        filePath += ".bedpe";
        FILE *outFile = NULL;
        outFile = std::fopen(filePath.c_str(), "w");
        if (outFile) {
            if (_VERBOSE_)
                printf("start writing tree into %s\n", filePath.c_str());
            else
                printf("write tree into %s\n", filePath.c_str());

            fprintf(outFile, "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\n");
            fprintf(outFile, "#%s\n", utils::version().c_str());

//            if (_RESOLUTION_==1) {
//                _CHROM1_START_ = 1;
//                _CHROM2_START_ = 1;
//            }

            int bin1ChrIdx, bin1Start, bin1End, bin2ChrIdx, bin2Start, bin2End;

            for (int i = 0; i < nodeList.size(); i++) {
                bin1ChrIdx = nodeList[i]->_val[0];
                bin2ChrIdx = nodeList[i]->_val[1];

                bin1Start = _CHROM1_START_ + bin1ChrIdx * _RESOLUTION_;
                bin1End = _CHROM1_START_ + bin2ChrIdx * _RESOLUTION_;

                bin2Start = _CHROM2_START_ + bin1ChrIdx * _RESOLUTION_;
                bin2End = _CHROM2_START_ + bin2ChrIdx * _RESOLUTION_;

                fprintf(outFile, "%s\t%d\t%d\t%s\t%d\t%d\tnode%d\n",
                        _CHROM1_.c_str(), bin1Start, bin1End, _CHROM2_.c_str(), bin2Start, bin2End, i + 1);
            }
            fclose(outFile);

            if (_VERBOSE_)
                printf("finish writing tree\n");

        }
        else
            fprintf(stderr, "cannot open file: %s\n", filePath);
    }

    template<class T>
    static void writeTreeInShort(std::string filePath, std::vector<T *> &nodeList)
    {
        filePath += ".short.txt";
        FILE *outFile = NULL;
        outFile = std::fopen(filePath.c_str(), "w");
        if (outFile) {
            if (_VERBOSE_)
                printf("start writing tree into %s\n", filePath.c_str());
            else
                printf("write tree into %s\n", filePath.c_str());

            fprintf(outFile, "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\n");
            fprintf(outFile, "#%s\n", utils::version().c_str());

            int bin1ChrIdx, bin1Start, bin1End, bin2ChrIdx, bin2Start, bin2End;

            for (int i = 0; i < nodeList.size(); i++) {
                bin1Start = _CHROM1_START_ + nodeList[i]->_val[0] * _RESOLUTION_;
                bin2Start = _CHROM2_START_ + nodeList[i]->_val[1] * _RESOLUTION_;

                fprintf(outFile, "0\t%s\t%d\t0\t0\t%s\t%d\t1\t%f\n",
                        _CHROM1_.c_str(), bin1Start, _CHROM2_.c_str(), bin2Start);
            }
            fclose(outFile);

            if (_VERBOSE_)
                printf("finish writing tree\n");

        }
        else
            fprintf(stderr, "cannot open file: %s\n", filePath);
    }

    template<class T>
    static void writeTree(std::string filePath, std::vector<T *> &nodeList)
    {
        if (_BEDPE_)
            writeTreeAsBedpe(filePath, nodeList);
        else if (_SHORT_)
            writeTreeInShort(filePath, nodeList);
        else if (_BIN_LIST_)
            writeTreeAsBinList(filePath, nodeList);
        else {
            if (_CHROM1_ != _CHROM2_) {
                fprintf(stderr, "chromosome indices are not the same, output will be written in BEDPE\n");
                writeTreeAsBedpe(filePath, nodeList);
            }
            writeTreeIn7Cols(filePath, nodeList);
        }
    }

    static void writeBoundaries(std::string path, std::vector<boundary> &boundaryList);

    static void dumpCoordinates(i2dMap &map, std::string path, std::ofstream *f=NULL);

    static void writeListOfCoordinates(str_2_i2dMap &map, std::string outPath);

    static void dumpListOfCoordinates(str_2_i2dMap &map, std::string outPath);
};

#endif //PROGRAM_INPUTANDOUTPUT_H
