//
// Created by wang mengbo on 2019-09-02.
//

#include "inputAndOutput.h"


std::string concatePath(std::string path1, std::string path2)
{
//  if (path1[path1.length ()-11] == '/')
//    path1 = path1[]
    return "";
}


bool pathExist(const std::string &s)
{
    struct stat buffer;
    return (stat(s.c_str(), &buffer) == 0);
}


//Reader::Reader(std::string fileName)
//{
//    _filePath = fileName;
//}


int Reader::parseMatrix(Eigen::MatrixXd &contactMat, std::string filePath)
{
    if (filePath == "") {
        std::cerr << "input must be provided\n";
        exit(1);
    }

    if (!pathExist(filePath)) {
        printf("input file doesn't exist\n");
        exit(1);
    }

    std::ifstream inFile;
    inFile.exceptions(std::ifstream::badbit);
    try {
        inFile.open(filePath);
        int count = 0;
        if (inFile.is_open()) {
            if (_VERBOSE_)
                printf("start parsing input: %s\n", filePath.c_str());
            std::string line;
            double c;
            bool init = false;
            int pos1 = 0;
            while (getline(inFile, line)) {
                std::istringstream iss(line);
                if (!init) {
                    std::string lineBack = line;
                    std::istringstream issBack(lineBack);
                    while (issBack >> c) {
                        count++;
                    }
                    contactMat.resize(count, count);
                    init = true;
                }
                int pos2 = 0;
                while (iss >> c) {
                    contactMat(pos1, pos2) = c;
                    pos2++;
                }
                pos1++;
            }
            inFile.close();
        }

        if (_VERBOSE_)
            printf("finish parsing input\n");
        else
            printf("parse input: %s\n", filePath.c_str());

        return count;
    }
    catch (const std::ifstream::failure& e) {
        printf("exception reading file\n");
        exit(1);
    }

    return 0;
}


//int Reader::parseTree(Eigen::MatrixXd &contactMat, std::vector<std::string> &fileNames)
//{
//    for (int i=0; i<fileNames.size(); i++) {
//        _infile.open(fileNames[i]);
//
//        _infile.close();
//    }
//}


void Writer::writeTree(std::string filePath, std::vector<binary::TreeNode *> &nodeList)
{
    std::ofstream outFile;
    outFile.open(filePath);
    if (outFile.is_open()) {
        if (_VERBOSE_)
            printf("start writing binary tree into: %s\n", filePath.c_str());

        for (int i = 0; i < nodeList.size(); i++) {
            for (int j = nodeList[i]->_val[0]; j <= nodeList[i]->_val[1]; j++)
                outFile << std::to_string(j + 1) << " ";
            outFile << "\n";
        }
        outFile.close();

        if (_VERBOSE_)
            std::cout << "finish writing binary tree\n";
        else
            printf("write binary tree into: %s\n", filePath.c_str());
    }
    else
        std::cerr << "cannot open file: " << filePath << "\n";
}


void Writer::writeTree(std::string filePath, std::vector<multi::TreeNode *> &nodeList)
{
    std::ofstream outFile;
    outFile.open(filePath);
    if (outFile.is_open()) {
        if (_VERBOSE_)
            printf("start writing multi-nary tree into %s\n", filePath.c_str());

        for (int i = 0; i < nodeList.size(); i++) {
            for (int j = nodeList[i]->_val[0]; j <= nodeList[i]->_val[1]; j++)
                outFile << std::to_string(j + 1) << " ";
            outFile << "\n";
        }
        outFile.close();

        if (_VERBOSE_)
            std::cout << "finish writing multi-nary tree\n";
        else
            printf("write tree into: %s\n", filePath.c_str());
    }
    else
        std::cerr << "cannot open file: " << filePath << "\n";
}


template<class T>
void Writer::writeTree(std::string filePath, std::vector<T *> &nodeList)
{
    if (_BEDPE_)
        writeTreeAsBedpe(filePath+".bedpe", nodeList);
    else
        writeTreeAsBinList(filePath+".txt", nodeList);
}


template<class T>
void Writer::writeTreeAsBinList(std::string filePath, std::vector<T *> &nodeList)
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
void Writer::writeTreeAsBedpe(std::string filePath, std::vector<T *> &nodeList, Data *data)
{
    FILE *outFile = NULL;
    outFile = std::fopen(filePath.c_str(), "w");
    if (outFile) {
        if (_VERBOSE_)
            printf("start writing tree into %s\n", filePath.c_str());

        fprintf(outFile, "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\n");

        int chr1Start, chr1End, chr2Start, chr2End;
        for (int i = 0; i < nodeList.size(); i++) {
            chr1Start = _CHROM1_START_ + (nodeList[i]->_val[0]-1) * _RESOLUTION_;
            chr1End =   _CHROM1_START_ + (nodeList[i]->_val[1]-1) * _RESOLUTION_;
            chr2Start = _CHROM2_START_ + (nodeList[i]->_val[0]-1) * _RESOLUTION_;
            chr2End =   _CHROM2_START_ + (nodeList[i]->_val[1]-1) * _RESOLUTION_;
            fprintf(outFile, "%s\t%d\t%d\t%s\t%d\t%d\t%node%d\n",
                _CHROM1_, chr1Start, chr1End, _CHROM2_, chr2Start, chr2End, i);
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


void Writer::writeBoundaries(std::string filePath, std::vector<utils::boundary> &boundaryList)
{
    std::ofstream file;
    file.open(filePath);
    if (file.is_open())
    {
        if (_VERBOSE_)
            printf("start writing boundaries into: %s\n", filePath.c_str());

        for (std::vector<utils::boundary>::iterator it=boundaryList.begin(); it!=boundaryList.end(); it++) {
            for (int i=it->first; i<=it->second; i++)
                file << i << " ";
            file << "\n";
        }
        file.close();

        if (_VERBOSE_)
            printf("finish writing boundaries\n");
        else
            printf("write boundaries into: %s\n", filePath.c_str());
    }
    else
        std::cerr << "cannot open file: " << filePath << "\n";
}


void Writer::dumpMatrix(Eigen::MatrixXd &mat, std::string outPath)
{
    std::ofstream file(outPath);

    if (file.is_open()) {
        printf("start dumping matrix into: %s\n", outPath.c_str());

        file << mat << '\n';
        if (_VERBOSE_)
            printf("finish dumping matrix\n");
        else
            printf("dump matrix into: %s\n", outPath.c_str());

        file.close();
    }
    else
        std::cerr << "cannot open file: " << outPath << "\n";
}


void Writer::dumpCoordinates(i2dMap &map, std::string outPath, std::ofstream *f)
{
    bool append = false;
    if (f) {
        append = true;
    } else {
        std::ofstream file(outPath);
        f = &file;
    }
    if (f->is_open()) {
        if (_VERBOSE_)
            printf("start dumping coordinates into: %s\n", outPath.c_str());

        for (auto it=map.begin(); it!=map.end(); it++) {
            *f << it->first << "\t" << it->second << "\n";
        }
        if (!append)
            f->close();

        if (_VERBOSE_)
            printf("finish dumping coordinates\n");
        else
            printf("dump coordinates into: %s", outPath.c_str());
    }
    else
        std::cerr << "cannot open file: " << outPath << "\n";
}


void Writer::writeListOfCoordinates(str_2_i2dMap &map, std::string outPath)
{
    std::ofstream file(outPath);
    if (file.is_open()) {
        if (_VERBOSE_)
            std::cout << "start dumping list of coordinates into: " << outPath << "\n";
        file << "{";
        for (auto it=map.begin(); it!=map.end(); it++) {
            file << it->first << ":{";
            for (auto it2=it->second.begin(); it2!=it->second.end(); it2++) {
                file << it2->first << ":" << it2->second;
                if (it2->first!=it->second.rbegin()->first)
                    file << ",";
            }
            file << "}";
            if (it->first!=map.rbegin()->first)
                file << ",";
            file << "\n";
        }
        file << "}\n";
        file.close();

        if (_VERBOSE_)
            std::cout << "finish dumping list of coordinates\n";
        else
            std::cout << "dump list of coordinates into: " << outPath << "\n";
    }
    else
        std::cerr << "cannot open file: " << outPath << "\n";
}


void Writer::dumpListOfCoordinates(str_2_i2dMap &map, std::string outPath)
{
    std::ofstream file(outPath);
    if (file.is_open()) {
        if (_VERBOSE_)
            std::cout << "start dumping list of coordinates into: " << outPath << "\n";

        file << "{";
        for (auto it=map.begin(); it!=map.end(); it++) {
            file << "\"" << it->first << "\"" << ":{";
            for (auto it2=it->second.begin(); it2!=it->second.end(); it2++) {
                file << "\"" << it2->first << "\":\"" << it2->second << "\"";
                if (it2->first!=it->second.rbegin()->first)
                    file << ",";
            }
            file << "}";
            if (it->first!=map.rbegin()->first)
                file << ",";
        }
        file << "}";
        file.close();

        if (_VERBOSE_)
            std::cout << "start dumping list of coordinates\n";
        else
            std::cout << "dump list of coordinates into: " << outPath << "\n";
    }
    else
        std::cerr << "cannot dump list of coordinates to " << outPath << "\n";
}
