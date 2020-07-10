//
// Created by wang mengbo on 2019-09-02.
//

#include "inputAndOutput.h"


bool pathExist(const std::string &s)
{
    struct stat buffer;
    return (stat(s.c_str(), &buffer) == 0);
}


void Reader::parseMatrix(Eigen::MatrixXd &contactMat, std::string filePath)
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
        _N_ = 0;
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
                        _N_++;
                    }
                    contactMat.resize(_N_, _N_);
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
    }
    catch (const std::ifstream::failure& e) {
        printf("exception reading file\n");
        exit(1);
    }

}


void Reader::parseMatrix2Table(double **&table, std::string path)
{
    if (path == "") {
        fprintf(stderr, "input must be provided\n");
        exit(1);
    }

    if (!pathExist(path)) {
        fprintf(stderr, "input file not exist\n");
        exit(1);
    }

    std::ifstream file;
    file.exceptions(std::ifstream::badbit);
    try {
        file.open(path);
        if (file.is_open()) {
            if (_VERBOSE_)
                printf("start parsing input: %s\n", path.c_str());
            _N_ = 0;
            std::string line;
            double c;
            getline(file, line);
            std::istringstream iss(line);
            while (iss >> c)
                _N_++;
            iss.clear();
            table = new double *[_N_];

            iss.str(line);
            int i=0, j=0;
            table[i] = new double [_N_]{};
            for (;j<_N_;j++) {
                iss >> c;
                table[i][j] = c;
            }
            iss.clear();

            while (getline(file, line)) {
                iss.str(line);
                table[++i] = new double [_N_]{};
                for (j=0; j<_N_; j++) {
                    iss >> c;
                    table[i][j] = c;
                }
                iss.clear();
            }
        }
    }
    catch (const std::ifstream::failure& e) {
        printf("exception reading file\n");
        exit(1);
    }
}


void Writer::writeBoundaries(std::string filePath, std::vector<boundary> &boundaryList)
{
    std::ofstream file;
    file.open(filePath);
    if (file.is_open())
    {
        if (_VERBOSE_)
            printf("start writing boundaries into: %s\n", filePath.c_str());
        else
            printf("write boundaries into: %s\n", filePath.c_str());

        for (std::vector<boundary>::iterator it=boundaryList.begin(); it!=boundaryList.end(); it++) {
            for (int i=it->first; i<=it->second; i++)
                file << i << " ";
            file << "\n";
        }
        file.close();

        if (_VERBOSE_)
            printf("finish writing boundaries\n");
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
