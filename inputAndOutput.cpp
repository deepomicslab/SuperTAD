//
// Created by wang mengbo on 2019-09-02.
//

#include "inputAndOutput.h"


void SuperTAD::Reader::parseInput(double **&table, std::string path)
{
    if (path == "") {
        fprintf(stderr, "input must be provided\n");
        exit(1);
    }

    if (!SuperTAD::pathExist(path)) {
        fprintf(stderr, "input file not exist\n");
        exit(1);
    }

    std::ifstream file;
    file.exceptions(std::ifstream::badbit);
    try {
        file.open(path);
        _N_ = 0;
        int m = 0;
        std::string line;
        std::string ct;
        getline(file, line);
        m++;
        std::istringstream iss(line);
        while (iss >> ct) {
            _N_++;
        }
        iss.clear();
        while (getline(file, line)) {
            if (line.size()>0) {
                m++;
            }
        }
        file.close();

        if (_N_==m) {
            file.open(path);
            if (file.is_open()) {
                if (SuperTAD::_VERBOSE_)
                    printf("start parsing input from %s\n", path.c_str());

                SuperTAD::_N_ = 0;
                std::string line;
                std::string ct;
                getline(file, line);
                std::istringstream iss(line);
                while (iss >> ct) {
                    SuperTAD::_N_++;
                }
                iss.clear();
                table = new double *[SuperTAD::_N_];

                iss.str(line);
                int i = 0, j = 0;
                table[i] = new double[SuperTAD::_N_]{};
                for (; j < SuperTAD::_N_; j++) {
                    iss >> ct;
                    if (std::isnormal(std::stod(ct)))
                        table[i][j] = std::stod(ct);
                }
                iss.clear();

                while (getline(file, line)) {
                    iss.str(line);
                    table[++i] = new double[SuperTAD::_N_]{};
                    for (j = 0; j < SuperTAD::_N_; j++) {
                        iss >> ct;
                        if (std::isnormal(std::stod(ct)))
                            table[i][j] = std::stod(ct);;
                    }
                    iss.clear();
                }
            }
        } else {
            std::string line;
            int id;
            m = 0;
            std::set<int> labels;
            file.open(path);
            if (file.is_open()) {
                while (std::getline(file, line)) {
                    if (line[0] == '#') {
                        continue;
                    }
                    iss.str(line);
                    for (int i = 0; i < 2; i++) {
                        iss >> id;
                        labels.emplace(id);
                    }
                    m++;
                    iss.clear();
                }
                file.close();
                _N_ = labels.size();
                printf("found %d records, %d samples\n", m, _N_);
            } else {
                fprintf(stderr, "cannot open file %s\n", path.c_str());
                exit(1);
            }

            int minus = 0;
            if (*labels.begin() == 1) {
                minus = 1;
            }

            table = new double *[_N_]{};
            for (int i = 0; i < _N_; i++) {
                table[i] = new double [_N_]{};
            }

            file.open(path);
            int i, j;
            double v;
            if (file.is_open()) {
                while (std::getline(file, line)) {
                    iss.str(line);
                    iss >> i >> j >> v;
                    table[i - minus][j - minus] = table[j - minus][i - minus] = v;
                    iss.clear();
                }
                file.close();
            } else {
                fprintf(stderr, "cannot open file %s\n", path.c_str());
                exit(1);
            }
        }
        if (SuperTAD::_VERBOSE_)
            printf("finish parsing input\n");
    }
    catch (const std::ifstream::failure& e) {
        printf("exception reading file\n");
        exit(1);
    }
}


void SuperTAD::Reader::readBoundariesIntoGraph(std::string path1, std::string path2, std::vector<Boundary> &boundaries1,
                                               std::vector<Boundary> &boundaries2, int **&graph)
{
    if (path1=="" || path2=="") {
        fprintf(stderr, "input must be provided\n");
        exit(1);
    }

    if (!SuperTAD::pathExist(path1)) {
        fprintf(stderr, "input file %s not exist\n", path1.c_str());
        exit(1);
    }
    if (!SuperTAD::pathExist(path2)) {
        fprintf(stderr, "input file %s not exist\n", path2.c_str());
        exit(1);
    }

    Reader::parseBoundariesIn8ColsFormat(boundaries1, path1);
    Reader::parseBoundariesIn8ColsFormat(boundaries2, path2);
    int n1 = boundaries1.size();
    int n2 = boundaries2.size();
    int n = n1 + n2 + 2;
    graph = new int *[n];

    for (int i=0; i<n; i++)
        graph[i] = new int [n]{0};

    for (int i=1; i<n1+1; i++) {
        graph[0][i] = std::numeric_limits<int>::max();

    }

    for (int i=n1+1; i<n-1; i++) {
        graph[i][n-1] = std::numeric_limits<int>::max();
    }

    for (int i=0; i<boundaries1.size(); i++) {
        for (int j=0; j<boundaries2.size(); j++) {
            int d = utils::boundariesIntersection(boundaries1[i], boundaries2[j]);
            graph[i+1][n1+1+j] = d;
//            if (d > 0) {
//                printf("graph[%d][%d]=%d\n", i+1, n1+1+j, d); fflush(stdout);
//            }
        }
    }

//    utils::print2Darray(graph, n, n);
}


int SuperTAD::Reader::parseBoundariesIn8ColsFormat(std::vector<Boundary> &boundaries, std::string path)
{
    std::ifstream file;
    file.exceptions(std::ifstream::badbit);
    try {
        file.open(path);
        if (file.is_open()) {
            if (SuperTAD::_VERBOSE_)
                printf("start parsing input from %s\n", path.c_str());
//            else
//                printf("parse input\n");
            std::string line, token;
            std::istringstream iss;
            Boundary boundary;
            bool determine_chrom = false;
            long long int c;
            while (getline(file, line)) {
                iss.str(line);
                c = 0;
                while (getline(iss, token, '\t')) {
                    if (SuperTAD::_COMPARE_) {
                        if (c == 2)
                            boundary.first = atoll(token.c_str());
                        if (c == 7) {
                            boundary.second = atoll(token.c_str()) - 1;
                            boundary.size = boundary.second - boundary.first + 1;
                            break;
                        }
                        c++;
                    } else if (SuperTAD::_FILTER_ || SuperTAD::_MULTI_H_) {
                        if (c == 1)
                            boundary.first = atoll(token.c_str());
                        if (c == 5) {
                            boundary.second = atoll(token.c_str());
                            boundary.size = boundary.second - boundary.first + 1;
                            break;
                        }
                        c++;
                    }
                }
//                printf("s=%d, e=%d, size=%d\n", boundary.first, boundary.second, boundary.size);
                boundaries.push_back(boundary);

                if (determine_chrom == false) {
                    iss.str(line);
                    int c = 0;
                    int pos = 0;
                    while (getline(iss, token, '\t')) {
                        if (c == 0)
                            SuperTAD::_CHROM1_ = token.c_str();
                        else if (c == 1)
                            pos = atoi(token.c_str());
                        else if (c == 2)
                            SuperTAD::_CHROM1_START_ = atoi(token.c_str());
                        else if (c == 3) {
                            SuperTAD::_RESOLUTION_ = atoi(token.c_str()) - SuperTAD::_CHROM1_START_;
                            SuperTAD::_CHROM1_START_ = atoi(token.c_str()) - SuperTAD::_RESOLUTION_ * pos;
                        } else if (c == 4)
                            SuperTAD::_CHROM2_ = token.c_str();
                        else if (c == 5)
                            pos = atoi(token.c_str());
                        else if (c == 7) {
                            SuperTAD::_CHROM2_START_ = atoi(token.c_str()) - SuperTAD::_RESOLUTION_ * pos;
                        }
                        c++;
                    }
                    determine_chrom = true;
//                    printf("%s, %d, %s, %d, %d\n", _CHROM1_.c_str(), _CHROM1_START_, _CHROM2_.c_str(), _CHROM2_START_,_RESOLUTION_);
                }

                iss.clear();
            }
            if (SuperTAD::_VERBOSE_)
                printf("finish parsing input\n");
            return SuperTAD::_RESOLUTION_;
            }
        else throw "exception reading file";
    }
//    catch (const std::ifstream::failure& e) {
    catch (...) {
        printf("exception reading file\n");
        exit(1);
    }
}


void SuperTAD::Writer::writeBoundaries(std::string path, std::vector<Boundary> &boundaryList)
{
    std::ofstream file;
    file.open(path);
    if (file.is_open())
    {
        if (SuperTAD::_VERBOSE_)
            printf("start writing boundaries into: %s\n", path.c_str());
        else
            printf("write boundaries into: %s\n", path.c_str());

        for (std::vector<Boundary>::iterator it=boundaryList.begin(); it != boundaryList.end(); it++) {
            for (int i=it->first; i<=it->second; i++)
                file << i << " ";
            file << "\n";
        }
        file.close();

        if (SuperTAD::_VERBOSE_)
            printf("finish writing boundaries\n");
    }
    else
        std::cerr << "cannot open file: " << path << "\n";
}


void SuperTAD::Writer::writeBoundIn8Cols(std::string path, std::vector<Boundary> &boundaryList) {
    if (_SPARSE_)
        path += "_sparse.tsv";
    else
        path += ".tsv";
    FILE *outFile = NULL;
    outFile = std::fopen(path.c_str(), "w");
    if (outFile)
    {
        if (SuperTAD::_VERBOSE_)
            printf("start writing boundaries into: %s\n", path.c_str());
        else
            printf("write boundaries into: %s\n", path.c_str());
        int bin1Idx, bin1Start, bin1End, bin2Idx, bin2Start, bin2End;
        for (std::vector<Boundary>::iterator it=boundaryList.begin(); it != boundaryList.end(); it++) {
            bin1Idx = it->first;
            bin1Start = SuperTAD::_CHROM1_START_ + (bin1Idx - 1) * SuperTAD::_RESOLUTION_;
            bin1End = SuperTAD::_CHROM1_START_ + bin1Idx * SuperTAD::_RESOLUTION_;
            bin2Idx = it->second;
            bin2Start = SuperTAD::_CHROM1_START_ + (bin2Idx - 1) * SuperTAD::_RESOLUTION_;
            bin2End = SuperTAD::_CHROM1_START_ + bin2Idx * SuperTAD::_RESOLUTION_;
            fprintf(outFile, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n",
                    SuperTAD::_CHROM1_.c_str(), bin1Idx, bin1Start, bin1End, SuperTAD::_CHROM2_.c_str(), bin2Idx, bin2Start, bin2End);
        }
        fclose(outFile);

        if (SuperTAD::_VERBOSE_)
            printf("finish writing boundaries\n");
    }
    else
        std::cerr << "cannot open file: " << path << "\n";
}


void SuperTAD::Writer::dumpCoordinates(Int2DoubleMap &map, std::string path, std::ofstream *f)
{
    bool append = false;
    if (f) {
        append = true;
    } else {
        std::ofstream file(path);
        f = &file;
    }
    if (f->is_open()) {
        if (SuperTAD::_VERBOSE_)
            printf("start dumping coordinates into: %s\n", path.c_str());

        for (auto it=map.begin(); it!=map.end(); it++) {
            *f << it->first << "\t" << it->second << "\n";
        }
        if (!append)
            f->close();

        if (SuperTAD::_VERBOSE_)
            printf("finish dumping coordinates\n");
        else
            printf("dump coordinates into: %s", path.c_str());
    }
    else
        std::cerr << "cannot open file: " << path << "\n";
}


void SuperTAD::Writer::writeListOfCoordinates(Str_2_Int2DoubleMap &map, std::string outPath)
{
    std::ofstream file(outPath);
    if (file.is_open()) {
        if (SuperTAD::_VERBOSE_)
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

        if (SuperTAD::_VERBOSE_)
            std::cout << "finish dumping list of coordinates\n";
        else
            std::cout << "dump list of coordinates into: " << outPath << "\n";
    }
    else
        std::cerr << "cannot open file: " << outPath << "\n";
}


void SuperTAD::Writer::dumpListOfCoordinates(Str_2_Int2DoubleMap &map, std::string outPath)
{
    std::ofstream file(outPath);
    if (file.is_open()) {
        if (SuperTAD::_VERBOSE_)
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

        if (SuperTAD::_VERBOSE_)
            std::cout << "start dumping list of coordinates\n";
        else
            std::cout << "dump list of coordinates into: " << outPath << "\n";
    }
    else
        std::cerr << "cannot dump list of coordinates to " << outPath << "\n";
}
