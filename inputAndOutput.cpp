//
// Created by wang mengbo on 2019-09-02.
//

#include "inputAndOutput.h"


bool pathExist(const std::string &s)
{
    struct stat buffer;
    return (stat(s.c_str(), &buffer) == 0);
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
                printf("start parsing input from %s\n", path.c_str());
            else
                printf("parse input\n");

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

            if (_VERBOSE_)
                printf("finish parsing input\n");

        }
    }
    catch (const std::ifstream::failure& e) {
        printf("exception reading file\n");
        exit(1);
    }
}


void Reader::readBoundariesIntoGraph(std::string path1, std::string path2, std::vector<Boundary> &boundaries1,
                                std::vector<Boundary> &boundaries2, int **&graph)
{
    if (path1=="" || path2=="") {
        fprintf(stderr, "input must be provided\n");
        exit(1);
    }

    if (!pathExist(path1)) {
        fprintf(stderr, "input file %s not exist\n", path1.c_str());
        exit(1);
    }
    if (!pathExist(path2)) {
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


int Reader::parseBoundariesIn8ColsFormat(std::vector<Boundary> &boundaries, std::string path)
{
    std::ifstream file;
    file.exceptions(std::ifstream::badbit);
    try {
        file.open(path);
        if (file.is_open()) {
            if (_VERBOSE_)
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
                    if (c==2)
                        boundary.first = atoll(token.c_str());
                    if (c==7) {
                        boundary.second = atoll(token.c_str())-1;
                        boundary.size = boundary.second - boundary.first+1;
                        break;
                    }
                    c++;
                }
//                printf("s=%d, e=%d, size=%d\n", boundary.first, boundary.second, boundary.size);
                boundaries.push_back(boundary);

                if (determine_chrom == false){
                    iss.str(line);
                    int c = 0;
                    int pos = 0;
                    while (getline(iss, token, '\t')) {
                        if (c==0)
                            _CHROM1_ = token.c_str();
                        else if (c==1)
                            pos = atoi(token.c_str());
                        else if (c==2)
                            _CHROM1_START_ = atoi(token.c_str());
                        else if (c==3) {
                            _RESOLUTION_ = atoi(token.c_str()) - _CHROM1_START_;
                            _CHROM1_START_ = atoi(token.c_str()) - _RESOLUTION_ * pos;
                        }
                        else if (c==4)
                            _CHROM2_ = token.c_str();
                        else if (c==5)
                            pos = atoi(token.c_str());
                        else if (c==7) {
                            _CHROM2_START_ = atoi(token.c_str()) - _RESOLUTION_ * pos;
                        }
                        c++;
                    }
                    determine_chrom = true;
//                    printf("%s, %d, %s, %d, %d\n", _CHROM1_.c_str(), _CHROM1_START_, _CHROM2_.c_str(), _CHROM2_START_, _RESOLUTION_);
                }

                iss.clear();
            }

            iss.str(line);
            c = 0;
            int s, e, res;
            while (getline(iss, token, '\t')) {
                if (c==2)
                    s = atoi(token.c_str());
                if (c==3) {
                    e = atoi(token.c_str());
                    break;
                }
                c++;
            }
            res = e-s;

            if (_VERBOSE_)
                printf("finish parsing input\n");

            return res;
        }
        return 0;
    }
    catch (const std::ifstream::failure& e) {
        printf("exception reading file\n");
        exit(1);
    }
}


void Writer::writeBoundaries(std::string path, std::vector<Boundary> &boundaryList)
{
    std::ofstream file;
    file.open(path);
    if (file.is_open())
    {
        if (_VERBOSE_)
            printf("start writing boundaries into: %s\n", path.c_str());
        else
            printf("write boundaries into: %s\n", path.c_str());

        for (std::vector<Boundary>::iterator it=boundaryList.begin(); it != boundaryList.end(); it++) {
            for (int i=it->first; i<=it->second; i++)
                file << i << " ";
            file << "\n";
        }
        file.close();

        if (_VERBOSE_)
            printf("finish writing boundaries\n");
    }
    else
        std::cerr << "cannot open file: " << path << "\n";
}


void Writer::dumpCoordinates(Int2DoubleMap &map, std::string path, std::ofstream *f)
{
    bool append = false;
    if (f) {
        append = true;
    } else {
        std::ofstream file(path);
        f = &file;
    }
    if (f->is_open()) {
        if (_VERBOSE_)
            printf("start dumping coordinates into: %s\n", path.c_str());

        for (auto it=map.begin(); it!=map.end(); it++) {
            *f << it->first << "\t" << it->second << "\n";
        }
        if (!append)
            f->close();

        if (_VERBOSE_)
            printf("finish dumping coordinates\n");
        else
            printf("dump coordinates into: %s", path.c_str());
    }
    else
        std::cerr << "cannot open file: " << path << "\n";
}


void Writer::writeListOfCoordinates(Str_2_Int2DoubleMap &map, std::string outPath)
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


void Writer::dumpListOfCoordinates(Str_2_Int2DoubleMap &map, std::string outPath)
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
