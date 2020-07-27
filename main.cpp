//
// Created by wang mengbo on 2019-09-01.
//

#include <iostream>
#include "params.h"
#include "detectorBinary.h"
#include "binaryTree.h"
#include "data.h"
#include <ctime>
#include "detectorMulti.h"
#include <cstdlib>
#include "detectorH.h"


int printUsage(char *argv[], int err)
{

    std::string info;
    info +=  "****************************************************************************************\n"
             "* SuperTAD: [Super]-fast [T]opological [A]ssociating [D]omain package for Hi-C dataset *\n"
             "* version: 1.0                                                                         *\n"
             "* Bug report to mbwang2016@gmail.com                                                   *\n"
             "*                                                                                      *\n"
             "* THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS              *\n"
             "* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,          *\n"
             "* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE          *\n"
             "* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER               *\n"
             "* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING              *\n"
             "* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER                  *\n"
             "* DEALINGS IN THE SOFTWARE.                                                            *\n"
             "****************************************************************************************\n"
             "USAGE: " + std::string(argv[0]) + " <input Hi-C matrix> [command] [-option value]\n";
    info += "COMMANDS:\n"
            "\tbinary \tThe first mode requires no user-defined parameters, run the nodes filtering by default\n"
            "\t         --no-filter \t if given, do not filter TADs\n"
            "\tfilter \tThe nodes filter for optimal coding tree:\n"
            "\t         ./SuperTAD <input Hi-C matrix> filter -i <original result> [-option values]\n"
            "\t         -i <string> \t The list of TAD candidates\n"
            "\tmulti \tThe second mode requires a parameter h to determine the number of layers\n"
            "\t         -h <int> \t The height of coding tree, default: 2\n"
            "GLOBAL COMMAND OPTIONS:\n"
            "\t-w <string>: working directory, default: the directory where the input file is located\n"
            "\t-K <int>: The number of leaves in the coding tree, default: nan (determined by the algorithm)\n"
            "\t--chrom1 <string>: chrom1 label, default: chr1\n"
            "\t--chrom2 <string>: chrom2 label, default: the same as chrom1\n"
            "\t--chrom1-start <int>: start pos on chrom1, default: 0\n"
            "\t--chrom2-start <int>: start pos on chrom2, default: the same as --chrom1-start\n"
            "\t-r/--resolution <int>: bin resolution, default: 10000\n"
            "\t-v/--verbose: print verbose\n";

    if (err)
        fprintf(stderr, "%s", info.c_str());
    else
        fprintf(stdout, "%s", info.c_str());

    return 0;
}


int parseArg(int argc, char *argv[])
{
    if (argc < 2) {
        printUsage(argv, 1);
        return 1;
    }

    _INPUT_ = std::string(*(argv+1));
    printf("input file is %s\n", _INPUT_.c_str());

    int i = 2;
    while (i < argc) {
        if (std::string(*(argv + i)) == std::string("--help")) {
            printUsage(argv, 0);
        }

        if (std::string(*(argv + i)) == std::string("-w")) {
            _WORK_DIR_ = std::string (*(argv + ++i));
            printf("working dir is %s\n", _WORK_DIR_.c_str());
        }

        if (std::string(*(argv+i))==std::string("--bedpe")) {
            _BEDPE_ = true;
            printf("output will be written in BEDPE format\n");
        }

        if (std::string(*(argv+i))==std::string("--short")) {
            _SHORT_ = true;
            printf("output will be written in Juicer short with score format\n");
        }

        if (std::string(*(argv+i))==std::string("--bin-list")) {
            _BIN_LIST_ = true;
            printf("output will be written as list of bins\n");
        }

        if (std::string(*(argv + i)) == std::string("-v") || std::string(*(argv + i)) == std::string("--verbose")) {
            _VERBOSE_ = true;
            setbuf(stdout, NULL);
            printf("print verbose\n");
        }

        if (std::string(*(argv + i)) == std::string("binary")) {
            _BINARY_ = true;
            _MULTI_ = false;
            printf("do binary\n");
        }

        if (std::string(*(argv + i)) == std::string("multi")) {
            _MULTI_ = true;
            _BINARY_ = false;
            printf("do multi\n");
        }

        if (std::string(*(argv + i)) == std::string("filter")) {
            _MULTI_ = false;
            _BINARY_ = false;
            _FILTER_ = true;
            printf("do filtering\n");
        }

        if (std::string(*(argv + i)) == std::string("-K")) {
            _K_ = atoi(*(argv + ++i));
            _DETERMINE_K_ = false;
            printf("set K to %d\n", _K_);
        }

        if (std::string(*(argv + i)) == std::string("-k")) {
            _K_ = atoi(*(argv + ++i));
            printf("set max K to %d\n", _K_);
        }

        if (std::string(*(argv + i)) == std::string("-h")) {
            _H_ = atoi(*(argv + ++i));
            printf("set H to %d\n", _H_);
        }

        if (std::string(*(argv + i)) == std::string("--no-filter")) {
            _FILTERING_ = false;
            printf("disable filtering\n");
        }

        if (std::string(*(argv+i))==std::string("--no-fast")) {
            _FAST_ = false;
            printf("disable fast mode\n");
        }

        if (std::string(*(argv+i))==std::string("--penalty")) {
            _PENALTY_ = atoi(*(argv + ++i));
            printf("fast mode penalty is set to %d\n", _PENALTY_);
        }

        if (std::string(*(argv + i)) == std::string("--chrom1")) {
            _CHROM1_ = std::string(*(argv + ++i));
            _CHROM2_ = _CHROM1_;
            printf("chrom1 lable is %s\n", _CHROM1_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom2")) {
            _CHROM2_ = std::string(*(argv + ++i));
            printf("chrom2 lable is %s\n", _CHROM2_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom1-start")) {
            _CHROM1_START_ = atol(*(argv + ++i));
            _CHROM2_START_ = _CHROM1_START_;
            printf("starting pos on chrom1 is %ld\n", _CHROM1_START_);
        }

        if (std::string(*(argv + i)) == std::string("--chrom2-start")) {
            _CHROM2_START_ = atol(*(argv + ++i));
            printf("starting pos on chrom2 is %ld\n", _CHROM2_START_);
        }

        if (std::string(*(argv + i)) == std::string("-r") || std::string(*(argv + i)) == std::string("--resolution") ) {
            _RESOLUTION_ = atoi(*(argv + ++i));
            printf("set resolution to %d bp\n", _RESOLUTION_);
        }

        // debug
        if (std::string(*(argv + i)) == std::string("--no-pre-log")) {
            _PRE_LOG_ = false;
            printf("do not pre-calculate log-volume table; longer execution time\n");
        }

        i++;
    }

    if (_WORK_DIR_ == "")
        _OUTPUT_ = _INPUT_;
    else {
        int pos = _INPUT_.rfind("/");
        _WORK_DIR_ = _INPUT_.substr(0, pos);
        _OUTPUT_ = _WORK_DIR_ + "/" + _INPUT_.substr(pos + 1);
//        printf("_OUTPUT_=%s\n", _OUTPUT_.c_str());
    }

    if (_FAST_ && _BINARY_)
        printf("enable fast mode for binary mode\n");

    if (_DEBUG_)
        _VERBOSE_ = true;

    return 0;
}


int main (int argc, char *argv[])
{
    std::clock_t t;

    if (parseArg(argc, argv))
        exit(1);

    if (_VERBOSE_)
        t = std::clock();

    Data data(_INPUT_);
    data.init();

    if (_BINARY_) {
        binary::Detector db(data);
        db.execute();
    }
    else if (_MULTI_ && _H_ == 1) {
        multi::DetectorH1 dm(data);
        dm.execute();
    }
    else if (_MULTI_ ) {
        multi::Detector dm(data);
        dm.execute();
    }
    else if (_FILTER_) {
        //????????
    }

    if (_VERBOSE_)
        printf("task finished in %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);

    return 0;
}
