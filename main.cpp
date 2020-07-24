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
    info =  "****************************************************************************************\n";
    info += "* SuperTAD: [Super]-fast [T]opological [A]ssociating [D]omain package for Hi-C dataset *\n";
    info += "* version: 1.0                                                                        *\n";
    info += "* Bug report to mbwang2016@gmail.com                                                   *\n";
    info += "*                                                                                      *\n";
    info += "* THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS            *\n";
    info += "* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,          *\n";
    info += "* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE          *\n";
    info += "* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER               *\n";
    info += "* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING              *\n";
    info += "* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER                  *\n";
    info += "* DEALINGS IN THE SOFTWARE.                                                            *\n";
    info += "****************************************************************************************\n";
    info += "USAGE: " + std::string(argv[0]) + " <input-path> [-option value]\n";
    info += "OPTIONS:\n";
    info += "\t-w <string>: working directory; If not given, use current input directory)\n";
    info += "\t-b: binary tree version\n";
    info += "\t-m: multiple tree version\n";
    info += "\t-K <int>: number of clusters in candidate coding tree\n";
    info += "\t-k <int>: max number of clusters in candidate coding tree\n";
    info += "\t-h <int>: hierarchy number (default 2)\n";
    info += "\t--no-filter: do not filter TADs\n";
//    info += "\t--no-fast: disable fast mode for binary mode\n";
    info += "\t--chrom1 <string>: chrom1 label\n";
    info += "\t--chrom2 <string>: chrom2 label (if only chrom1 is given, assume chrom1 and 2 are identical)\n";
    info += "\t--chrom1-start <int>: start pos on chrom1\n";
    info += "\t--chrom2-start <int>: start pos on chrom2\n";
    info += "\t-r/--resolution <int>: resolution\n";
    info += "\t-v/--verbose: print verbose\n";

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
            std::cout << "working dir: " << _WORK_DIR_ << std::endl;
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
            std::cout << "print verbose\n";
        }

        if (std::string(*(argv + i)) == std::string("-b")) {
            _BINARY_ = true;
            _MULTI_ = false;
            std::cout << "do binary\n";
        }

        if (std::string(*(argv + i)) == std::string("-m")) {
            _MULTI_ = true;
            _BINARY_ = false;
            std::cout << "do multi\n";
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
            std::cout << "H=" << _H_ << "\n";
        }

        if (std::string(*(argv + i)) == std::string("--no-filter")) {
            _FILTERING_ = false;
            std::cout << "disable filtering\n";
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
            printf("chrom2 lable is %s\n", _CHROM1_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom1-start")) {
            _CHROM1_START_ = atol(*(argv + ++i));
            printf("starting pos on chrom1 is %ld\n", _CHROM1_START_);
        }

        if (std::string(*(argv + i)) == std::string("--chrom2-start")) {
            _CHROM2_START_ = atol(*(argv + ++i));
            printf("starting pos on chrom2 is %ld\n", _CHROM2_START_);
        }

        if (std::string(*(argv + i)) == std::string("-r") || std::string(*(argv + i)) == std::string("--resolution") ) {
            _RESOLUTION_ = atoi(*(argv + ++i));
            printf("set resolution to %dbp\n", _RESOLUTION_);
        }

        // debug
        if (std::string(*(argv + i)) == std::string("--no-pre-log")) {
            _PRE_LOG_ = false;
            printf("test pre-calculate log volume table time\n");
        }

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
    else if (_MULTI_) {
        multi::Detector dm(data);
        dm.execute();
    }
    else if (_MULTI_ && _H_==1) {
        multi::DetectorH1 dm(data);
        dm.execute();
    }

    if (_VERBOSE_)
        printf("task finished in %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);

    return 0;
}
