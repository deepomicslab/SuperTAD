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
    info += "* version: 1.1                                                                         *\n";
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
    info += "\t-H: multiple version (h1 fast mode)\n";
    info += "\t-k <int>: number of leaves in candidate coding tree (default NAN)\n";
    info += "\t-h <int>: hierarchy number (default 2)\n";
    info += "\t--no-filter: do not filter TADs\n";
    info += "\t--no-bold: disable bold mode\n";
    info += "\t--no-fast: disable fast mode for binary mode\n";
    info += "\t--bedpe: write output in BEDPE format\n";
    info += "\t--chrom1 <string>: chrom1 label\n";
    info += "\t--chrom2 <string>: chrom2 label (if only chrom1 is given, assume chrom1 and 2 are identical)\n";
    info += "\t--chrom1-start <int>: start pos on chrom1\n";
    info += "\t--chrom2-start <int>: start pos on chrom2\n";
    info += "\t-r/--resolution <int>: resolution\n";
    info += "\t-v/--verbose: print verbose\n";

    if (err)
        fprintf(stderr, info.c_str());
    else
        fprintf(stdout, info.c_str());

    return 0;
}


int parseArg(int argc, char *argv[])
{
    if (argc < 2) {
        printUsage(argv, 1);
        return 1;
    }

    _INPUT_ = std::string(*(argv+1));
    printf("input path: %s\n", _INPUT_.c_str());

    int i = 2;
    while (i < argc) {
        if (std::string(*(argv + i)) == std::string("--help")) {
            printUsage(argv, 0);
        }

        if (std::string(*(argv + i)) == std::string("-w")) {
            _WORK_DIR_ = std::string (*(argv + ++i));
            std::cout << "working dir: " << _WORK_DIR_ << std::endl;
        }

        if (std::string(*(argv + i)) == std::string("-v") || std::string(*(argv + i)) == std::string("--verbose")) {
            _VERBOSE_ = true;
            setbuf(stdout, NULL);
            std::cout << "print verbose\n";
        }

        if (std::string(*(argv + i)) == std::string("-b")) {
            _BINARY_ = true;
            std::cout << "do binary\n";
        }

        if (std::string(*(argv + i)) == std::string("-m")) {
            _MULTI_ = true;
            std::cout << "do multi\n";
        }

        if (std::string(*(argv+i)) == std::string("-H")) {
            _MULTI_H_ = true;
            std::cout << "do new_multi\n";
        }

        if (std::string(*(argv + i)) == std::string("-k")) {
            _K_ = atoi(*(argv + ++i));
            _DETERMINE_K_ = false;
            std::cout << "K=" << _K_ << "\n";
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

        if (std::string(*(argv+i))==std::string("--bedpe")) {
            _BEDPE_ = true;
            printf("output will be written in BEDPE format\n");
        }

        if (std::string(*(argv + i)) == std::string("--chrom1")) {
            _CHROM1_ = std::string(*(argv + ++i));
            _CHROM2_ = _CHROM1_;
            printf("chrom1 lable is given as: %s", _CHROM1_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom2")) {
            _CHROM2_ = std::string(*(argv + ++i));
            printf("chrom1 lable is given as: %s", _CHROM1_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom1-start")) {
            _CHROM1_START_ = atoi(*(argv + ++i));
            printf("starting pos on chrom1: %d", _CHROM1_START_);
        }

        if (std::string(*(argv + i)) == std::string("--chrom2-start")) {
            _CHROM2_START_ = atoi(*(argv + ++i));
            printf("starting pos on chrom2: %d", _CHROM2_START_);
        }

        // debug
        if (std::string(*(argv + i)) == std::string("--test-log-time")) {
            _TEST_LOG2_TIME_ = true;
            printf("test log2 execution time\n");
        }

        if (std::string(*(argv + i)) == std::string("--test-pre-log")) {
            _TEST_LOG_VOL_TABLE_ = true;
            printf("test pre-calculate log volume table time\n");
        }

        if (std::string(*(argv + i)) == std::string("--tmp-path")) {
            std::string tmp = std::string(*(argv + ++i));
            std::cout << "tmp_path=" << tmp << "\n";
            _TMP_PATH_ = tmp;
        }

        i++;
    }

    if (_WORK_DIR_ == "")
        _OUTPUT_ = _INPUT_;
    else {
        int pos = _INPUT_.rfind("/");
        _WORK_DIR_ = _INPUT_.substr(0, pos);
        _OUTPUT_ = _WORK_DIR_ + "/" + _INPUT_.substr(pos + 1);
    }

    if (_FAST_)
        printf("enable fast mode\n");

    return 0;
}


int main (int argc, char *argv[])
{
    std::clock_t t = std::clock();

    if (parseArg(argc, argv))
        exit(1);

    if (_DEBUG_)
        _VERBOSE_ = true;

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
    else if (_MULTI_H_) {
        multi::DetectorH1 dH(data);
        dH.execute();
    }

    printf("running time: %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);

    return 0;
}
