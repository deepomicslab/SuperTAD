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
#include "compare.h"


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
             "****************************************************************************************\n";
//             "USAGE: " + std::string(argv[0]) + " <input Hi-C matrix> [command] [-option value]\n";
    info += "COMMANDS:\n"

            "\tbinary\tThe first mode requires no user-defined parameters, run the nodes filtering by default\n"
            "\t\t./SuperTAD binary <input Hi-C matrix> [-option values]\n"
            "\t\tOPTIONS:\n"
            "\t\t\t--no-filter: If given, do not filter TADs after TAD detection\n"

            "\tmulti\tThe second mode requires a parameter h to determine the number of layers\n"
            "\t\t./SuperTAD multi <input Hi-C matrix> -h <height> [-option values]\n"
            "\t\tOPTIONS:\n"
            "\t\t\t-h <int>: The height of coding tree, default: 2\n"

            "\tmulti_2d\tThe third mode requires two parameter h1 and h2 to determine the iteractions for dividing and merging\n"
            "\t\t./SuperTAD multi_2d <input Hi-C matrix> [-option values]\n"
            "\t\tOPTIONS:\n"
            "\t\t\t--hd <int>: The height of layers for dividing (go down), default: 2\n"
            "\t\t\t--hu <int>: The height of layers for merging (go up), default: 1\n"
            "\t\t\t--pre <string>: The pre-detected result file\n"

            "\t    SHARED OPTIONS for binary and multi COMMAND:\n"
            "\t\t-K <int>: The number of leaves in the coding tree, default: nan (determined by the algorithm)\n"
            "\t\t--chrom1 <string>: chrom1 label, default: chr1\n"
            "\t\t--chrom2 <string>: chrom2 label, default: the same as chrom1\n"
            "\t\t--chrom1-start <int>: start pos on chrom1, default: 0\n"
            "\t\t--chrom2-start <int>: start pos on chrom2, default: the same as --chrom1-start\n"
            "\t\t-r/--resolution <int>: bin resolution, default: 10000\n"

            "\tfilter\tThe nodes filter for optimal coding tree:\n"
            "\t\t./SuperTAD filter <input Hi-C matrix> -i <original result> \n"
            "\t\tOPTIONS:\n"
            "\t\t\t-i <string>: The list of TAD candidates\n"

            "\tcompare\tThe symmetric metric overlapping ratio to assess the agreement between two results\n"
            "\t\t./SuperTAD compare <result1> <result2>\n"

            "GLOBAL OPTIONS:\n"
            "\t-w <string>: Working directory path, default: the directory where the input file is located\n"
            "\t-v/--verbose: Print verbose\n";

    if (err)
        fprintf(stderr, "%s", info.c_str());
    else
        fprintf(stdout, "%s", info.c_str());

    return 0;
}


int parseArg(int argc, char *argv[], int i)
{
    if (SuperTAD::_BINARY_ || SuperTAD::_MULTI_ || SuperTAD::_FILTER_ || SuperTAD::_MULTI_H_) {
        SuperTAD::_INPUT_ = std::string(*(argv + i));
        printf("input file is %s\n", SuperTAD::_INPUT_.c_str());
    }

    while (i < argc) {
        if (std::string(*(argv + i)) == std::string("--help")) {
            printUsage(argv, 0);
        }

        if (std::string(*(argv + i)) == std::string("-w")) {
            SuperTAD::_WORK_DIR_ = std::string (*(argv + ++i));
            printf("working dir is %s\n", SuperTAD::_WORK_DIR_.c_str());
        }

        if (std::string(*(argv+i))==std::string("--bedpe")) {
            SuperTAD::_BEDPE_ = true;
            printf("output will be written in BEDPE format\n");
        }

        if (std::string(*(argv+i))==std::string("--short")) {
            SuperTAD::_SHORT_ = true;
            printf("output will be written in Juicer short with score format\n");
        }

        if (std::string(*(argv+i))==std::string("--bin-list")) {
            SuperTAD::_BIN_LIST_ = true;
            printf("output will be written as list of bins\n");
        }

        if (std::string(*(argv + i)) == std::string("-v") || std::string(*(argv + i)) == std::string("--verbose")) {
            SuperTAD::_VERBOSE_ = true;
            setbuf(stdout, NULL);
            printf("print verbose\n");
        }

        if (std::string(*(argv + i)) == std::string("-K")) {
            SuperTAD::_K_ = atoi(*(argv + ++i));
            SuperTAD::_DETERMINE_K_ = false;
            printf("set K to %d\n", SuperTAD::_K_);
        }

        if (std::string(*(argv + i)) == std::string("-k")) {
            SuperTAD::_K_ = atoi(*(argv + ++i));
            printf("set max K to %d\n", SuperTAD::_K_);
        }

        if (std::string(*(argv + i)) == std::string("-h")) {
            SuperTAD::_H_ = atoi(*(argv + ++i));
            printf("set H to %d\n", SuperTAD::_H_);
        }

        if (std::string(*(argv + i)) == std::string("--no-filter")) {
            SuperTAD::_FILTERING_ = false;
            printf("disable filtering\n");
        }

        if (std::string(*(argv + i)) == std::string("--hd")) {
            SuperTAD::_HD_ = atoi(*(argv + ++i));
            printf("set the height for going down to %d\n", SuperTAD::_HD_);
        }

        if (std::string(*(argv + i)) == std::string("--hu")) {
            SuperTAD::_HU_ = atoi(*(argv + ++i));
            printf("set the height for going up to %d\n", SuperTAD::_HU_);
        }

        if (std::string(*(argv + i)) == std::string("--pre")) {
            SuperTAD::_PRE_ = std::string(*(argv + ++i));
            printf("the pre-detected result: %s\n", SuperTAD::_PRE_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("-i")) {
            SuperTAD::_RESULT_ = std::string(*(argv + ++i));
            printf("the input result for nodes filtering: %s\n", SuperTAD::_RESULT_.c_str());
        }

        if (std::string(*(argv+i))==std::string("--no-fast")) {
            SuperTAD::_FAST_ = false;
            printf("disable fast mode\n");
        }

        if (std::string(*(argv+i))==std::string("--penalty")) {
            SuperTAD::_PENALTY_ = atoi(*(argv + ++i));
            printf("fast mode penalty is set to %d\n", SuperTAD::_PENALTY_);
        }

        if (std::string(*(argv + i)) == std::string("--chrom1")) {
            SuperTAD::_CHROM1_ = std::string(*(argv + ++i));
            SuperTAD::_CHROM2_ = SuperTAD::_CHROM1_;
            printf("chrom1 lable is %s\n", SuperTAD::_CHROM1_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom2")) {
            SuperTAD::_CHROM2_ = std::string(*(argv + ++i));
            printf("chrom2 lable is %s\n", SuperTAD::_CHROM2_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom1-start")) {
            SuperTAD::_CHROM1_START_ = atol(*(argv + ++i));
            SuperTAD::_CHROM2_START_ = SuperTAD::_CHROM1_START_;
            printf("starting pos on chrom1 is %ld\n", SuperTAD::_CHROM1_START_);
        }

        if (std::string(*(argv + i)) == std::string("--chrom2-start")) {
            SuperTAD::_CHROM2_START_ = atol(*(argv + ++i));
            printf("starting pos on chrom2 is %ld\n", SuperTAD::_CHROM2_START_);
        }

        if (std::string(*(argv + i)) == std::string("-r") || std::string(*(argv + i)) == std::string("--resolution") ) {
            SuperTAD::_RESOLUTION_ = atoi(*(argv + ++i));
            printf("set resolution to %d bp\n", SuperTAD::_RESOLUTION_);
        }

        // debug
        if (std::string(*(argv + i)) == std::string("--no-pre-log")) {
            SuperTAD::_PRE_LOG_ = false;
            printf("do not pre-calculate log-volume table; longer execution time\n");
        }

        i++;
    }

    if (SuperTAD::_WORK_DIR_ == "")
        SuperTAD::_OUTPUT_ = SuperTAD::_INPUT_;
    else {
        int pos = SuperTAD::_INPUT_.rfind("/");
        SuperTAD::_WORK_DIR_ = SuperTAD::_INPUT_.substr(0, pos);
        SuperTAD::_OUTPUT_ = SuperTAD::_WORK_DIR_ + "/" + SuperTAD::_INPUT_.substr(pos + 1);
    }

    if (SuperTAD::_FAST_ && SuperTAD::_BINARY_)
        printf("enable fast mode for binary mode\n");

    if (SuperTAD::_DEBUG_)
        SuperTAD::_VERBOSE_ = true;

    return 0;
}


int parseCommands(int argc, char *argv[])
{
    if (argc < 2) {
        printUsage(argv, 1);
        return 1;
    }
    int i;
    if (std::string(*(argv + 1)) == std::string("compare")) {
        SuperTAD::_COMPARE_ = true;
        SuperTAD::_BINARY_ = false;
        SuperTAD::_MULTI_ = false;
        SuperTAD::_FILTER_ = false;
        printf("calculate overlapping ratio\n");
        SuperTAD::_RESULT_1_ = std::string(*(argv + 2));
        SuperTAD::_RESULT_2_ = std::string(*(argv + 3));
        printf("result 1 is %s\nresult 2 is %s\n", SuperTAD::_RESULT_1_.c_str(), SuperTAD::_RESULT_2_.c_str());
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("binary")) {
        SuperTAD::_BINARY_ = true;
        SuperTAD::_MULTI_ = false;
        SuperTAD::_FILTER_ = false;
        printf("do binary\n");
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("multi")) {
        SuperTAD::_MULTI_ = true;
        SuperTAD::_BINARY_ = false;
        SuperTAD::_FILTER_ = false;
        printf("do multi\n");
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("filter")) {
        SuperTAD::_MULTI_ = false;
        SuperTAD::_BINARY_ = false;
        SuperTAD::_FILTER_ = true;
        printf("do filtering\n");
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("multi_2d")){
        SuperTAD::_MULTI_H_ = true;
        printf("do multi_2d\n");
        i = 2;
    }

    return parseArg(argc, argv, i);
}


int main (int argc, char *argv[])
{
    std::clock_t t;

    if (parseCommands(argc, argv))
        exit(1);

    if (SuperTAD::_VERBOSE_)
        t = std::clock();

    if (SuperTAD::_BINARY_ || SuperTAD::_MULTI_ || SuperTAD::_FILTER_ || SuperTAD::_MULTI_H_) {
        SuperTAD::Data data(SuperTAD::_INPUT_);
        data.init();

        if (SuperTAD::_BINARY_) {
            binary::Detector db(data);
            db.execute();
        }
        else if (SuperTAD::_MULTI_) {
            if (SuperTAD::_H_ == 1) {
                multi::DetectorH1 dm(data);
                dm.execute();
            } else {
                multi::Detector dm(data);
                dm.execute();
            }
        }
        else if (SuperTAD::_FILTER_){
            binary::Detector db(data);
            db.executeFILTER(SuperTAD::_RESULT_);
        }
        else if (SuperTAD::_MULTI_H_){
            multi::detectorH dh(data);
            dh.pipeline(SuperTAD::_PRE_);
        }
    }

    else if (SuperTAD::_COMPARE_) {
        Comparator comparator(SuperTAD::_RESULT_1_, SuperTAD::_RESULT_2_);
        comparator.execute();
    }

    if (SuperTAD::_VERBOSE_)
        printf("task finished in %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);

    return 0;
}
