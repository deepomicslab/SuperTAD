//
// Created by wang mengbo on 2019-09-01.
//

#include <iostream>
#include <ctime>
#include <cstdlib>
#include "params.h"
#include "detectorBinary.h"
#include "detectorBinaryV2.h"
#include "binaryTree.h"
#include "data.h"
#include "detectorMulti.h"
#include "detectorMultiV2.h"
#include "detectorH.h"
#include "compare.h"
#include "detectorDeepBinary.h"
#include "detectorMultiFast.h"

using namespace SuperTAD;


int printUsage(char *argv[], int err)
{
    std::string info;
    info +=  "****************************************************************************************\n"
             "* SuperTAD: [Super]-fast [T]opological [A]ssociating [D]omain package for Hi-C dataset *\n"
             "* version: 2.0                                                                         *\n"
             "* Bug report to ywzhang224@gmail.com                                                   *\n"
             "*                                                                                      *\n"
             "* THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS            *\n"
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
            "\t\t\t--fast: If given, run a more efficient implementation of the second mode with discretization and neighbor searching, SuperTAD-Fast\n"
            "\t\t\t--step <int>: The number of steps for discretization in SuperTAD-Fast, default: the number of bins\n"
            "\t\t\t--window <int> : The size of the searching window in SuperTAD-Fast, default: 5 (bin)\n"

            "\tmulti_2d\tThe third mode requires two parameter h1 and h2 to determine the iteractions for dividing and merging\n"
            "\t\t./SuperTAD multi_2d <input Hi-C matrix> [-option values]\n"
            "\t\tOPTIONS:\n"
            "\t\t\t--hd <int>: The height of layers for dividing (go down), default: 2\n"
            "\t\t\t--hu <int>: The height of layers for merging (go up), default: 1\n"
            "\t\t\t--pre <string>: The pre-detected result file\n"

            "\tdeepbinary\tThe fouth mode requires no user-defined parameters, an updated version of binary\n"
            "\t\t./SuperTAD deepbinary <input Hi-C matrix> [-option values]\n"
            "\t\tOPTIONS:\n"
            "\t\t\t-p/--prune: whether prune binary tree into subtrees; must be set along -k\n"
            "\t\t\t-k <int>: number of subtrees to be pruned into; must be set along --prune\n"

            "\t    SHARED OPTIONS for binary and multi COMMAND:\n"
            "\t\t-K <int>: The number of leaves in the coding tree; default: nan (determined by the algorithm)\n"
            "\t\t--maxbin <int>: The maximum number of bins in a TAD; default: 10000000 (10 Mb)/resolution\n"
            "\t\t--chrom1 <string>: chrom1 label, default: chr1\n"
            "\t\t--chrom2 <string>: chrom2 label, default: the same as chrom1\n"
            "\t\t--chrom1-start <int>: start pos on chrom1, default: 0\n"
            "\t\t--chrom2-start <int>: start pos on chrom2, default: the same as --chrom1-start\n"
            "\t\t-r/--resolution <int>: bin resolution, default: 10000\n"
            "\t\t-s/--sparse: If given, apply the Bayesian corrected version for the sparse input matrix\n"
            "\t\t--bayes: the factor of pseudo count, default: 1 * log(N); must be set along -s\n"

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
    if (_BINARY_ || _MULTI_ || _FILTER_ || _MULTI_H_ || _DEEPBINARY_) {
        _INPUT_ = std::string(*(argv + i));
        printf("Input file is %s\n", _INPUT_.c_str());
    }

    while (i < argc) {
        if (std::string(*(argv + i)) == std::string("--help")) {
            printUsage(argv, 0);
        }

        if (std::string(*(argv + i)) == std::string("-w")) {
            _WORK_DIR_ = std::string (*(argv + ++i));
            printf("Working dir is %s\n", _WORK_DIR_.c_str());
        }

        if (std::string(*(argv+i))==std::string("--bedpe")) {   // three different output formats, default: 8Cols
            _BEDPE_ = true;
            printf("Output will be written in BEDPE format\n");
        }

        if (std::string(*(argv+i))==std::string("--short")) {
            _SHORT_ = true;
            printf("Output will be written in Juicer short with score format\n");
        }

        if (std::string(*(argv+i))==std::string("--bin-list")) {
            _BIN_LIST_ = true;
            printf("Output will be written as list of bins\n");
        }

        if (std::string(*(argv + i)) == std::string("-v") || std::string(*(argv + i)) == std::string("--verbose")) {
            _VERBOSE_ = true;
            setbuf(stdout, NULL);
            printf("Print verbose\n");
        }

        if (std::string(*(argv + i)) == std::string("-s") || std::string(*(argv + i)) == std::string("--sparse")) {
            _SPARSE_ = true;
            printf("Apply the sparse input modification\n");
        }

        if (std::string(*(argv + i)) == std::string("--bayes")) {
            _BAYESFACTOR_ = stod(*(argv + ++i));
            printf("Change the pseudo count to %f * log(N)\n", _BAYESFACTOR_);
        }

        if (std::string(*(argv + i)) == std::string("-K")) {    // for the v1.0, where k is an dimension in dp table
            _K_ = atoi(*(argv + ++i));
            _DETERMINE_K_ = false;
            printf("Set K to %d\n", _K_);
        }

        if (std::string(*(argv + i)) == std::string("-h")) {    // the restriction of tree height
            _H_ = atoi(*(argv + ++i));
            printf("Set H to %d\n", _H_);
        }

        if (std::string(*(argv + i)) == std::string("--no-filter")) {   // for binary mode
            _PRUNE_ = false;
            printf("Disable tree pruning\n");
        }

        if (std::string(*(argv + i)) == std::string("--hd")) {  // for multi_2d mode
            _HD_ = atoi(*(argv + ++i));
            printf("Set the height for going down to %d\n", _HD_);
        }

        if (std::string(*(argv + i)) == std::string("--hu")) {
            _HU_ = atoi(*(argv + ++i));
            printf("Set the height for going up to %d\n", _HU_);
        }

        if (std::string(*(argv + i)) == std::string("--pre")) {
            _PRE_RESULT_ = std::string(*(argv + ++i));
            printf("The pre-detected result: %s\n", _PRE_RESULT_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("-i")) {    // for filter mode
            _PRE_RESULT_ = std::string(*(argv + ++i));
            printf("The input result for nodes filtering: %s\n", _PRE_RESULT_.c_str());
        }

        if (std::string(*(argv+i))==std::string("--fast")) { // for multi mode, run multi-fast
            _FAST_ = true;
            printf("Run SuperTAD-Fast\n");
        }

        if (std::string(*(argv + i)) == std::string("--step")) {    // for SuperTAD-Fast
            _STEP_ = atoi(*(argv + ++i));
            printf("Set step for the SuperTAD-Fast to %d\n", _STEP_);
        }

        if (std::string(*(argv + i)) == std::string("--window")) {  // for SuperTAD-Fast
            _WINDOW_ = atoi(*(argv + ++i));
            _WINDOW_ = int(_WINDOW_ / 2);
            printf("Set window size for the SuperTAD-Fast mode to %d\n", 2*_WINDOW_+1);
        }

        if (std::string(*(argv+i))==std::string("--penalty")) { // to add an extra searching area for the boundary if
            // we restrict the boundary is within the higher structure, under testing ###
            _PENALTY_ = atoi(*(argv + ++i));
            printf("Fast mode penalty is set to %d\n", _PENALTY_);
        }

        if (std::string(*(argv + i)) == std::string("--chrom1")) {  // for real Hi-C data
            _CHROM1_ = std::string(*(argv + ++i));
            _CHROM2_ = _CHROM1_;
            printf("Chrom1 lable is %s\n", _CHROM1_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom2")) {
            _CHROM2_ = std::string(*(argv + ++i));
            printf("Chrom2 lable is %s\n", _CHROM2_.c_str());
        }

        if (std::string(*(argv + i)) == std::string("--chrom1-start")) {
            _CHROM1_START_ = atol(*(argv + ++i));
            _CHROM2_START_ = _CHROM1_START_;
            printf("Starting pos on chrom1 is %ld\n", _CHROM1_START_);
        }

        if (std::string(*(argv + i)) == std::string("--chrom2-start")) {
            _CHROM2_START_ = atol(*(argv + ++i));
            printf("Starting pos on chrom2 is %ld\n", _CHROM2_START_);
        }

        if (std::string(*(argv + i)) == std::string("-r") || std::string(*(argv + i)) == std::string("--resolution") ) {
            _RESOLUTION_ = atoi(*(argv + ++i));
            _MAX_SIZE_ = int(10000000 / _RESOLUTION_);   // the max domain size is set to 10 Mb
            printf("Set resolution to %d bp, max size of domain to %d bins (10 Mb)\n", _RESOLUTION_, _MAX_SIZE_);
        }

        if (std::string(*(argv + i)) == std::string("--maxbin")){
            _MAX_SIZE_ = atoi(*(argv + ++i));
            printf("Set the maximum number of bins within a TAD to %d bins\n", _MAX_SIZE_);
        }

        if (std::string(*(argv + i)) == std::string("-p") || std::string(*(argv + i)) == std::string("--prune") ) { // for deepbinary mode
            _PRUNE_ = true;
            printf("Deep binary tree will be pruned for best fit\n");
        }

        if (std::string(*(argv + i)) == std::string("--turbo-prune")) { // under testing ###
            _TURBO_PRUNE_ = true;
            printf("Enable turbo prune\n");
        }

        if (std::string(*(argv + i)) == std::string("--prune_method") ) {   // for deepbinary mode
            _PRUNE_METHOD_ = atoi(*(argv+ ++i));
            switch (_PRUNE_METHOD_) {
                case binary::PruneMethod1 :
                    printf("prune method 1\n");
                    break;
                case binary::PruneMethod2 :
                    printf("prune method 2\n");
                    break;
                default:
                    printf("undefined prune method; use default");
            }
        }

        if (std::string(*(argv + i)) == std::string("--v1") ) { // for running version1.0 where K is kept
            _V1_ = true;
            printf("Run version1.0 where the number of leaves K is kept\n");
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

    if (_DEBUG_)
        _VERBOSE_ = true;

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
        _COMPARE_ = true;
        printf("Calculate overlapping ratio\n");
        _RESULT_1_ = std::string(*(argv + 2));
        _RESULT_2_ = std::string(*(argv + 3));
        printf("result 1 is %s\nresult 2 is %s\n", _RESULT_1_.c_str(), _RESULT_2_.c_str());
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("binary")) {
        _BINARY_ = true;
        printf("Do binary\n");
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("multi")) {
        _MULTI_ = true;
        printf("Do multi\n");
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("filter")) {
        _FILTER_ = true;
        printf("Do filtering\n");
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("multi_2d")){
        _MULTI_H_ = true;
        printf("Do multi_2d\n");
        i = 2;
    }

    else if (std::string(*(argv + 1)) == std::string("deepbinary")){
        _DEEPBINARY_ = true;
        printf("Do deepbinary\n");
        i = 2;
    }

    return parseArg(argc, argv, i);
}


int main (int argc, char *argv[])
{
    std::clock_t t;

    if (parseCommands(argc, argv))
        exit(1);

    if (_VERBOSE_)
        t = std::clock();
    Data data(_INPUT_);
    data.init();

    if (_BINARY_ || _MULTI_ || _FILTER_ || _MULTI_H_ || _DEEPBINARY_){
        Data data(_INPUT_);
        data.init();

        if (_BINARY_) {
            if (true)   // default: run V1, as V2 derive the deepbinary tree
            {
                binary::Detector db(data);
                db.execute();
            } else
            {
                binary::BinaryV2 db(data);
                db.execute();
            }
        }else if (_MULTI_) {
            if (_H_ == 1){
                multi::PartitionV2 dm(data);
                dm.execute(-1);
            } else if (_FAST_){
                multi::Discretization dm(data);
                dm.execute();

                multi::NeighborSearch dn(data, &dm._multiTree, dm._tableF);
                dn.execute();


            } else if (_V1_ or not _DETERMINE_K_){
                multi::Detector dm(data);
                dm.execute();
            } else{
                multi::MultiV2 dm(data);
                dm.execute();
            }
        }else if (_FILTER_){
            binary::Detector db(data);
            db.executeFilter(_PRE_RESULT_);
        }else if (_MULTI_H_){
            multi::detectorH dh(data);
            dh.pipeline(_PRE_RESULT_);
        }else if (_DEEPBINARY_){
            deepBinary::Detector ddb(data);
            ddb.execute();
        }
    }else if (_COMPARE_)
    {
        Comparator comparator(_RESULT_1_, _RESULT_2_);
        comparator.execute();
    }

    if (_VERBOSE_)
        printf("Task finished in %fs\n", (float)(std::clock() - t)/CLOCKS_PER_SEC);

    return 0;
}
