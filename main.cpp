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


int main (int argc, char *argv[])
{
    std::clock_t t = std::clock ();

    if (argc < 2) {
        std::cerr
            << "***************************************************************************************\n"
            <<"* SuperTAD: [Super]-fast [T]opological [A]ssociating [D]omain package for Hi-C dataset *\n"
            <<"* version: 0.1                                                                         *\n"
            <<"* Bug report to mbwang2016@gmail.com                                                   *\n"
            <<"*                                                                                      *\n"
            <<"* THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS              *\n"
            <<"* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,          *\n"
            <<"* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE          *\n"
            <<"* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER               *\n"
            <<"* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING              *\n"
            <<"* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER                  *\n"
            <<"* DEALINGS IN THE SOFTWARE.                                                            *\n"
            <<"****************************************************************************************\n";
        std::cerr << "USAGE: " << argv[0] << " [-option value]\n";
        std::cerr << "OPTIONS:\n";
//        std::cerr << "\t-f <input path>: Input contact matrix file path\n";
//        std::cerr << "\t-w <working directory path>: Working directory path (default current working directory)\n";
        std::cerr << "\t-b: Binary tree version\n";
        std::cerr << "\t-m: Multiple tree version\n";
        std::cerr << "\t-H: Multiple version (h1 fast mode)\n";
        std::cerr << "\t-k <int>: Number of leaves in candidate coding tree (default NAN)\n";
        std::cerr << "\t-h <int>: Hierarchy number (default 2)\n";
        std::cerr << "\t--filter <true/True/TRUE/false/False/FALSE>: Filter TADs or not (default: true)\n";
        std::cerr << "\t-v/--verbose: Print verbose\n";
        return 0;
    }

    _INPUT_ = std::string(*(argv+1));
    printf("input file path: %s\n", _INPUT_.c_str());

    int i = 2;
    while (i < argc) {
//        if (std::string(*(argv + i)) == std::string("-f")) {
//            _INPUT_ = std::string (*(argv + ++i));
//            std::cout << "input=" << _INPUT_ << std::endl;
//        }

//        if (std::string(*(argv + i)) == std::string("-w")) {
//            _WORK_DIR = std::string (*(argv + ++i));
//            std::cout << "work dir:" << _WORK_DIR << std::endl;
//        }

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
            int h = atoi(*(argv + ++i));
            std::cout << "H=" << h << "\n";
            _H_ = h;
        }

        if (std::string(*(argv + i)) == std::string("--filter")) {
            std::string tmp = std::string (*(argv + ++i));
            if (tmp == "true" or  tmp == "True" or tmp == "TRUE") {
                _FILTERING_ = true;
                std::cout << "enable filtering\n";
            }
            else {
                _FILTERING_ = false;
                std::cout << "disable filtering\n";
            }
        }

        if (std::string(*(argv + i)) == std::string("--tmp-path")) {
            std::string tmp = std::string (*(argv + ++i));
            std::cout << "tmp_path=" << tmp << "\n";
            _TMP_PATH_ = tmp;
        }

        if (std::string(*(argv+i))==std::string("--no-bold")) {
            _BOLD_ = false;
            printf("disable bold mode may cause extra execution time\n");
        }

        i++;
    }

    if (_DEBUG_)
        _VERBOSE_ = true;

    if (_INPUT_ == "") {
        std::cerr << "input must be provided\n";
        exit(1);
    }
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
    t = std::clock() - t;
    std::cout << "running time: " << (float)t/CLOCKS_PER_SEC << "s\n";

    return 0;
}
