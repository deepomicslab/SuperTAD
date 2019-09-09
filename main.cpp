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
    std::cerr << "\t-f: Input Contact Matrix File Path\n";
    std::cerr << "\t-w: Working Directory Path\n";
    std::cerr << "\t-b: Binary Version\n";
    std::cerr << "\t-m: Multiple Version\n";
    std::cerr << "\t-k: \n";
    std::cerr << "\t-h: \n";
    std::cerr << "\t--filter: Do Filtering Or Not (True By Default)";
    return 0;
  }
  
  int i = 0;
  while (i < argc) {
    if (std::string (*(argv + i)) == std::string ("-f")) {
      _INPUT = std::string (*(argv + ++i));
      std::cout << "input: " << _INPUT << std::endl;
    }
  
    if (std::string (*(argv + i)) == std::string ("-w")) {
      _WORK_DIR = std::string (*(argv + ++i));
      std::cout << "work dir:" << _WORK_DIR << std::endl;
    }
  
    if (std::string (*(argv + i)) == std::string ("-b")) {
      _BINARY = true;
      std::cout << "do binary\n";
    }
  
    if (std::string (*(argv + i)) == std::string ("-m")) {
      _MULTI = true;
      std::cout << "do multi\n";
    }
    
    if (std::string (*(argv + i)) == std::string ("-k")) {
      _K = atoi (*(argv + ++i));
      std::cout << "k=" << _K << std::endl;
    }
  
    if (std::string (*(argv + i)) == std::string ("-h")) {
      _H = atoi (*(argv + ++i));
      std::cout << "h=" << _H << std::endl;
    }
    
    if (std::string (*(argv + i)) == std::string ("--filter")) {
      std::string tmp = std::string (*(argv + ++i));
      if (tmp == "true" or  tmp == "True" or tmp == "TRUE") {
        _FILTERING = true;
        std::cout << "do filtering\n";
      }
      else {
        _FILTERING = false;
        std::cout << "no filtering\n";
      }
    }
    
    i++;
  }
  
  if (_INPUT == "") {
    std::cerr << "input must be provided\n";
    exit (1);
  }
  
  Data data (_INPUT);
  data.init ();
  
  if (_BINARY) {
    binary::Detector db (data);
    db.execute ();
  }
  else if (_MULTI) {
    multi::Detector dm (data);
    dm.execute ();
  }
  
  t = std::clock() - t;
  std::cout << "running time: " << (float)t/CLOCKS_PER_SEC << "s\n";
  
  return 0;
}