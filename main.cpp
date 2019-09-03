#include <iostream>
#include "params.h"
#include "detectorBinary.h"
#include "binaryTree.h"
#include "data.h"
#include "utils.h"


int main (int argc, char *argv[])
{
  int i = 0;
  while (i < argc) {
    if (std::string (*(argv + i)) == std::string ("-k")) {
      _K = atoi (*(argv + ++i));
      std::cout << "k=" << _K << std::endl;
    }
    
    if (std::string (*(argv + i)) == std::string ("-f")) {
      _INPUT = std::string (*(argv + ++i));
      std::cout << "input: " << _INPUT << std::endl;
    }
    
    if (std::string (*(argv + i)) == std::string ("-w")) {
      _WORK_DIR = std::string (*(argv + ++i));
      std::cout << "work dir:" << _WORK_DIR << std::endl;
    }
  
    if (std::string (*(argv + i)) == std::string ("--filter")) {
      std::cout << "do filtering\n";
    }
    
    i++;
  }
  
  if (_INPUT == "")
    std::cerr << "input must be provided\n";
  
  Data data (_INPUT);
  data.init ();
  
  DetectorBinary db (data);
  db.execute ();
  
  return 0;
}