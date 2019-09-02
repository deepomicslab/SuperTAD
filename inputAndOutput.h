//
// Created by wang mengbo on 2019-09-02.
//

#ifndef PROGRAM_INPUTANDOUTPUT_H
#define PROGRAM_INPUTANDOUTPUT_H

#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "binaryTree.h"


std::string concatePath (std::string path1, std::string path2);


class Reader {

};


class Writer {
private:
  std::ofstream _outfile;

public:
  Writer ();

  ~Writer ();
  
  void writeTree (std::string workDir, std::string fileName, std::vector<TreeNode *> &nodeList);
  
};

#endif //PROGRAM_INPUTANDOUTPUT_H
