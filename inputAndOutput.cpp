//
// Created by wang mengbo on 2019-09-02.
//

#include "inputAndOutput.h"


std::string concatePath (std::string path1, std::string path2)
{
//  if (path1[path1.length ()-11] == '/')
//    path1 = path1[]
  return "";
}


void Writer::writeTree (std::string workDir, std::string fileName, std::vector<TreeNode *> &nodeList)
{
  std::string path = workDir + fileName;
  _outfile.open (path);
  for (int i = 0; i < nodeList.size (); i++) {
    for (int j = nodeList[i]->_val[0]; j <= nodeList[i]->_val[1]; j++)
      _outfile << std::to_string (j + 1) << " ";
    _outfile << "\n";
  }
  _outfile.close ();
}

