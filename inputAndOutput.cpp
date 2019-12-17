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


bool isPathExist (const std::string &s)
{
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}


Reader::Reader (std::string fileName)
{
  _fileName = fileName;
}


int Reader::parse (Eigen::MatrixXd &contactMat, std::string fileName)
{
  std::string line;
  double c;
  
  if (fileName == "")
    fileName = _fileName;
  
  if (!isPathExist (fileName)) {
    std::cerr << "file '" << fileName << "'not exist\n";
    exit (1);
  }
    
  _infile.open (fileName);
  
  bool init = false;
  int count = 0;
  int pos1 = 0;
  while (getline (_infile, line)) {
    std::istringstream iss (line);
    if (!init) {
      std::string lineBack = line;
      std::istringstream issBack (lineBack);
      while (issBack >> c) {
        count++;
      }
      contactMat.resize (count, count);
      init = true;
    }
    int pos2 = 0;
    while (iss >> c) {
      contactMat(pos1, pos2) = c;
      pos2++;
    }
    pos1++;
  }
  _infile.close ();
  
  return count;
}


void Writer::writeTree (std::string workDir, std::string fileName, std::vector<binary::TreeNode *> &nodeList)
{
  std::string path = workDir + fileName;
  std::cout << "path=" << path << "\n";
  _outfile.open (path);
  for (int i = 0; i < nodeList.size (); i++) {
    for (int j = nodeList[i]->_val[0]; j <= nodeList[i]->_val[1]; j++)
      _outfile << std::to_string (j + 1) << " ";
    _outfile << "\n";
  }
  _outfile.close ();
}


void Writer::writeTree (std::string workDir, std::string fileName, std::vector<multi::TreeNode *> &nodeList)
{
  std::string path = workDir + fileName;
  std::cout << "path=" << path << "\n";
  _outfile.open (path);
  for (int i = 0; i < nodeList.size (); i++) {
    for (int j = nodeList[i]->_val[0]; j <= nodeList[i]->_val[1]; j++)
      _outfile << std::to_string (j + 1) << " ";
    _outfile << "\n";
  }
  _outfile.close ();
}