#pragma once
#ifndef IN_HPP
#define IN_HPP

#include <fstream>
#include <string>
#include <vector>

struct Interaction {
  size_t a;
  size_t b;
  double d0;
};

class Structure {
public:
  Structure(std::ifstream &source);
  
  void printPDB(std::ofstream & output, int modelno = 0) const;

  std::vector<double> positions;
  
  std::vector<std::string> res3;
};

#endif
