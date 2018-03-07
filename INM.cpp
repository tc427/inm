#include <cstdio>

#include <iostream>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "INM.hpp"

Structure::Structure(std::ifstream &source) {
  using boost::lexical_cast;

  if(positions.size() > 0) {
    return;
  }

  double x,y,z;
  std::stringstream posstr;

  std::string line;
  while(source.good()) {
    getline(source, line);
    if(line.substr(0, 4) != "ATOM" ||
        line.substr(13, 2) != "CA") {
      continue;
    }
    posstr.str(line.substr(31, 24));
 
    posstr >> x >> y >> z;

    positions.push_back(x);
    positions.push_back(y);
    positions.push_back(z);
    
    res3.push_back(line.substr(17, 3));
  }
}

void Structure::printPDB(std::ofstream & output, int modelno) const {
  output << "MODEL " << modelno << std::endl;
  
  char line [80];
  for(size_t idx = 0; idx < positions.size(); idx += 3) {
    int resid = static_cast<int>(idx/3);
    sprintf(line,
    "%-6s%5d %-4s%c%3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
    "ATOM",
    resid+1, // atom no
    "CA",// atom->rawName,
    ' ',// atom->alternate,
    res3[resid].c_str(),// residue (3 letters)
    'A',// protein chain
    resid+1, // residue no
    ' ',
    positions[idx],
    positions[idx+1],
    positions[idx+2],
    0.0, //atom->occupancy,
    0.0,//atom->bFactor,
    "C",//atom->rawElement,
    ""//atom->charge
    );
    output << line;
  
    //~ output << "ATOM      " << idx+1 << "  CA  ARG A   " << idx+1 << "     ";
    //~ for( ; idx % 3 != 0; ++idx ) {
      //~ output << setprecision(3) << setw(7) << right << positions[idx];
    //~ }
  }
  output << "TER" << std::endl;
}
