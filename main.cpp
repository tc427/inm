#include <cmath>

#include "INM.hpp"
#include "grad.hpp"


void min(Structure & prot, ENM & enm) {
  std::ofstream outPDB;
  outPDB.open("bla.pdb", std::ofstream::out);

  CG minimizer(enm, prot.positions, 0.1);
  
  int mod = 0;
  
  while(true) {
    //minimizer.info();
    prot.positions = minimizer.getpos();
    prot.printPDB(outPDB, mod++);
    
    if(minimizer.hasConverged() || minimizer.itercount() > 100) {
      break;
    }
    
    minimizer.iter();    
  }
  minimizer.info();
  outPDB.close();
}

void interMin(Structure const& begin, Structure const& end) {
  Structure inter = begin;
  
  ENM enmbegin(begin, 10);
  ENM enmend(end, 10);
  
  inter.positions = begin.positions;
  
  INM inm(inter, {enmbegin, enmend});
  
  std::ofstream outPDB;
  outPDB.open("bla.pdb", std::ofstream::out);
  
  int mod = 0;
  
  for(double delta = 1; delta > 0; delta -= 0.01) {
    inm.setDelta(delta);
    std::cout << delta << std::endl;
    
    CGPR minimizer(inm, inter.positions, 0.1);
    minimizer.info();
    //std::cout << "mod : " << mod << std::endl;
  
  
    while(true) {
      //minimizer.info();
      
      if(minimizer.hasConverged() || minimizer.itercount() > 300) {
        break;
      }
      
      minimizer.iter();
    }
    
    inter.positions = minimizer.getpos();
    inter.printPDB(outPDB, mod++);
    minimizer.info();
  }
  outPDB.close();
}

void one(int argc, char **argv) {

  if(argc < 3) {
    std::cout << "Usage : " << argv[0] << " pdb1 pdb2" << std::endl;
    return;
  }

  std::ifstream beginfile;
  std::ifstream endfile;

  beginfile.open(argv[1], std::ifstream::in);
  endfile.open(argv[2], std::ifstream::in);


  Structure beginprot(beginfile);
  Structure endprot(endfile);

  beginfile.close();
  endfile.close();

  ENM enm(beginprot, 10);
  
  beginprot.positions[7] += 5;
  //beginprot.positions = endprot.positions;
  
  min(beginprot, enm);
}

void two(int argc, char **argv) {

  if(argc < 3) {
    std::cout << "Usage : " << argv[0] << " pdb1 pdb2" << std::endl;
    return;
  }

  std::ifstream beginfile;
  std::ifstream endfile;

  beginfile.open(argv[1], std::ifstream::in);
  endfile.open(argv[2], std::ifstream::in);


  Structure beginprot(beginfile);
  Structure endprot(endfile);

  beginfile.close();
  endfile.close();
  
  interMin(beginprot, endprot);
}

class TestFunc : public Function {
public:
  TestFunc() : Function(2) {
  }
  
protected:
  virtual double calc_(std::vector<double> const& params) {
    double p0 = params[0] - 4;
    double p1 = params[1] - 2;
    return 6*p0*p0 + 1*p1*p1;
  }
  
  virtual std::vector<double> const& grad_(std::vector<double> const& params) {
    grad = {12*(params[0] - 4), 2*(params[1] - 2) };
    return grad;
  }
};

void test() {
  TestFunc f;
  
  std::vector<double> init = {53, 25};
  
  CGPR min(f, init);
  
  min.fullMin(20);
}

int main(int argc, char **argv) {
 
 //one(argc, argv);
 two(argc, argv);
 //test();
 
 return 0; 
}
