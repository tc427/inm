#include <iostream>

#include "grad.hpp"
#include "min.hpp"

Grad::Grad(Function & func, std::vector<double> pos, double conval_) :
  fn(func),
  conval(conval_),
  counter(0),
  lastpos(pos),
  lastdir(pos.size(), 0),
  lastgrad(fn.gradient(pos)) {
}

void Grad::fullMin(size_t iterlimit) {
  while(true) {
    info();
    
    if(hasConverged() || itercount() >= iterlimit)
      break;
    
    iter();
  }
}

std::vector<double> const& Grad::iter() {
  ++counter;
  linesearchDot();
  return lastpos;
}

void Grad::linesearchDot() {
  static double at = 1e-5;
  
  getDir();
  
  double dirGE = std::inner_product(lastdir.begin(), lastdir.end(), lastgrad.begin(), 0.0f);
  
  for(size_t idx = 0; idx < lastpos.size(); ++idx) {
    lastpos[idx] += lastdir[idx]*at;
  }
  updateGrad();
  
  double dirGEt = std::inner_product(lastdir.begin(), lastdir.end(), lastgrad.begin(), 0.0f);
  double diffGE = dirGEt - dirGE;
  
  double a = -at * dirGE / diffGE;
    if(fabs(a) > 100) {
    a = 1e-3;
  }

  for(size_t idx = 0; idx < lastdir.size(); ++idx) {
    lastpos[idx] += lastdir[idx]*(a - at);
  }
}

std::vector<double> const& Grad::getDir() {
  setDir();
  double ndir = norm(lastdir);
  if(ndir != 1) {
    for(size_t idx = 0; idx < lastdir.size(); ++idx) {
      lastdir[idx] /= ndir;
    }
  }
  return lastdir;
}

std::vector<double> const& Grad::getpos() const {
  return lastpos;
}

std::vector<double> const& Grad::getGrad() const {
  return lastgrad;
}

void Grad::updateGrad() {
  lastgrad = fn.gradient(lastpos);
}
  
bool Grad::hasConverged() const {
  return hasConverged(conval);
}

bool Grad::hasConverged(double val) const {
  if(grms(getGrad()) < val) {
    return true;
  }
  return false;
}

size_t Grad::itercount() const {
  return counter;
}

void Grad::info(std::ostream & out) const {
  out << "grms = " << grms(getGrad()) << " i : " << itercount() <<  " converged : " << hasConverged() << std::endl;
}

SD::SD(Function & func, std::vector<double> pos, double conval_) :
  Grad(func, pos, conval_) {
}

void SD::setDir() {
  for(size_t idx = 0; idx < lastdir.size(); ++idx) {
    lastdir[idx] = -lastgrad[idx];
  }
}

CG::CG(Function & func, std::vector<double> pos, double conval_) :
  Grad(func, pos, conval_),
  beta(0.0) {
}

void CG::setDir() {
  //std::cout << "beta = " << beta << std::endl;  
  for(size_t idx = 0; idx < lastdir.size(); ++idx) {
    lastdir[idx] = -lastgrad[idx] + lastdir[idx]*beta;
  }
  
  updateParams();
}

void CG::updateParams() {
  double lastD2 = normsq(lastgrad);

  lastgrad = fn.gradient(lastpos);
  // Fletcher-Reeves CG
  beta = normsq(lastgrad) / lastD2;
  if(beta < 0 || beta > 1) {
    beta = 0;
  }
}

CGPR::CGPR(Function & func, std::vector<double> pos, double conval_) :
  CG(func, pos, conval_) {
}

void CGPR::updateParams() {
  double lastD2 = normsq(lastgrad);

  std::vector<double> Delta(lastgrad);

  lastgrad = fn.gradient(lastpos);
  // Polak-Ribiere CG
  for(size_t idx = 0; idx < lastdir.size(); ++idx) {
    Delta[idx] = lastgrad[idx] - Delta[idx];
  }
  beta = std::inner_product(lastgrad.begin(), lastgrad.end(), Delta.begin(), 0.0f) / lastD2;
  
  if(beta < 0 || beta > 1) {
    beta = 0;
  }
}

void printpos(std::vector<double> pos, std::ostream & out) {
  for(size_t idx = 0; idx < pos.size(); ++idx) {
    out << pos[idx++] << " ";
    out << pos[idx++] << " ";
    out << pos[idx] << std::endl;
  }
  out << std::endl;
}




/*
void linesearchSep(Function &func, std::vector<double> dir, std::vector<double> & pos) {
  
  double stop = 1e-4;
  double init = 0;
  double step = 1e-1;
  
  std::vector<double> grad(func.gradient(pos));
  double orig = std::inner_product(dir.begin(), dir.end(), grad.begin(), 0.0f);
  std::cout << "orig : " << orig << std::endl;
  
  if(orig > 0) {
    step *= -1;
  }
  
  double dirGE = orig;
  
  while(orig*dirGE > 0) {
    init += step;
    for(size_t idx = 0; idx < pos.size(); ++idx) {
      pos[idx] += dir[idx]*init;
    }
    
    grad = func.gradient(pos);
    dirGE = std::inner_product(dir.begin(), dir.end(), grad.begin(), 0.0f);
  }
  
  std::cout << "orig2 : " << dirGE << std::endl;
  
  orig = dirGE;
  
  init /= -2;
  
  int mid = 0;
  while(fabs(dirGE) > stop && fabs(init) > 1e-8) {
    for(size_t idx = 0; idx < pos.size(); ++idx) {
      pos[idx] += dir[idx]*init;
    }
    
    grad = func.gradient(pos);
    dirGE = std::inner_product(dir.begin(), dir.end(), grad.begin(), 0.0f);
    //std::cout << "mid : " << dirGE << std::endl;
    
    if(orig*dirGE < 0) {
      init /= -2;
    }
    else {
      init /= 2;
    }
    orig = dirGE;
    ++mid;
  }
  std::cout << "mid : " << dirGE << std::endl;
  std::cout << "mid iter : " << mid << std::endl;
}









void linesearchDot(Function &func, std::vector<double> dir, std::vector<double> & pos) {
  static double at = 1e-5;
  
  double ndir = norm(dir);
  //std::cout << "direction norm : " << ndir << std::endl;
  if(ndir != 1) {
    for(size_t idx = 0; idx < dir.size(); ++idx) {
      dir[idx] /= ndir;
    }
  }
  //std::cout << "direction norm : " << norm(dir) << std::endl;
  
  std::vector<double> grad(func.gradient(pos));
  //printpos(grad);
  std::vector<double> gradold(grad);
  double dirGE = std::inner_product(dir.begin(), dir.end(), grad.begin(), 0.0f);
  
  for(size_t idx = 0; idx < pos.size(); ++idx) {
    pos[idx] += dir[idx]*at;
  }
  
  grad = func.gradient(pos);
  double dirGEt = std::inner_product(dir.begin(), dir.end(), grad.begin(), 0.0f);
  double diffGE = dirGEt - dirGE;
  
  double a = -at * dirGE / diffGE;
    if(fabs(a) > 100) {
    a = 1e-3;
  }

  //std::cout << "at = " << at << "\tdirGE =  " << dirGE << "\tdirGEt " << dirGEt << "\tdiffGE = " << diffGE << std::endl;
  //std::cout << "a = " << a << std::endl;
  
  //std::cout << pos[0] << std::endl;
  for(size_t idx = 0; idx < dir.size(); ++idx) {
    //if(gradold[idx] > 1e-4) {
      //std::cout << "idx : " << idx << " dirGE : " << dirGE << " gradold : " << gradold[idx] << " grad : " << grad[idx] << std::endl;
      pos[idx] += dir[idx]*(a - at);
    //}
  }
  //std::cout << dir[0]*(a - at) << std::endl;
  //std::cout << pos[0] << std::endl;
  
  grad = func.gradient(pos);
  //std::cout << "dir.grad = " << std::inner_product(dir.begin(), dir.end(), grad.begin(), 0.0f) << std::endl;
}

*/
