#pragma once
#ifndef GRAD_H
#define GRAD_H

#include <iostream>

#include "min.hpp"


class Grad {
public:
  Grad(Function & func, std::vector<double> pos, double conval_ = 1e-5);
  
  void fullMin(size_t iterlimit = 200);
  std::vector<double> const& iter();
  
  std::vector<double> const& getpos() const;
  std::vector<double> & setPos();
  
  std::vector<double> const& getDir();
  
  std::vector<double> const& getGrad() const;
  void updateGrad();
  
  virtual bool hasConverged(double val) const;
  virtual bool hasConverged() const;
  
  size_t itercount() const;
  
  void info(std::ostream & out = std::cout) const;
  
protected:
  virtual void setDir() = 0;
  void linesearchDot();
  
  Function & fn;
  
  double conval;
  
  size_t counter;

  std::vector<double> lastpos;
  std::vector<double> lastdir;
  std::vector<double> lastgrad;
};

class SD : public Grad {
public:
  SD(Function & func, std::vector<double> pos, double conval_ = 1e-5);
  
protected:
  virtual void setDir();
};

class CG : public Grad {
public:
  CG(Function & func, std::vector<double> pos, double conval_ = 1e-5);
  
protected:
  virtual void setDir();
  
  virtual void updateParams();
  
  double beta;
};

class CGPR : public CG {
public:
  CGPR(Function & func, std::vector<double> pos, double conval_ = 1e-5);
protected:
  virtual void updateParams();
};

void printpos(std::vector<double> pos, std::ostream & out = std::cout);

#endif
