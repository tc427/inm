#pragma once
#ifndef MIN_H
#define MIN_H

#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

#include "INM.hpp"

class Function {
public:
  Function(size_t dim_) :
    dim(dim_),
    grad(dim_, 0) {
  }

  double operator () (std::vector<double> const& params) {
    return calculate(params);
  }
  
  double calculate(std::vector<double> const& params) {
    assert(params.size() == dim);
    return calc_(params);
  }
  
  size_t ndim() const {
    return dim;
  }
  
  virtual std::vector<double> const& gradient(std::vector<double> const& params) {
    assert(params.size() == dim);
    return grad_(params);
  }
  
protected:
  size_t dim;
  std::vector<double> grad;
  
private:
  virtual double calc_(std::vector<double> const& params) = 0;
  virtual std::vector<double> const& grad_(std::vector<double> const& params) = 0;
};

class ENM : public Function {
public:
  ENM(Structure const& molstruct, double max = 10);
  
private:
  virtual double calc_(std::vector<double> const& params);
  virtual std::vector<double> const& grad_(std::vector<double> const& params);
  
  void gradInt(Interaction const& inter, std::vector<double> const& positions);

  double cutoff;
  std::vector<Interaction> interactions;
};

class INM : public Function {
public:
  INM(Structure const& molstruct, std::vector<ENM> nrgs);
  
  void setDelta(double newdelta);

private:
  virtual double calc_(std::vector<double> const& params);
  virtual std::vector<double> const& grad_(std::vector<double> const& params);

  double delta;
  double qdelta;

  std::vector<ENM> enms;
};




template<typename T>
double dist(T const * id1, T const * id2, size_t range) {
  T val = T();
  T temp = T();
  for(size_t i = 0; i < range; ++i) {
    temp = *(id1+i) - *(id2+i);
    val += temp*temp;
  }
  return sqrt(val);
}

template<typename T>
T normsq(std::vector<T> const& vec) {
  return std::inner_product(vec.begin(), vec.end(), vec.begin(), T());
}

template<typename T>
double norm(std::vector<T> const& vec) {
  return sqrt(normsq(vec));
}

template<typename T>
double grms(std::vector<T> const& vec) {
  return sqrt(normsq(vec) / vec.size());
}


#endif
