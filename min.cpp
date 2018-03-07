#include <cmath>
#include <iostream>

#include "min.hpp"

ENM::ENM(Structure const& molstruct, double max) :
  Function(molstruct.positions.size()),
  cutoff(max) {

  for(size_t idx1 = 0; idx1 < dim; idx1 += 3) {
    for(size_t idx2 = idx1+3; idx2 < dim; idx2 += 3) {
      double dist0 = dist(&molstruct.positions[idx1], &molstruct.positions[idx2], 3);
      if(dist0 > max) {
        continue;
      }
      interactions.push_back(Interaction{idx1, idx2, dist0});
    }
  }
}

double ENM::calc_(std::vector<double> const& params) {
  return 0.0f;
}

void ENM::gradInt(Interaction const& inter, std::vector<double> const& positions) {
  const size_t a = inter.a;
  const size_t b = inter.b;
  
  //std::cout << "a : " << a << "\tb : " << b << "\t d0 : " << inter.d0 << std::endl;

  double coeff = (b-a < 4)? 10 : 1;

  double n [3] = {0, 0, 0};
  for(size_t i = 0; i < 3; ++i) {
    n[i] = (positions[a+i] - positions[b+i]);
  }
  double d1 = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  double dist = coeff*(d1 - inter.d0)/d1;
  
  //std::cout << "diff : " << d1 - inter.d0 << "\tcoeff : " << coeff << "\tk : " << dist << std::endl;

  for(size_t i = 0; i < 3; ++i) {
    grad[a+i] += n[i] * dist;
    grad[b+i] -= n[i] * dist;
    //std::cout << a+i << " : " << n[i] * dist << std::endl;
  }
}

std::vector<double> const& ENM::grad_(std::vector<double> const& params) {
  grad.assign(grad.size(), 0);
  for(size_t idx = 0; idx < interactions.size(); ++idx) {
    gradInt(interactions[idx], params);
 }

  //std::cout << "Energy gradient = " << norm(grad) << std::endl;
  
  return grad;
}

INM::INM(Structure const& molstruct, std::vector<ENM> nrgs) :
  Function(molstruct.positions.size()),
  delta(1),
  qdelta(0),
  enms(nrgs) {
}

double INM::calc_(std::vector<double> const& params) {
  return enms[0](params)*delta + enms[1](params)*qdelta;
}

std::vector<double> const& INM::grad_(std::vector<double> const& params) {
  this->grad = enms[0].gradient(params);
  std::vector<double> grad1(enms[1].gradient(params));
  
  for(size_t idx = 0; idx < grad.size(); ++idx) {
    grad[idx] *= delta;
    grad[idx] += grad1[idx]*qdelta;
  }
  return grad;
}

void INM::setDelta(double newdelta) {
  delta = newdelta;
  qdelta = 1 - delta;
}
