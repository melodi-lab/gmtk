#pragma once

#include <boost/date_time/posix_time/posix_time.hpp>
using namespace boost::posix_time;

#include <boost/random.hpp>

class Random {
  typedef boost::mt19937 RandBase;
  typedef boost::variate_generator<RandBase&, boost::uniform_real<double>> RandGen;
  typedef boost::variate_generator<RandBase&, boost::normal_distribution<double>> NormGen;

  RandBase _base;
  RandGen _rand;
  NormGen _norm;

public:
  Random(unsigned int seed = -1) 
    : _base((seed == -1) ? (unsigned int)microsec_clock::local_time().time_of_day().fractional_seconds() : seed),
    _rand(_base, boost::uniform_real<double>(0, 1.0)),
    _norm(_base, boost::normal_distribution<double>(0, 1.0))
  {
  }

  double Normal() {
    return _norm();
  }

  double Uniform() {
    return _rand();
  }

  int Rand(int numValues) {
    return (int)floor(_rand() * numValues);
  }
};

extern Random globalRandom;
