#ifndef _RandomGen_h_
#define _RandomGen_h_

#include "FunzioneBase.h"
#include <cmath>
using namespace std;

class RandomGen {

public:
  RandomGen(unsigned int) {
    m_a = 1664525;
    m_c = 1013904223;
    m_m = 1 << 31; //=2^31; max intero rappresentabile è con 32 bit più uno per
                   //il segno
  }

  void SetA(unsigned int a) { m_a = a; }
  void SetC(unsigned int c) { m_c = c; }
  void SetM(unsigned int m) { m_m = m; }

  double Rand() {
    m_seed = (m_a * m_seed + m_c) % m_m;
    m_r = (double)m_seed / m_m;
    return m_r;
  }

  double Unif(double xmin, double xmax) {
    double y = xmin + (xmax - xmin) * Rand();
    return y;
  }

  double Exp(double lambda) {
    double y = -(1 / lambda) * log(1 - Rand());
    return y;
  }

  double Gauss(double mean, double sigma) {
    double s = Rand();
    double t = Rand();
    double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);
    return mean + x * sigma;
  }

  double GaussAR(double xmin, double xmax, double max, FunzioneBase &F) {

    double x, y;
    do {
      double s = Rand();
      double t = Rand();
      x = xmin + (xmax - xmin) * s;
      y = max * t;
    } while (y > F.Eval(x));

    return x;
  }

private:
  unsigned int m_a, m_c, m_m;
  double m_r;
  unsigned int m_seed;
};

#endif