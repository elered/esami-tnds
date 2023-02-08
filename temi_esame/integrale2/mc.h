#ifndef __randomgen_h__
#define __randomgen_h__
#include <cmath>
#include <iostream>
#include "funzione_base.h"

using namespace std;

class RandomGen{

  public:
     RandomGen(unsigned int n) {n = m_seed, m_a = 1664525, m_c = 1013904223, m_m = 1<<31;}; // 1<<31 è due alla 31, operatore binario shift

     void seta(unsigned int a) {m_a = a;}
     void setc(unsigned int c) {m_c = c;}
     void setm(unsigned int m) {m_m = m;}

     double Rand() { //distribuiti uniformemente in 0,1

        m_seed = (m_a*m_seed+m_c)%m_m;
        return double(m_seed)/double(m_m);

     };
     double Unif(double xmin, double xmax) { //distribuiti uniformememte in a,b

        return xmin + (xmax -xmin) * Rand();

     };
     double Exp(double mean) {

        return -(1/mean)*log(1.-Rand());

     };
     double Gaus(double mean, double sigma) {
      
      double s = Rand();
      double t = Rand();
      double x = sqrt(-2.*log(1.-s))*cos(2*M_PI*t);

      return mean+x*sigma;

     };
     
     double GausAR(double mean, double sigma) {

      double s = Rand();

      double t = Rand();

      double xmin = mean - 5*sigma;

      double xmax = mean + 5*sigma;

      double x = xmin + (xmax-xmin)*s;

      double max = 1/(sigma*sqrt(2 * M_PI));

      double f = max * 1/pow(M_E,(pow((x-mean), 2)/(2 * pow(sigma, 2))));

      double y = max*t;

      while (y > f) {

         return GausAR(mean, sigma);

      }

      return x;

     };

  private:
     unsigned int m_a, m_c, m_m;
     unsigned int m_seed;
};


class integralMC {
    
    public:

    integralMC(unsigned int seed) : m_myrand(seed) {};
    ~integralMC() {;} ;

    double integralehom(double xmin, double xmax, double fmax, funzionebase &f, int punti) {

        m_punti = punti;

        double nhit = 0;

        double x = 0;

        double y = 0;

        for(int i=0; i<m_punti; i++) {

            x = m_myrand.Unif(xmin,xmax);

            y = m_myrand.Unif(0,fmax);

            if(y<f.eval(x)) {

            nhit ++;

            }
        }

        return (xmax-xmin)*fmax*(nhit/m_punti);

    };

    double integraleave(double xmin, double xmax, funzionebase &f, int punti) {

        m_punti = punti;

        double sum = 0;

        for(int i = 0; i<m_punti; i++) {
            
            double x = m_myrand.Unif(xmin,xmax);

            sum += f.eval(x);
        }

        return (xmax-xmin)*(sum/m_punti);
    };

    private:

    RandomGen m_myrand;
    unsigned int m_punti;
};

#endif