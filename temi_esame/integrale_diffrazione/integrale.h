#ifndef __integrali_h__
#define __integrali_h__

#include <iostream>
#include <cmath>
#include "funzione_base.h"

using namespace std;

// Base class : generic integrator

class integral {

 public:
  
  integral (double a, double b){
    checkinterval (a,b);
    m_nstep = 0;
    m_h = 0; 
    m_sum = 0;
    m_integral =0;
  }

  virtual double integra(unsigned int nstep, funzionebase& f) {return 0;};
  virtual double integra(double prec, funzionebase& f) {return 0;};

 protected:

  void checkinterval (double a, double b){

    m_a = min(a,b);
    m_b = max(a,b);

    if(a>b) {
      m_sign = -1;
    } else {
      m_sign =1;
    };

  }

  unsigned int m_nstep;
  double m_a, m_b;
  double m_sum, m_integral, m_h;
  int m_sign;
};

// derived class : Midpoint integrator

class midpoint : public integral {

 public:

  midpoint (double a, double b) : integral (a, b){;};

  double integra(unsigned int nstep, funzionebase& f) override {

    m_nstep = nstep;

    m_sum = 0;

    double m_h = (m_b-m_a)/m_nstep;

    for (int k = 0; k<(m_nstep); k++) {

    double xk = m_a + (k+0.5)*m_h;

    m_sum += f.eval(xk);

    }

    m_integral = m_sign*m_sum*m_h;

    return m_integral;

  };

};

class simpson : public integral {

 public:

  simpson (double a, double b) : integral (a, b){;};

  double integra(unsigned int nstep, funzionebase& f) override{

    m_nstep = nstep;

    m_sum = 0;

    if((m_nstep%2)!=0) {

        cout << "Inserire n pari" << endl;
        exit(-1);

    }

    double m_h = (m_b-m_a)/m_nstep;

    for (int k = 1; k<m_nstep; k++) {

    double xk = m_a + k*m_h;

    double fxk = f.eval(xk);

    if(k%2==0) {

    m_sum += 2*fxk;

    } else {

    m_sum += 4*fxk;

    }

 }

  m_integral = m_sign*((m_h/3)*(f.eval(m_a) + m_sum + f.eval(m_b)));

  return m_integral;



  }

};

class trapezoidi : public integral {

 public:

  trapezoidi (double a, double b) : integral (a, b){;};

  double integra(double prec, funzionebase& f) override{

    double fa = f.eval(m_a);
    double fb = f.eval(m_b);
    m_nstep = 2;
    double fk = 0.5*fa + 0.5*fb;

        do {

            m_h = (m_b-m_a)/m_nstep;

            for(int k=1; k<m_nstep; k+=2) { 

            fk += f.eval(m_a+k*m_h);

            };

            m_sum = m_integral;

            m_integral = m_sign*fk*m_h;

            m_nstep *= 2;
    

        } while(prec<(4./3.)*fabs(m_integral-m_sum));

        return m_integral;
 

    }

};

#endif 