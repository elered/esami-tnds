#ifndef __Integral_h__
#define __Integral_h__

#include "funzgiuli.h"
#include <iostream>
#include <cfloat>

using namespace std;

//Base class : generic integrator

class Integral {

 public:
  
  Integral (double a, double b){
    checkInterval (a,b);
    m_prec = 0;
    m_h = 0; 
    m_sum = 0;
    m_sum_e = 0;
    m_sum_o = 0;
    m_integral =0;
    m_nstep = 0;
  };

  virtual double Integra(double prec, FunzioneBase& f) = 0;
  double GetNstep() {return m_nstep;}

 protected:

  void checkInterval (double a, double b){
    m_a = min(a,b);
    m_b = max(a,b);
    if(a>b){
      m_sign = -1;
    }
    else{
      m_sign = 1;
    }
  }

  unsigned int m_nstep;
  double m_prec;
  double m_a, m_b;
  double m_sum, m_sum_e, m_sum_o, m_integral, m_h;
  int m_sign;
};

class Trapezoide : public Integral {

 public:

  Trapezoide (double a, double b) : Integral (a,b){;};

  double Integra(double prec, FunzioneBase& f) override{

    double fa = f.Eval(m_a);
    double fb = f.Eval(m_b);
    m_nstep = 2;
    double In;
    double I2n;
    double fk = (1./2.)*fa+(1./2.)*fb;

    do{
      double m_h = (m_b-m_a)/m_nstep;
      for(int i=1; i<m_nstep; i+=2) { 
        fk += f.Eval(m_a+i*m_h);
      }
      In = I2n; //??
      I2n = m_sign*fk*m_h;
      m_nstep *= 2; //Perchè? Consegna?
    }while(prec<(4./3.)*fabs(I2n-In));

      return I2n;
  };

};

#endif