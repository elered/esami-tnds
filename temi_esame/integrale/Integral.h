#ifndef __Integral_h__
#define __Integral_h__

#include "FunzioneBase.h"
#include "RandomGen.h"
#include "funzioni.h"
#include <iostream>
#include <vector>
#include <cfloat>

using namespace std;

//Base class : generic integrator

class Integral {

 public:
  
  Integral (double a, double b){
    checkInterval (a,b);
    m_nstep = 0;
    m_h = 0; 
    m_sum = 0;
    m_integral =0;
  };

  virtual double Integra(unsigned int nstep, FunzioneBase& f) = 0;

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
  double m_a, m_b;
  double m_sum, m_integral, m_h;
  int m_sign;
};

//Derived class : Midpoint integrator

class Midpoint : public Integral {

 public:

  Midpoint (double a, double b) : Integral (a,b){;};

  double Integra(unsigned int nstep, FunzioneBase &f) override{
    
    if(nstep == 0){
      cout << "Error number of steps is zero" << endl;
      exit(-1);
    }

    m_nstep = nstep;
    m_h = (m_b - m_a)/m_nstep;

    m_sum = 0.;
    for(unsigned int i=0; i<m_nstep; i++){
      m_sum += f.Eval(m_a+(i+0.5)*m_h);
    }
    m_integral = m_sign*m_sum*m_h;
    return m_integral;
  };

};

//Derived class : Midright integrator

class Midright : public Integral {

 public:

  Midright (double a, double b) : Integral (a,b){;};

  double Integra(unsigned int nstep, FunzioneBase &f) override{
    
    if(nstep == 0){
      cout << "Error number of steps is zero" << endl;
      exit(-1);
    }

    m_nstep = nstep;
    m_h = (m_b - m_a)/m_nstep;

    m_sum = 0.;
    for(unsigned int i=0; i<m_nstep; i++){
      m_sum += f.Eval(m_a+i*m_h);
    }
    m_integral = m_sign*m_sum*m_h;
    return m_integral;
  };

};

//Derived class : IntegralMC integrator

class IntegralMC{

 public:

  IntegralMC (unsigned int seed): m_myrand(seed){};
  ~IntegralMC() {;};

double IntegraleAVE(double xmin, double xmax, FunzioneBase &f, int punti){

    vector<double> v;

    for(int i=0; i<punti; i++){
      double x = m_myrand.Unif(xmin,xmax);
      double y = f.Eval(x);
      v.push_back(y);
    }
  
    double med = media<double>(v);

    return (xmax-xmin)*med;
  }

  private:
  RandomGen m_myrand;

};

#endif