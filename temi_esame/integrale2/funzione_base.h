#ifndef __funzione_base_h__
#define __funzione_base_h__
#include <cmath>

class funzionebase {

  public:

  virtual double eval (double) const = 0;

  virtual ~funzionebase () {}

};

class segno {

  public: 

  double eval(double x) const {
  if(x<0)
    return -1;
  else 
    return 1;
  }

};

class seno : public funzionebase {
    
    public:

    double eval(double x) const {

        return sin(x);
    }
};

class gaussianina : public funzionebase {

    public:

    double eval(double x) const {

        double tmp = - (pow(x,2))/2;

        return pow(M_E,tmp);
    }


};

class gauss2 : public funzionebase {

    public:

    double eval(double x, double y) const {

      double rho = 0;

      double numeratore = -(pow(x,2)-2*rho*x*y+pow(y,2));

      double denominatore = 2*(1-pow(rho,2));

      return pow(M_E,numeratore/denominatore);
      
    }


};

#endif 