#ifndef __FunzioneBase__
#define __FunzioneBase__

#include <cmath>
#include <iostream>

class FunzioneBase {

  public:

  virtual double Eval (double) const = 0;
  virtual ~FunzioneBase () {}

};


class gaussianina : public FunzioneBase {

    public:

    double eval(double x) const {

        double tmp = - (pow(x,2))/2;

        return pow(M_E,tmp);
    }


};
#endif