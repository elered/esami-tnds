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

class diffrazione : public funzionebase {

    public:

    double eval(double x) const {

        return
    }


};

#endif 