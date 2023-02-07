#ifndef __equazionidifferenziali_h__
#define __equazionidifferenziali_h__

#include "vectoroperations.h"
#include<vector>
#include <cmath>

using namespace std;

// ===============================================================
// Classe astratta, restituisce la derivata nel punto x al tempo t
// ===============================================================

class FunzioneVettorialeBase {

public:
  virtual vector<double> Eval(double t, const vector<double> & x) const = 0;
};

// ===============================================================
// Caso fisico concreto
// ===============================================================

class OscillatoreArmonico : public FunzioneVettorialeBase {

public:

  OscillatoreArmonico(double omega0) { m_omega0 = omega0; };

  vector<double> Eval(double t, const vector<double> &x) const override {

  vector<double> v { x[1], -pow(m_omega0, 2)*x[0] };

  return v;

  }

private:
  double m_omega0;
};

// ===============================================================
// Classe astratta per un integratore di equazioni differenziali
// ===============================================================

class EquazioneDifferenzialeBase {
public:
  virtual vector<double> Passo(double t, const vector<double> &x, double h, const FunzioneVettorialeBase &f) const = 0;
};

// ===============================================================
// Integratore concreto, metodo di Eulero
// ===============================================================

class Eulero : public EquazioneDifferenzialeBase {
public:

  vector<double> Passo(double t, const vector<double> &x, double h, const FunzioneVettorialeBase &f) const override {

    return x+f.Eval(t,x)*h;

  }
};

class rungekutta : public EquazioneDifferenzialeBase {
public:

  vector<double> Passo(double t, const vector<double> &x, double h, const FunzioneVettorialeBase &f) const override {

    vector<double> k1 = f.Eval(t,x);
    vector<double> k2 = f.Eval(t+(h/2.), x+(k1*(h/2.)));
    vector<double> k3 = f.Eval(t+(h/2.), x+(k2*(h/2.)));
    vector<double> k4 = f.Eval(t+h, x+(k3*h));

    return x+(k1+(2.*k2)+(2.*k3)+k4)*(h/6.);

  }
};

class pendolo : public FunzioneVettorialeBase {

  public:

  pendolo(double lunghezza) {m_lunghezza = lunghezza;};

  double g = 9.806;

  vector<double> Eval(double t, const vector<double> &x) const override {

  vector<double> v { x[1], -pow(sqrt(g/m_lunghezza), 2)*sin(x[0]) };

  return v;

  }

private:
  double m_lunghezza;
};

class smorzforz : public FunzioneVettorialeBase {

  public:

  smorzforz(double omega0, double alpha, double omega) {m_omega0 = omega0, m_alpha = alpha, m_omega = omega;};

  vector<double> Eval(double t, const vector<double> &x) const override {

  vector<double> v { x[1], -pow(m_omega0, 2)*x[0]-m_alpha*x[1]+sin(m_omega*t) };

  return v;

  }

private:

  double m_omega0, m_alpha, m_omega;

};

class Protone : public FunzioneVettorialeBase {

    public:

    Protone(double E0, double k, double w) {m_omega = w, m_E0 = E0, m_k = k;};

    vector<double> Eval(double t, const vector<double> &x) const override {

    double mp = 1.67*pow(10,-27);
    double qp = 1.6*pow(10,-19);

    vector<double> v { x[1], -((qp/mp)*m_E0*sin((m_k*x[0])-(m_omega*t))) };

    return v;

    }

    private:

    double m_omega, m_k, m_E0;

};

class Razzo : public FunzioneVettorialeBase {

    public:

    Razzo(int M, double theta0, int S) {m_M = M, m_S = S, m_theta0 = theta0;};

    vector<double> Eval(double t, const vector<double> &x) const override {

    double g = 9.80665;


    if(t<=1) {

      vector <double> v { x[2], x[3], (m_S/m_M)*cos(m_theta0), -g+(m_S/m_M)*sin(m_theta0) };
      return v;

    } else {

        vector <double> v { x[2], x[3], 0, -g };

        return v;

    }

    }

    private:

    double m_M, m_S, m_theta0;

};

double passoh(vector<double> &x, double h, Razzo R) {

  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  x[3] = 0;

  rungekutta myrk;

  double t = 0;
  double tmpx = 0;
  double tmpy = 0;

  while(x[1]>=0){
    tmpx = x[0];
    tmpy = x[1];
    x = myrk.Passo(t,x,h,R);
    t=t+h;
  }

  if ((h/tmpx>=1E-4) && (h/tmpy>=1E-4)) { 

  return h;

  }

  h = h/2; 

  return passoh(x,h,R);

}

#endif 
