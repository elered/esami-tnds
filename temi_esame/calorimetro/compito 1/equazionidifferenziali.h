#ifndef __equazionidifferenziali_h__
#define __equazionidifferenziali_h__

#include "vectoroperations.h"
#include <vector>
#include <cmath>

using namespace std;

// ===============================================================
// Classe astratta, restituisce la derivata nel punto x al tempo t
// ===============================================================

class FunzioneVettorialeBase {

public:
  virtual vector<double> Eval(double t, const vector<double> & x) const = 0;
};


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


class calorimetro: public FunzioneVettorialeBase {

  public:

  calorimetro(double k1, double k2, double T10, double T20) {m_k1 = k1, m_k2 = k2, m_T10 = T10, m_T20 = T20;};

  vector<double> Eval(double t, const vector<double> &x) const override {
		
  	vector <double> v {(m_k1*m_T20+m_k2*m_T10)-x[0]*(m_k1+m_k2)};

  	return v;
  }

private:

  double m_k1, m_k2, m_T10, m_T20 ;

};

class calorimetroreale: public FunzioneVettorialeBase {

  public:

  calorimetroreale(double k1, double k2, double k3, double T10, double T20) {m_k1 = k1, m_k2 = k2, m_k3 = k3, m_T10 = T10, m_T20 = T20;};

	

  vector<double> Eval(double t, const vector<double> &x) const override {

		const double rho = (m_k1*m_k3);
		const double gamma =  0.5*(m_k1+m_k2+m_k3);
		
		
  	vector <double> v(2);
		v[0] =x[1];
		v[1] = -2*gamma*v[0]-rho*x[0];
		
  	return v;
  }

private:

  double m_k1, m_k2, m_k3, m_T10, m_T20 ;

};



#endif 


