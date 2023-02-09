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


//CALORIMETRO IDEALE


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

//CALORIMETRO REALE (AVEVI LA DERIVATA PRIMA NEL MAIN)

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

//MOTO GRAVITAZIONALE ROSETTA

class MotoGravitazionalerosetta : public FunzioneVettorialeBase {

public:

  MotoGravitazionalerosetta(double G, double M) { 
    m_G = G;
    m_M = M;
  };

  vector<double> Eval(double t, const vector<double> &w) const override {

    double raggio = 6*m_M/pow(w[0]*w[0]+w[1]*w[1],2);

    vector<double> r {w[2],w[3],-m_G*m_M*w[0]/(pow(w[0]*w[0]+w[1]*w[1],3./2.))+raggio*w[0],-m_G*m_M*w[1]/(pow(w[0]*w[0]+w[1]*w[1],3./2.))+raggio*w[1]};

    return r;
  }

private:
  double m_G, m_M;
};

//MOTO GRAVITAZIONALE NORMALE

class MotoGravitazionale : public FunzioneVettorialeBase {

public:

  MotoGravitazionale(double G, double M) { 
    m_G = G;
    m_M = M;
  };

  vector<double> Eval(double t, const vector<double> &w) const override {

    vector<double> r {w[2],w[3],-m_G*m_M*w[0]/(pow(w[0]*w[0]+w[1]*w[1],3./2.)),-m_G*m_M*w[1]/(pow(w[0]*w[0]+w[1]*w[1],3./2.))};

    return r;
  }

private:
  double m_G, m_M;
};

//PENDOLO

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

//OSCILLATORE ARMONOCO SMORZATO FORZATO

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

//ESAME PROTONE

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

//ESAME RAZZO (CLASSE CON TEMPO SEPARATO MAGGIORE E MINORE DI UNO)

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

//ESAME PUNTO MATERIALE

class puntomat: public FunzioneVettorialeBase{
	public:
	puntomat(double alpha) {m_alpha = alpha;}

	vector<double> Eval(double t, const vector<double> &x) const override {
		double B = pow(x[0]* x[0] + x[1]* x[1], m_alpha);
		
		vector<double> v { x[2], x[3], B * x[3], -B * x[2]};
    return v;
		}

		private:
		double m_alpha;
};

//ERRORE RUNTIME RUNGEKUTTA (RICORDATI CHE HAI FATTO IL RELATIVO)

double Passo_prec(double t, const vector<double> &x, double prec, double h, const FunzioneVettorialeBase &f) {
    unsigned int N = 2;
    rungekutta myk;
    vector <double> Xn(x.size());
    vector<double> X2n = myk.Passo(t,x,h,f);

    do{
      Xn = X2n;
      N *= 2.;
      h = (t-0)/N;

      X2n = myk.Passo(t,x,h,f);

    } while (prec < 16./15. * Distanza(Xn, X2n)/(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3])));
    
    return h;

  }

//ALTRO ERRORE RUNGEKUTTA

double errorrunge(double t, double tmax, vector<double> &x, double h, const FunzioneVettorialeBase &f) {

    rungekutta myk;
    vector <double> In (x.size());
    vector <double> I2n (x.size());
    for(int i = 0; i<x.size(); i++) {
      I2n = x;
    }
    double dist = 0;

    while(t <= tmax){
      In = I2n;
      I2n = myk.Passo(t,I2n,h,f);
      t = t+h;
      dist = Distanza(In,I2n);
    }

    return dist;
    
};

//PASSO DA PRECISIONE FISSATA

  double PassoPrec(double t, const vector<double> &x, double prec, const FunzioneVettorialeBase &f) {

    double nstep = 2;
    vector <double> In (x.size());
    vector <double> I2n (x.size());
		double h = 0;
    I2n = x;

    do{
      In = I2n;
      nstep *= 2;
      h = (t-0)/nstep;
      I2n = Passo(t,x,h,f);
    }while(prec < distanza(In,I2n));
    //Con la norma trovo l'errore relativo

    return h;
    
  };



#endif 


