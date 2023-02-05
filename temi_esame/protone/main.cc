#include <iostream>
#include <cmath>
#include <vector>
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TApplication.h"
#include "TAxis.h"
#include "equazionidifferenziali.h"

using namespace std;

int main(){

  double E0 = pow(10,7);
  double k = 10.;
  double w = pow(10,9);
  double x0 = 0.;
  double v0 = pow(10,8);
  double t = 0.;
  double tmax = pow(10,-8);
  double h = pow(10,-11);

  vector<double> x {x0,v0};

  Protone P (E0,k,w);

  rungekutta myk;



  TApplication myApp("myApp",0,0);
  TGraph retta;
  TGraph b;
  TGraph g;
  TGraph b2;
  TGraph g2;
  TGraph f;
  TCanvas c;
  c.Divide(2,3);

  c.cd(1);

  while(t<=tmax) {

    retta.AddPoint(t,x[1]);
    x = myk.Passo(t,x,h,P);
    t+=h;

  }

  retta.SetTitle("moto protone");
  retta.GetXaxis()->SetTitle("tempo [s]");
  retta.GetYaxis()->SetTitle("velocita [m/s]"); 
  retta.Draw("ALP");

  c.cd(2);

  x[0] = 0.;
  t = 0.;
  x[1] = 1.1*pow(10,8);
  tmax = pow(10,-7);

  while(t<=tmax) {
    b.AddPoint(t,x[1]);
    x = myk.Passo(t,x,h,P);
    t+=h;
  }

  double A = x[1]-pow(10,8);

  b.SetTitle("moto protone");
  b.GetXaxis()->SetTitle("tempo [s]");
  b.GetYaxis()->SetTitle("velocita [m/s]"); 
  b.Draw("ALP");

  c.cd(3);

  x[0] = 0.;
  t = 0.;
  x[1] = 1.1*pow(10,8);
  tmax = pow(10,-7);

  for(double i = pow(10,-11); i<pow(10,-10); i+=pow(10,-11)) {

    g.AddPoint(i,fabs(x[1]-A));
    x = myk.Passo(t,x,i,P);
    t = t+i;

  }

  g.SetTitle("andamento dell'errore");
  g.GetXaxis()->SetTitle("passo");
  g.GetYaxis()->SetTitle("errore"); 
  g.Draw("ALP");

  c.cd(4);

  x[0] = 0.;
  t = 0.;
  x[1] = pow(10,8);
  tmax = pow(10,-8);
  E0 = pow(-10,7);
  Protone P2 (E0,k,w);

  while(t<=tmax) {

    b2.AddPoint(t,x[1]);
    x = myk.Passo(t,x,h,P2);
    t+=h;

  }

  b2.SetTitle("moto protone campo elettrico negativo");
  b2.GetXaxis()->SetTitle("tempo [s]");
  b2.GetYaxis()->SetTitle("velocita [m/s]"); 
  b2.Draw("ALP");

  c.cd(5);

  x[0] = 0.;
  t = 0.;
  x[1] = 1.1*pow(10,8);
  tmax = pow(10,-7);

  while(t<=tmax) {
    g2.AddPoint(t,x[1]);
    x = myk.Passo(t,x,h,P2);
    t+=h;
  }

  g2.SetTitle("moto protone campo elettrico negativo v = 1.1*10^8");
  g2.GetXaxis()->SetTitle("tempo [s]");
  g2.GetYaxis()->SetTitle("velocita [m/s]"); 
  g2.Draw("ALP");

 c.cd(6);

  x[0] = 0.;
  t = 0.;
  x[1] = 0.9*pow(10,8);
  tmax = pow(10,-7);

  while(t<=tmax) {
    f.AddPoint(t,x[1]);
    x = myk.Passo(t,x,h,P2);
    t+=h;
  }

  f.SetTitle("moto protone campo elettrico negativo v = 0.9*10^8");
  f.GetXaxis()->SetTitle("tempo [s]");
  f.GetYaxis()->SetTitle("velocita [m/s]"); 
  f.Draw("ALP");

  myApp.Run();

return 0;

}