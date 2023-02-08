#include <cmath>
#include <iostream>
#include <iomanip>

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"

#include "equazionidifferenziali.h"

using namespace std;

int main() {

  TApplication myApp("myApp",0,0);
  double h = 10000;
  double t = 0;
  double x = 147098074000;
  double y = 0;
  double vx = 0;
  double vy = 30287;
  double G = 6.6742E-11;
  double M = 1.98844E30;
  vector <double> w {x,y,vx,vy};
  MotoGravitazionale Grav(G,M);
  MotoGravitazionalerosetta Rosa(G,M);
  rungekutta RungeKutta;

  TGraph traiettoria;
  TGraph rosetta;

  do{
    w = RungeKutta.Passo(t,w,h,Grav);
    t = t+h;
  }while(w[1]>=0);

  double afelio = fabs(w[0]);

  double periodo = t*2;

  cout << "Periodo in giorni: " << periodo/86400 << endl;

  double t2 = 0;
  w = {147098074000,0,0,30287};
  double tmp = 0;

  do{
    tmp = w[1];
    traiettoria.AddPoint(w[1],w[0]);
    w = RungeKutta.Passo(t,w,h,Grav);
    t2 = t2+h;
  } while (t2<=periodo);

  w = {147098074000,0,0,30287};

  cout << "Afelio: " << afelio << endl;
  cout << "Rapporto tra perielio e afelio: " << x/afelio << endl;

  t2 = 0;
  periodo = periodo*100;

  do{
    rosetta.AddPoint(w[1],w[0]);
    w = RungeKutta.Passo(t,w,h,Rosa);
    t2 = t2+h;
  } while (t2<=periodo);
  
  TCanvas c;
  c.Divide(1,2);

  c.cd(1);
  traiettoria.SetTitle("Moto gravitazionale");
  traiettoria.GetXaxis()->SetTitle("x");
  traiettoria.GetYaxis()->SetTitle("y");
  traiettoria.SetLineColor(1);
  traiettoria.SetMarkerStyle(10);
  traiettoria.Draw("ALP");
  c.cd(2);
  rosetta.SetTitle("Moto rosetta");
  rosetta.GetXaxis()->SetTitle("x");
  rosetta.GetYaxis()->SetTitle("y");
  rosetta.SetLineColor(1);
  rosetta.SetMarkerStyle(10);
  rosetta.Draw("ALP");

  myApp.Run();


  return 0;
}