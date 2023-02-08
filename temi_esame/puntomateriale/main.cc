#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TColor.h"
#include "TApplication.h"
#include "TAxis.h"
#include "equazionidifferenziali.h"

using namespace std;

int main(){

  double x = 1;
  double y = 0;
  double vx = 0;
  double vy = -1;
  double T = 2*M_PI;

  double h = 0.1;
  double tmax = 10*T;
  double t = 0;
  double alpha = 0;

  vector<double> vec {x,y,vx,vy};

  rungekutta E;

  puntomat P (alpha);

  double hnuovo = passoh(tmax,vec,h,P);
  cout << "Il passo giusto è: " << hnuovo << endl;

  TApplication myApp("myApp",0,0);
  TCanvas c;
  c.Divide(3,2);

  c.cd(1);

  TGraph f;
  vec = {1.,0.,0.,-1.};
  t = 0;
  while(t<=tmax) {
    f.AddPoint(vec[0],vec[1]);
    vec = E.Passo(t,vec,hnuovo,P);
    t = t+hnuovo;
  }
	
  f.SetTitle("Moto con alpha = 0 e runge kutta");
  f.GetXaxis()->SetTitle("x [m]");
  f.GetYaxis()->SetTitle("y [m]"); 
  f.Draw("ALP");

  cout << "Errore percentuale di runge kutta: " << setprecision(2) << fabs(vec[0]-1.)*100 << " %" << endl;

  vec = {1.,0.,0.,-1.};
  t = 0;
  alpha = 1.;
  puntomat P2(alpha);
  Eulero E2;

  c.cd(2);

  TGraph f1;
  while(t<=tmax) {
    f1.AddPoint(vec[0],vec[1]);
    vec = E2.Passo(t,vec,hnuovo,P2);
    t = t+hnuovo;
  }
	
  f1.SetTitle("Moto con alpha = 1 e eulero");
  f1.GetXaxis()->SetTitle("x [m]");
  f1.GetYaxis()->SetTitle("y [m]"); 
  f1.Draw("ALP");

  c.cd(3);

  vec = {1.,0.,0.,-1.};
  t = 0;
  alpha = -1.;
  puntomat P3(alpha);

  TGraph f2;
  while(t<=tmax) {
    f2.AddPoint(vec[0],vec[1]);
    vec = E2.Passo(t,vec,hnuovo,P3);
    t = t+hnuovo;
  }
	
  f2.SetTitle("Moto con alpha = -1 e eulero");
  f2.GetXaxis()->SetTitle("x [m]");
  f2.GetYaxis()->SetTitle("y [m]"); 
  f2.Draw("ALP");

  vec = {1.1,0.,0.,-1.};
  t = 0;

  c.cd(4);

  TGraph f3;
  while(t<=tmax) {
    f3.AddPoint(vec[0],vec[1]);
    vec = E.Passo(t,vec,hnuovo,P2);
    t = t+hnuovo;
  }
	
  f3.SetTitle("Moto con alpha = 1 e rungekutta");
  f3.GetXaxis()->SetTitle("x [m]");
  f3.GetYaxis()->SetTitle("y [m]"); 
  f3.Draw("ALP");

  vec = {1.1,0.,0.,-1.};
  t = 0;

  c.cd(5);

  TGraph f4;
  while(t<=tmax) {
    f4.AddPoint(vec[0],vec[1]);
    vec = E.Passo(t,vec,hnuovo,P3);
    t = t+hnuovo;
  }
	
  f4.SetTitle("Moto con alpha = -1 e rungekutta");
  f4.GetXaxis()->SetTitle("x [m]");
  f4.GetYaxis()->SetTitle("y [m]"); 
  f4.Draw("ALP");
  myApp.Run();


return 0;

}