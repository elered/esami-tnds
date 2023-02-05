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

  double x = 1;
  double y = 0;
  double vx = 0;
  double vy = -1;
  double T = 2*M_PI;

  double h = 0.003;
  double tmax = 10*T;
  double t = 0;
  double alpha = 0;

  vector<double> vec {x,y,vx,vy};

  rungekutta E;

  puntomat P (alpha);

  TApplication myApp("myApp",0,0);

  TGraph f;

  while(t<=tmax) {
    f.AddPoint(vec[0],vec[1]);
    vec = E.Passo(t,vec,h,P);
    t = t+h;
  }
	double temp = vec[0];
	double dep = vec[1];
  double dep5 = vec[2];
  double dep6 = vec[3];
	
  f.SetTitle("Moto");
  f.GetXaxis()->SetTitle("x [m]");
  f.GetYaxis()->SetTitle("y [m]"); 
  f.Draw("ALP");

	t = 0;
	vec[0] = 1;
  vec[1] = 0;
  vec[2] = 0;
  vec[3] = -1;
  while(t<=2*tmax) {
    vec = E.Passo(t,vec,h,P);
    t = t+h;
  }
	double tmp2 = vec[0];
	double dep2 = vec[1];
  double dep3 = vec[2];
  double dep4 = vec[3];

	double errore = (temp-tmp2)*(16./15.);
	double errore2 = (dep-dep2)*(16./15.);
  double errore3 = (dep5-dep3)*(16./15.);
  double errore4 = (dep6-dep4)*(16./15.);

  cout << (errore+errore2+errore3+errore4)/4 << endl;

  //double hnuovo = passoh(tmax,vec,h,P);
  //cout << "Il passo giusto è: " << hnuovo << endl;

  /*rungekutta myk;

  TApplication myApp("myApp",0,0);
  TGraph primo;
  TGraph secondo;

  TCanvas c;
  c.Divide(1,2);

  c.cd(1);

  while(t<=1) {

    vec = myk.Passo(t,vec,h1,R);
    t+=h1;
  }

  cout << "condizioni al termine della spinta: "<< endl;
  cout <<  "x: " << vec[0] <<endl;
  cout << "y: " << vec[1] << endl;
  cout << "vx: " << vec[2] << endl;
  cout << "vy: " << vec[3] << endl;

  while(vec[1]!=0){

    //primo.AddPoint(vec[0],vec[1]);
    vec = myk.Passo(t,vec,h1,R2);
    t+=h1;
    cout << vec[0] << endl;

  }

  primo.SetTitle("moto delrazzo");
  primo.GetXaxis()->SetTitle("x [s]");
  primo.GetYaxis()->SetTitle("y [m/s]"); 
  primo.Draw("ALP");*/

  myApp.Run();


return 0;

}