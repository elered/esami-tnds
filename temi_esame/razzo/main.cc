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

  int M = 500;
  double theta0 = M_PI/3;
  double x = 0;
  double y = 0;
  int S = 80000; //agisce per un secondo
  double vx = 0;
  double vy = 0;
  double h1 = 0.0001;
  double h2 = 0.01;
  double tmax = 100;
  double t = 0;

  vector<double> vec {x,y,vx,vy};

  Razzo R (M,theta0,S);

  Razzo2 R2 (M,theta0,S);

  rungekutta myk;

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
  primo.Draw("ALP");

  myApp.Run();


return 0;

}