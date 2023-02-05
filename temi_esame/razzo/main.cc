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
  double h1 = 0.1;
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

    primo.AddPoint(t,vec[1]);
    vec = myk.Passo(t,vec,h1,R);
    t+=h1;
  }

  cout << "condizioni al termine della spinta: "<< endl;
  cout <<  "x: " << vec[0] <<endl;
  cout << "y: " << vec[1] << endl;
  cout << "vx: " << vec[2] << endl;
  cout << "vy: " << vec[3] << endl;

  while(vec[1]>0){

    primo.AddPoint(t,vec[1]);
    vec = myk.Passo(t,vec,h1,R2);
    t+=h1;

  }

  primo.SetTitle("moto del razzo con h = 0.1");
  primo.GetXaxis()->SetTitle("x [m]");
  primo.GetYaxis()->SetTitle("y [m]"); 
  primo.Draw("ALP");

  c.cd(2);

  vec[0] = 0;
  vec[1] = 0;
  vec[2] = 0;
  vec[3] = 0;
  t = 0;

  while(t<=1) {

    secondo.AddPoint(vec[0],vec[1]);
    vec = myk.Passo(t,vec,h2,R);
    t+=h2;
  }

  cout << "condizioni al termine della spinta2: "<< endl;
  cout <<  "x2: " << vec[0] <<endl;
  cout << "y2: " << vec[1] << endl;
  cout << "vx2: " << vec[2] << endl;
  cout << "vy2: " << vec[3] << endl;

  while(vec[1]>0){

    secondo.AddPoint(vec[0],vec[1]);
    vec = myk.Passo(t,vec,h2,R2);
    t+=h2;

  }

  secondo.SetTitle("moto del razzo con h = 0.01");
  secondo.GetXaxis()->SetTitle("x [m]");
  secondo.GetYaxis()->SetTitle("y [m]"); 
  secondo.Draw("ALP");


  cout << "Gittata max: " << vec[0] << endl;

  myApp.Run();


return 0;

}