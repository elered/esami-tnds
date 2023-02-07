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

  rungekutta myk;

  TApplication myApp("myApp",0,0);
  TGraph primo;
  TGraph secondo;

  TCanvas c;
  c.Divide(1,2);

  c.cd(1);

  while(vec[1]>=0){

    primo.AddPoint(t,vec[1]);
    vec = myk.Passo(t,vec,h1,R);
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

  while(vec[1]>=0){

    secondo.AddPoint(t,vec[1]);
    vec = myk.Passo(t,vec,h2,R);
    t+=h1;

  }

  secondo.SetTitle("moto del razzo con h = 0.01");
  secondo.GetXaxis()->SetTitle("x [m]");
  secondo.GetYaxis()->SetTitle("y [m]"); 
  secondo.Draw("ALP");

  double hnuovo = passoh(vec,h1,R);

  cout << "il nuovo passo di integrazione è: " << hnuovo << endl;

  vec[0] = 0;
  vec[1] = 0;
  vec[2] = 0;
  vec[3] = 0;
  t = 0;
  double tmpx = 0;
  double tmpy = 0;

  while(vec[1]>=0){

    if(t==1) {

      cout << "La x dopo la spinta è: " << vec[0] << endl;
      cout << "La y dopo la spinta è: " << vec[1] << endl;
      cout << "La vx dopo la spinta è: " << vec[2] << endl;
      cout << "La vy dopo la spinta è: " << vec[3] << endl;
    }

    tmpx = vec[0];
    tmpy = vec[1];
    vec = myk.Passo(t,vec,hnuovo,R);
    t+=hnuovo;

  }


  double gittata = (vec[0]*tmpy-tmpx*vec[1])/(tmpy-vec[1]);
  cout << "Gittata max: " << gittata << endl;

  myApp.Run();


return 0;

}