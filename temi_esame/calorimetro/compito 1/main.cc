#include <iostream>
#include <string>
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include <vector>
#include <algorithm>
#include "equazionidifferenziali.h"
using namespace std;

int main() {

	TApplication myApp("myApp",0,0);
	
	double tmax = 600;
	double k1 = 0.002;
	double k2 = 0.01;
	double k3 = 0.0004;
	double T10 = 289.37;
	double T20 = 323.15;
	double t = 0;
	double h = 0.001;

	calorimetro C (k1,k2,T10,T20);
	calorimetroreale CR (k1,k2,k3,T10,T20);

	vector<double> temp {T10};

	TGraph mygraph;
  TCanvas can;
	
	rungekutta kuttino;

	while(t<=tmax) {
		mygraph.AddPoint(t, temp[0]);
		temp = kuttino.Passo(t,temp,h,C);
		t = t+h;	
	}

	
	
	cout << "La temperatura raggiunta del sistema ideale è: " << temp[0] << endl;

	t = 0;
	double dT_0 = (k1*T20 + k2 * T10) - (k1 + k2) * T10;
	vector <double> temp2{T10, dT_0};
	
	TGraph mygraph2;
	
	while(t<=tmax) {
		mygraph2.AddPoint(t, temp2[0]);
		temp2 = kuttino.Passo(t,temp2,h,CR);
		t = t+h;	
	}

	cout << "La temperatura raggiunta del sistema reale è: " << temp2[0] << endl;

  can.cd();
  TLegend *legenda = new TLegend(0.7,0.7,0.9,0.9);
  legenda->AddEntry(&mygraph, "calorimetro ideale", "LP");
  legenda->AddEntry(&mygraph2, "calorimetro reale", "LP");
  string title = "Calorimetro (RK h = " + to_string(h) + ")";
  mygraph.SetTitle(title.c_str());
  mygraph.GetXaxis()->SetTitle("Tempo [s]");
  mygraph.GetYaxis()->SetTitle("Temperatura T [K]");
  mygraph.SetLineColor(9);
  mygraph.SetMarkerColor(9);
  mygraph.Draw("ALP");
  mygraph2.GetXaxis()->SetTitle("Tempo [s]");
  mygraph2.GetYaxis()->SetTitle("Temperatura T [K]");
  mygraph2.SetLineColor(6);
  mygraph2.SetMarkerColor(6);
  mygraph2.Draw("sameLP");
  legenda->Draw("sameLP");

  myApp.Run("myApp");

	return 0;
}