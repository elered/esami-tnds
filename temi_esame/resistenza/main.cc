#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "resistenza.h"
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

int main() {

  TApplication myapp("myapp", 0, 0);
  int N = 1000;

	vector <double> r = {200, 500, 900, 1500, 3000};
  vector <TH1F*> histo(r.size());
	vector <double> stddev(r.size());
  

  TCanvas c;
  c.Divide(2,3);

	for(int k=0; k< r.size(); k++) {

		Resistenza Ressy(r[k]);
		string histoname = "Grafico X con R = " + to_string(r[k]);
        histo[k] = new TH1F(histoname.c_str(), histoname.c_str(), 100, 0, 0);	
		
		for (int i = 0; i < N; i++) {
		
   		Ressy.Esegui();
    	Ressy.Analizza();
		histo[k]->Fill(Ressy.getXmis());
  	}
		c.cd(k+1);
    	histo[k]->GetXaxis()->SetTitle("valore integrale");
    	histo[k]->GetYaxis()->SetTitle("occorrenze");
   		histo[k]->SetLineColor(9);
    	histo[k]->Draw();
   		histo[k]->StatOverflows(kTRUE);
		stddev[k] = histo[k]->GetStdDev();
	}

	TGraph error;
	TCanvas canvassino;
	double x = 0;
	double y = 0;

	for(int i = 0; i<r.size(); i++) {
		x = r[i];
		y = stddev[i];
		error.AddPoint(x,y);
	}

	int indmin;
	double min=stddev[0];
	for(int i = 0; i<r.size(); i++) {
		if(stddev[i]<min) {
			min = stddev[i];
			indmin = i;
		}
	}

	cout << "il valore della resistenza che minimzza l'errore è: " << r[indmin] << endl;
	
  canvassino.cd();
  string title = "Errore in funzione di R";
  error.SetTitle(title.c_str());
  error.GetXaxis()->SetTitle("R [Ohm]");
  error.GetYaxis()->SetTitle("Errore");
  error.SetLineColor(6);
  error.SetMarkerColor(6);
  error.SetMarkerStyle(8);
  error.Draw("ALP");

  myapp.Run();
  
  return 0;
}



