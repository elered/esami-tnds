#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "circuito.h"
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

int main() {

  TApplication myapp("myapp", 0, 0);
  int N = 1000;

  vector <double> r = {1, 3, 5, 7, 9};
  vector <TH1F*> histo(r.size());
  vector <double> stddev(r.size());
  

  TCanvas c;
  c.Divide(2,3);

	for(int k=0; k< r.size(); k++) {

		Circuito C(r[k]);
		string histoname = "Grafico Rin con sigma R = " + to_string(r[k]);
        histo[k] = new TH1F(histoname.c_str(), histoname.c_str(), 100, 0, 0);	
		
		for (int i = 0; i < N; i++) {
		
   		C.Esegui();
    	C.Analizza();
		histo[k]->Fill(C.getrinmis());
  	}
		c.cd(k+1);
    	histo[k]->GetXaxis()->SetTitle("valore Rin");
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

	TGraph invg(error.GetN(), error.GetY(), error.GetX());
	double x2 = invg.Eval(3);

	cout << "sigma R per avere l'errore su Rin di tre: " << x2 << endl;
	
  canvassino.cd();
  string title = "Errore in funzione di R";
  error.SetTitle(title.c_str());
  error.GetXaxis()->SetTitle("R [Ohm]");
  error.GetYaxis()->SetTitle("Errore");
  error.SetLineColor(6);
  error.SetMarkerColor(6);
  error.SetMarkerStyle(8);
  error.Draw("ALP");

  TF1 *f1 = new TF1("f1","[0]+[1]*pow(x,[2])",0.,0.);

  error.Fit(f1);

  myapp.Run();
  
  return 0;
}



