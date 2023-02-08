#include "TApplication.h"
#include "TCanvas.h" 
#include "TH1F.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TAxis.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "randomgen.h"
#include "calorimetro.h"

using namespace std;

int main(){

    vector <double> std (5);

    TApplication app("app",0,0);

    TCanvas c1;

    int n = 1000;

    EsperimentoCal c;

    c1.cd();

    TH1F cal("cal","cal",100,0,0);

    for(int i = 0; i<n; i++) {

        c.Esegui();
        c.Analizza();
        cal.Fill(c.getmcxmisurato());
    }

    cal.StatOverflows(kTRUE);
    double media = cal.GetMean();
    double dev = cal.GetStdDev();


    cal.SetTitle("Calorimetro");
    cal.GetXaxis()->SetTitle("valori di cx");
    cal.GetYaxis()->SetTitle("occorrenze");

    cal.Draw();

  
    cout << "Media: " << media << endl;
    cout << "Errore: " << dev << endl;

    double valorevero = c.getmcxinput();

    cout << "Valore vero: " << valorevero << endl;

    app.Run();


    return 0;

}
