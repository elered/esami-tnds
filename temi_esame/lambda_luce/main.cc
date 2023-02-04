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
#include "luce.h"

using namespace std;

int main(){

    TApplication app("app",0,0);
    vector <double> dev(5);
    vector <double> media(5);

    TCanvas c;
    c.Divide(2,3);

    int n = 10000;

    EsperimentoLuce p;

    c.cd(1);

    TH1F l1("l1","l1",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l1.Fill(p.getlambda1());
    }

    media[0] = l1.GetMean();
    dev[0] = l1.GetStdDev();
    l1.StatOverflows(kTRUE);

    l1.SetTitle("lambda ordine 1");
    l1.GetXaxis()->SetTitle("valori di lambda");
    l1.GetYaxis()->SetTitle("occorrenze");

    l1.Draw();

    c.cd(2);

    TH1F l2("l2","l2",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l2.Fill(p.getlambda2());
    }

    media[1] = l2.GetMean();
    dev[1] = l2.GetStdDev();
    l2.StatOverflows(kTRUE);

    l2.SetTitle("lambda ordine 2");
    l2.GetXaxis()->SetTitle("valori di lambda");
    l2.GetYaxis()->SetTitle("occorrenze");

    l2.Draw();

    c.cd(3);

    TH1F l3("l3","l3",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l3.Fill(p.getlambda3());
    }
  
    media[2] = l3.GetMean();
    dev[2] = l3.GetStdDev();
    l3.StatOverflows(kTRUE);

    l3.SetTitle("lambda ordine 3");
    l3.GetXaxis()->SetTitle("valori di lambda");
    l3.GetYaxis()->SetTitle("occorrenze");

    l3.Draw();

    c.cd(4);

    TH1F l4("l4","l4",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l4.Fill(p.getlambda4());
    }

    media[3] = l4.GetMean();
    dev[3] = l4.GetStdDev();
    l4.StatOverflows(kTRUE);

    l4.SetTitle("lambda ordine 4");
    l4.GetXaxis()->SetTitle("valori di lambda");
    l4.GetYaxis()->SetTitle("occorrenze");

    l4.Draw();

    c.cd(5);

    TH1F l5("l5","l5",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l5.Fill(p.getlambda5());
    }

    media[4] = l5.GetMean();
    dev[4] = l5.GetStdDev();
    l5.StatOverflows(kTRUE);

    l5.SetTitle("lambda ordine 5");
    l5.GetXaxis()->SetTitle("valori di lambda");
    l5.GetYaxis()->SetTitle("occorrenze");

    l5.Draw();

    c.cd(6);

    TGraph err;

    for ( int k = 0 ; k < 5 ; k++ ) {
    double x = k+1;
    double y =  dev[k];
    err.SetPoint(k, x, y);
    }

    err.SetTitle("Andamento dell'errore");
    err.GetXaxis()->SetTitle("m");
    err.GetYaxis()->SetTitle("stddev");

    err.SetMarkerStyle(20);

    err.Draw("ALP");

    double sum = 0;
    double sumpesi = 0;

    for(int i = 0; i<5; i++) {

      sum += (1/pow(dev[i],2))*media[i];
      sumpesi += (1/pow(dev[i],2));
    }

    double mediapon = sum/sumpesi;
    double errmpon = 1/sqrt(sumpesi);

    cout << "Media pesata: " << mediapon << endl;
    cout << "Errore media pesata: " << errmpon << endl;

    vector <double> dev2(5);
    vector <double> media2(5);

    TCanvas c2;
    c2.Divide(2,3);

    c2.cd(1);

    TH1F l21("l21","l21",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l21.Fill(p.get2lambda1());
    }

    media2[0] = l21.GetMean();
    dev2[0] = l21.GetStdDev();
    l21.StatOverflows(kTRUE);

    l21.SetTitle("lambda ordine 1");
    l21.GetXaxis()->SetTitle("valori di lambda");
    l21.GetYaxis()->SetTitle("occorrenze");

    l21.Draw();

    c2.cd(2);

    TH1F l22("l22","l22",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l22.Fill(p.get2lambda2());
    }

    media2[1] = l22.GetMean();
    dev2[1] = l22.GetStdDev();
    l22.StatOverflows(kTRUE);

    l22.SetTitle("lambda ordine 2");
    l22.GetXaxis()->SetTitle("valori di lambda");
    l22.GetYaxis()->SetTitle("occorrenze");

    l22.Draw();

    c2.cd(3);

    TH1F l23("l23","l23",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l23.Fill(p.get2lambda3());
    }
  
    media2[2] = l23.GetMean();
    dev2[2] = l23.GetStdDev();
    l23.StatOverflows(kTRUE);

    l23.SetTitle("lambda ordine 3");
    l23.GetXaxis()->SetTitle("valori di lambda");
    l23.GetYaxis()->SetTitle("occorrenze");

    l23.Draw();

    c2.cd(4);

    TH1F l24("l24","l24",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l24.Fill(p.get2lambda4());
    }

    media2[3] = l24.GetMean();
    dev2[3] = l24.GetStdDev();
    l24.StatOverflows(kTRUE);

    l24.SetTitle("lambda ordine 4");
    l24.GetXaxis()->SetTitle("valori di lambda");
    l24.GetYaxis()->SetTitle("occorrenze");

    l24.Draw();

    c2.cd(5);

    TH1F l25("l25","l25",100,0,0);

    for(int i = 0; i<n; i++) {

        p.Esegui();
        p.Analizza();
        l25.Fill(p.get2lambda5());
    }

    media2[4] = l25.GetMean();
    dev2[4] = l25.GetStdDev();
    l25.StatOverflows(kTRUE);

    l25.SetTitle("lambda ordine 5");
    l25.GetXaxis()->SetTitle("valori di lambda");
    l25.GetYaxis()->SetTitle("occorrenze");

    l25.Draw();

    c2.cd(6);

    TGraph err2;

    for ( int k = 0 ; k < 5 ; k++ ) {
    double x = k+1;
    double y =  dev2[k];
    err2.SetPoint(k, x, y);
    }

    err2.SetTitle("Andamento dell'errore");
    err2.GetXaxis()->SetTitle("m");
    err2.GetYaxis()->SetTitle("stddev");

    err2.SetMarkerStyle(20);

    err2.Draw("ALP");

    double sum2 = 0;
    double sumpesi2 = 0;

    for(int i = 0; i<5; i++) {

      sum2 += (1/pow(dev2[i],2))*media2[i];
      sumpesi2 += (1/pow(dev2[i],2));
    }

    double mediapon2 = sum2/sumpesi2;
    double errmpon2 = 1/sqrt(sumpesi2);

    cout << "Media pesata 2: " << mediapon2 << endl;
    cout << "Errore media pesata 2: " << errmpon2 << endl;

    app.Run();

    return 0;

}
