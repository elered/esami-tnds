#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TPad.h"

#include "visco.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

int main() {

  TApplication myApp("myApp",0,0);
  stringstream print;
  unsigned int ripetizioni = 1000;
  Viscosita V;

  double eta1, eta2, sigma_eta1, sigma_eta2;
  double eta1A, eta2A, sigma_eta1A, sigma_eta2A;
  double eta1B, eta2B, sigma_eta1B, sigma_eta2B;
  double eta1C, eta2C, sigma_eta1C, sigma_eta2C;

  // Istogramma di eta 1, R1 = 0.01, errori su t,s,R

  TH1F E1("Misura di eta 1", "Misura di eta con R = 0.01 m", 100, 0, 0);
  for (int i = 0; i < ripetizioni; i++) {
    V.Esegui();
    V.Analizza();
    E1.Fill(V.geteta1mis());
  }
  E1.StatOverflows(kTRUE);

  eta1 = E1.GetMean();
  sigma_eta1 = E1.GetStdDev();

  // Istogramma di eta 1, R1 = 0.01, errori su t

  TH1F E1A("Eta 1, errore su t","Eta 1, errore su t", 100,0,0);
  for (int i = 0; i < ripetizioni; i++) {
    V.Esegui();
    V.Analizza();
    E1A.Fill(V.geteta1Amis());
  }
  E1A.StatOverflows(kTRUE);

  eta1A = E1A.GetMean();
  sigma_eta1A = E1A.GetStdDev();

  // Istogramma di eta 1, R1 = 0.01, errori su s

  TH1F E1B("Eta 1, errore su s","Eta 1, errore su s", 100,0,0);
  for (int i = 0; i < ripetizioni; i++) {
    V.Esegui();
    V.Analizza();
    E1B.Fill(V.geteta1Bmis());
  }
  E1B.StatOverflows(kTRUE);

  eta1B = E1B.GetMean();
  sigma_eta1B = E1B.GetStdDev();

  // Istogramma di eta 1, R1 = 0.01, errori su R

  TH1F E1C("Eta 1, errore su R","Eta 1, errore su R", 100,0,0);
  for (int i = 0; i < ripetizioni; i++) {
    V.Esegui();
    V.Analizza();
    E1C.Fill(V.geteta1Cmis());
  }
  E1C.StatOverflows(kTRUE);

  eta1C = E1C.GetMean();
  sigma_eta1C = E1C.GetStdDev();

  // Istogramma di eta 2, R2 = 0.005, errori su t,s,R

  TH1F E2("Misura di eta 2", "Misura di eta con R = 0.005 m", 100,0,0);
  for (int i = 0; i < ripetizioni; i++) {
    V.Esegui();
    V.Analizza();
    E2.Fill(V.geteta2mis());
  }
  E2.StatOverflows(kTRUE);

  eta2 = E2.GetMean();
  sigma_eta2 = E2.GetStdDev();

  // Istogramma di eta 2, R2 = 0.005, errori su t

  TH1F E2A("Eta 2, errore su t", "Eta 2, errore su t", 100,0,0);
  for (int i = 0; i < ripetizioni; i++) {
    V.Esegui();
    V.Analizza();
    E2A.Fill(V.geteta2Amis());
  }
  E2A.StatOverflows(kTRUE);

  eta2A = E2A.GetMean();
  sigma_eta2A = E2A.GetStdDev();

  // Istogramma di eta 2, R2 = 0.005, errori su s

  TH1F E2B("Eta 2, errore su s", "Eta 2, errore su s", 100,0,0);
  for (int i = 0; i < ripetizioni; i++) {
    V.Esegui();
    V.Analizza();
    E2B.Fill(V.geteta2Bmis());
  }
  E2B.StatOverflows(kTRUE);

  eta2B = E2B.GetMean();
  sigma_eta2B = E2B.GetStdDev();

  // Istogramma di eta 2, R2 = 0.005, errori su R

  TH1F E2C("Eta 2, errore su R", "Eta 2, errore su R", 100,0,0);
  for (int i = 0; i < ripetizioni; i++) {
    V.Esegui();
    V.Analizza();
    E2C.Fill(V.geteta2Cmis());
  }
  E2C.StatOverflows(kTRUE);

  eta2C = E2C.GetMean();
  sigma_eta2C = E2C.GetStdDev();

  // Stampo i grafici su Canvas

  TCanvas c;
  c.Divide(2,2);
  c.cd(1);
  E1.GetXaxis()->SetTitle("eta1");
  E1.GetYaxis()->SetTitle("N");
  E1.Draw();
  c.cd(2);
  E1A.GetXaxis()->SetTitle("eta1");
  E1A.GetYaxis()->SetTitle("N");
  E1A.Draw();
  c.cd(3);
  E1B.GetXaxis()->SetTitle("eta1");
  E1B.GetYaxis()->SetTitle("N");
  E1B.Draw();
  c.cd(4);
  E1C.GetXaxis()->SetTitle("eta1");
  E1C.GetYaxis()->SetTitle("N");
  E1C.Draw();

  TCanvas d;
  d.Divide(2,2);
  d.cd(1);
  E2.GetXaxis()->SetTitle("eta2");
  E2.GetYaxis()->SetTitle("N");
  E2.Draw();
  d.cd(2);
  E2A.GetXaxis()->SetTitle("eta2");
  E2A.GetYaxis()->SetTitle("N");
  E2A.Draw();
  d.cd(3);
  E2B.GetXaxis()->SetTitle("eta2");
  E2B.GetYaxis()->SetTitle("N");
  E2B.Draw();
  d.cd(4);
  E2C.GetXaxis()->SetTitle("eta2");
  E2C.GetYaxis()->SetTitle("N");
  E2C.Draw();

  // Stampo i valori finali su terminale
  
  cout << "Valori di eta misurati con R = 0.01m" << endl << endl;
  cout << "Errore su t,s,R: " << eta1 << " +/- " << sigma_eta1 << endl;
  cout << "Errore su t: " << eta1A << " +/- " << sigma_eta1A << endl;
  cout << "Errore su s: " << eta1B << " +/- " << sigma_eta1B << endl;
  cout << "Errore su R: " << eta1C << " +/- " << sigma_eta1C << endl << endl;

  cout << "Valori di eta misurati con R = 0.005m" << endl << endl;
  cout << "Errore su t,s,R: " << eta2 << " +/- " << sigma_eta2 << endl;
  cout << "Errore su t: " << eta2A << " +/- " << sigma_eta2A << endl;
  cout << "Errore su s: " << eta2B << " +/- " << sigma_eta2B << endl;
  cout << "Errore su R: " << eta2C << " +/- " << sigma_eta2C << endl;

  myApp.Run();

  return 0;
}