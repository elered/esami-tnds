#include "Integral.h"
#include "funzioni.h"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPad.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

int main() {

  double esatto = (3. / 16.) * pow(M_E, 2); // Valore esatto dell'integrale
  double a = 0.;
  double b = sqrt(M_E);
  int passi = 10;
  int nmax = 1024;
  stringstream print;
  Funz f;
  Function g;

  // Calcolo integrale della funzione con Midpoint

  Midpoint IntMP(a, b);
  double I_MP = IntMP.Integra(passi, f);

  print << "Calcolo dell'integrale con il metodo del Midpoint" << endl << endl;
  print << endl << setw(20) << "Integrale" << setw(20) << "Errore" << endl << endl;
  for (int i = 1; i <= passi; i++) {
    int k = pow(2, i);
    double integrale = IntMP.Integra(k, f);
    double errore = fabs(esatto - integrale);
    print << setw(20) << integrale << setw(20) << errore << endl;
  }
  print << endl;

  print << "Passi = " << passi << endl << "Integrale finale = " << I_MP << endl;

  // Calcolo integrale della funzione con Midright

  Midright IntMR(a, b);
  double I_MR = IntMR.Integra(passi, f);

  print << endl << "Calcolo dell'integrale con il metodo del Midright" << endl << endl;
  print << endl << setw(20) << "Integrale" << setw(20) << "Errore" << endl << endl;
  for (int i = 1; i <= passi; i++) {
    int k = pow(2, i);
    double integrale = IntMR.Integra(k, f);
    double errore = fabs(esatto - integrale);
    print << setw(20) << integrale << setw(20) << errore << endl;
  }
  print << endl;

  print << "Passi = " << passi << endl << "Integrale finale = " << I_MR << endl;

  // Calcolo integrale della funzione con Monte Carlo

  int seed = 1;
  int punti = 16;
  int nvolte = 1000;
  IntegralMC IntMC(seed);
  double I_MC;

  for (int i = 0; i < nvolte; i++) {
    I_MC = IntMC.IntegraleAVE(a, b, f, punti);
  }
  double errore_MC = fabs(esatto - I_MC);

  print << endl << "Calcolo dell'integrale con il metodo della Media" << endl << endl;
  print << "Integrale finale = " << I_MC << endl << "Errore commesso = " << errore_MC << endl << endl;

  // Calcolo n punti/errore

  print << "Calcolo dell'integrale con il metodo del Midpoint a 16 punti" << endl << endl;
  for (int i = 1; i <= 16; i++) {
    I_MP = IntMP.Integra(i, f);
  }
  double errore_MP = fabs(esatto - I_MP);

  print << "Passi = 16" << endl << "Integrale finale = " << I_MP << endl << "Errore commesso = " << errore_MP << endl << endl;

  double k = sqrt(pow(errore_MC, 2) * punti);
  double n_vero = pow((k / errore_MP), 2);

  print << "Numero di punti necessari per ottenere la precisione di Midpoint con Media = " << n_vero << endl;

  // Calcolo integrale della funzione con Midpoint

  double c = 0.;
  double d = 2.;
  Midpoint IntMP2(c, d);
  double I_MP2 = IntMP2.Integra(passi, g);

  print << endl << "Calcolo dell'integrale con il metodo del Midpoint" << endl;
  print << endl << setw(20) << "Integrale" << setw(20) << "Errore" << endl << endl;
  for (int i = 1; i <= passi; i++) {
    int k = pow(2, i);
    double integrale = IntMP2.Integra(k, g);
    double errore = fabs(esatto - integrale);
    print << setw(20) << integrale << setw(20) << errore << endl;
  }
  print << endl;

  print << "Passi = " << passi << endl << "Integrale finale = " << I_MP2 << endl;

  PrintAll(print);

  return 0;
}