#include <iostream>
#include "integrale.h"
#include "funzione_base.h"
#include <iomanip>

using namespace std;

int main(int argc, char * argv[] ) {

  double d = 1E-5;
  double lambda = 5.89E-7;
  double l = 1;

  double prec = 0.0001;

  seno f;
  trapezoidi myint(-d/2, d/2);


  double I = myint.integra(prec, f);

  cout << setw(5) << "Prec: " << setw(5) << prec << endl
       << setw(5) << "I: " << setw(5) << I << endl;
  
  return 0;

}