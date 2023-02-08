#include <iostream>
#include "integrale.h"
#include "funzione_base.h"
#include "mc.h"
#include <iomanip>

using namespace std;

int main(int argc, char * argv[] ) {

  double prec = 0.000001;

  gaussianina f;
  trapezoidi myint(-5, 5);

  double I = myint.integra(prec, f);

  cout << "Valore integrale: " << I << endl;

  cout << "Ha fatto questo numero di passi: " << myint.Getnstep() << endl;

  cout << "Errore: " << fabs(I-sqrt(2*M_PI)) << endl;

  integralMC myint2(1);

  double I2 = 0;

  I2 = myint2.integraleave(-5,5,f,100);

  cout << I2 << endl;




  
  return 0;

}