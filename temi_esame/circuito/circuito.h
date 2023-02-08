#ifndef _circuito_h_
#define _circuito_h_

#include "randomgen.h"

class Circuito {

public:
  Circuito(double);
  ~Circuito() { ; };

  void Esegui() {

    m_R1_misurato = m_rgen.Gaus(m_R1_input, m_errR);
    m_R2_misurato = m_rgen.Gaus(m_R2_input, m_errR);
    m_R3_misurato = m_rgen.Gaus(m_R3_input, m_errR);

		
  };

  void Analizza(){

		m_rin_misurato = (m_R1_misurato*m_R2_misurato+m_R2_misurato*m_R3_misurato+m_R3_misurato*m_R1_misurato)/(m_R1_misurato+m_R2_misurato);

  };

  double getrinmis() { return m_rin_misurato; };
  double getrininput () { return m_rin_input; };

private:
  // generatore di numeri casuali

  RandomGen m_rgen;

  // parametri dell'apparato sperimentale
  double m_errR;

  // valori delle quantita' misurabili :
  // input    : valori assunti come ipotesi nella simulazione
  // misurato : valore dopo la simulazione di misura
  double m_R1_input, m_R1_misurato;
  double m_R2_input, m_R2_misurato;
  double m_R3_input, m_R3_misurato;

  double m_rin_input, m_rin_misurato;


};

#endif
