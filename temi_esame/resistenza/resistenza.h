#ifndef _resistenza_h_
#define _resistenza_h_

#include "randomgen.h"

class Resistenza {

public:
  Resistenza(double);
  ~Resistenza() { ; };

  void Esegui() {

    m_R_misurato = m_rgen.Gaus(m_R_input, m_errR);
    m_l1_misurato = m_rgen.Gaus(m_l1_input, m_errl1);
		
  };

  void Analizza(){

		m_X_misurato = m_R_misurato*(m_L-m_l1_misurato)/m_l1_misurato;

  };

  double getRmis() { return m_R_misurato; };
  double getXmis() { return m_X_misurato; };
  double getRinput() { return m_R_input; };
  double getl1mis() { return m_l1_misurato; };
  double getl1input() { return m_l1_input; };

private:
  // generatore di numeri casuali

  RandomGen m_rgen;

  // parametri dell'apparato sperimentale
  double m_errR, m_errl1, m_L;

  // valori delle quantita' misurabili :
  // input    : valori assunti come ipotesi nella simulazione
  // misurato : valore dopo la simulazione di misura
  double m_R_input, m_R_misurato;
  double m_X_input, m_X_misurato;
  double m_l1_input, m_l1_misurato;
  double m_l2_input, m_l2_misurato;

};

#endif
