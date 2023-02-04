#ifndef _visco_h_
#define _visco_h_

#include "randomgen.h"

class EsperimentoVisco {

 public :

  EsperimentoVisco() ;
  ~EsperimentoVisco() {;} ;

  void Esegui() {

    m_th0_misurato = m_rgen.Gaus(m_th0_input, m_sigmat);
    m_th1_misurato = m_rgen.Gaus(m_th1_input, m_sigmat);
    m_th2_misurato = m_rgen.Gaus(m_th2_input, m_sigmat);

  } ;

  void Analizza() {

    m_dm1_misurato = m_th1_misurato - m_th0_misurato;
    m_dm2_misurato = m_th2_misurato - m_th0_misurato;

    m_n1_misurato = (sin((m_dm1_misurato+m_alpha)/2))/(sin(m_alpha/2));
    m_n2_misurato = (sin((m_dm2_misurato+m_alpha)/2))/(sin(m_alpha/2));
    
    m_A_misurato = ((m_lambda2*m_lambda2*m_n2_misurato*m_n2_misurato)-(m_lambda1*m_lambda1*m_n1_misurato*m_n1_misurato))/((m_lambda2*m_lambda2)-(m_lambda1*m_lambda1));
    m_B_misurato = ((m_n2_misurato*m_n2_misurato)-(m_n1_misurato*m_n1_misurato))/((1/(m_lambda2*m_lambda2))-(1/(m_lambda1*m_lambda1)));

  } ;

  double getAmis() { return m_A_misurato ; } ;
  double getBmis() { return m_B_misurato ; } ;
  double getn1mis() { return m_n1_misurato ; } ;
  double getn2mis() { return m_n2_misurato ; } ;
  double getdm1mis() { return m_dm1_misurato ; } ;
  double getdm2mis() { return m_dm2_misurato ; } ;
  double getth0mis() { return m_th0_misurato ; } ;
  double getth1mis() { return m_th1_misurato ; } ;
  double getth2mis() { return m_th2_misurato ; } ;
                                                           
 private:

  // generatore di numeri casuali                                           

  RandomGen m_rgen ;

  // parametri dell'apparato sperimentale                                                                                        
  double m_rho, m_rho0, m_alpha, m_sigmat;

  // valori delle quantita' misurabili :                                     
	// input    : valori assunti come ipotesi nella simulazione               
	// misurato : valore dopo la simulazione di misura                                                                                     
  double m_A_input, m_A_misurato;
  double m_B_input, m_B_misurato;
  double m_n1_input, m_n1_misurato;
  double m_n2_input, m_n2_misurato;
  double m_dm1_input, m_dm1_misurato;
  double m_dm2_input, m_dm2_misurato;
  double m_th0_input, m_th0_misurato;
  double m_th1_input, m_th1_misurato;
  double m_th2_input, m_th2_misurato;

};

#endif

