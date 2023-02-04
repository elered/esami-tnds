#ifndef _visco_h_
#define _visco_h_

#include "randomgen.h"

class EsperimentoLuce {

 public :

  EsperimentoLuce() ;
  ~EsperimentoLuce() {;} ;

  void Esegui() {

    m_th0_misurato = m_rgen.Gaus(m_th0_input, m_sigmat);
    m_thm1_misurato = m_rgen.Gaus(m_thm1_input, m_sigmat);
    m_thm2_misurato = m_rgen.Gaus(m_thm2_input, m_sigmat);
    m_thm3_misurato = m_rgen.Gaus(m_thm3_input, m_sigmat);
    m_thm4_misurato = m_rgen.Gaus(m_thm4_input, m_sigmat);
    m_thm5_misurato = m_rgen.Gaus(m_thm5_input, m_sigmat);

  } ;

  void Analizza() {

    m_lambda1 = sin(m_thm1_misurato-m_th0_misurato)*m_d;
    m_lambda2 = sin(m_thm2_misurato-m_th0_misurato)*(m_d/2);
    m_lambda3 = sin(m_thm3_misurato-m_th0_misurato)*(m_d/3);
    m_lambda4 = sin(m_thm4_misurato-m_th0_misurato)*(m_d/4);
    m_lambda5 = sin(m_thm5_misurato-m_th0_misurato)*(m_d/5);

    m_2lambda1 = sin(m_thm1_misurato-m_th0_input)*m_d;
    m_2lambda2 = sin(m_thm2_misurato-m_th0_input)*(m_d/2);
    m_2lambda3 = sin(m_thm3_misurato-m_th0_input)*(m_d/3);
    m_2lambda4 = sin(m_thm4_misurato-m_th0_input)*(m_d/4);
    m_2lambda5 = sin(m_thm5_misurato-m_th0_input)*(m_d/5);


  } ;

  double getlambda1() { return m_lambda1 ; } ;
  double getlambda2() { return m_lambda2 ; } ;
  double getlambda3() { return m_lambda3 ; } ;
  double getlambda4() { return m_lambda4 ; } ;
  double getlambda5() { return m_lambda5 ; } ;

  double get2lambda1() { return m_2lambda1 ; } ;
  double get2lambda2() { return m_2lambda2 ; } ;
  double get2lambda3() { return m_2lambda3 ; } ;
  double get2lambda4() { return m_2lambda4 ; } ;
  double get2lambda5() { return m_2lambda5 ; } ;
                                                      
 private:                                       

  RandomGen m_rgen ;
                                                                                     
  double m_sigmat, m_d, m_lambda;

  double m_th0_input, m_th0_misurato;
  double m_thm1_input, m_thm1_misurato;
  double m_thm2_input, m_thm2_misurato;
  double m_thm3_input, m_thm3_misurato;
  double m_thm4_input, m_thm4_misurato;
  double m_thm5_input, m_thm5_misurato;

  double m_lambda1, m_lambda2, m_lambda3, m_lambda4, m_lambda5;
  double m_2lambda1, m_2lambda2, m_2lambda3, m_2lambda4, m_2lambda5;


};

#endif

