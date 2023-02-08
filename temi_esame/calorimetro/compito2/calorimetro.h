#ifndef _calorimetro_h_
#define _calorimetro_h_

#include "randomgen.h"

class EsperimentoCal {

 public :

  EsperimentoCal() ;
  ~EsperimentoCal() {;} ;

  void Esegui() {

    m_MA_misurato = m_rgen.Gaus(m_MA_input, err_MA);
    m_TE_misurato = m_rgen.Gaus(m_TE_input, err_TE);
    m_MSTAR_misurato = m_rgen.Gaus(m_MSTAR_input, err_MSTAR);
    m_TA_misurato = m_rgen.Gaus(m_TA_input, err_TA);
    m_TC_misurato = m_rgen.Gaus(m_TC_input, err_TC);

  } ;

  void Analizza() {

    m_cx_misurato = ((m_MA_misurato+m_MSTAR_misurato)*m_Ca*(m_TE_misurato-m_TA_misurato))/(m_mC*(m_TC_misurato-m_TE_misurato));

    m_cx_mis_MA = ((m_MA_misurato+m_MSTAR_input)*m_Ca*(m_TE_input-m_TA_input))/(m_mC*(m_TC_input-m_TE_input));
    m_cx_mis_MSTAR = ((m_MA_input+m_MSTAR_misurato)*m_Ca*(m_TE_input-m_TA_input))/(m_mC*(m_TC_input-m_TE_input));
    m_cx_mis_TA = ((m_MA_input+m_MSTAR_input)*m_Ca*(m_TE_input-m_TA_misurato))/(m_mC*(m_TC_input-m_TE_input));
    m_cx_mis_TC = ((m_MA_input+m_MSTAR_input)*m_Ca*(m_TE_input-m_TA_input))/(m_mC*(m_TC_misurato-m_TE_input));
    m_cx_mis_TE = ((m_MA_input+m_MSTAR_input)*m_Ca*(m_TE_misurato-m_TA_input))/(m_mC*(m_TC_input-m_TE_misurato));

  } ;

  double getmcxmisurato() { return m_cx_misurato ; } ;
  double getmcxinput() { return m_cx_input ; } ;

  double getmcxinput() { return m_cx_input ; } ;
  double getmcxinput() { return m_cx_input ; } ;
  double getmcxinput() { return m_cx_input ; } ;
  double getmcxinput() { return m_cx_input ; } ;
  double getmcxinput() { return m_cx_input ; } ;
                                                      
 private:                                       

  RandomGen m_rgen ;
                                                                                     
  double m_Ca, m_mC;

  double m_TA_input, m_TA_misurato;
  double m_MA_input, m_MA_misurato;
  double m_TC_input, m_TC_misurato;
  double m_TE_input, m_TE_misurato;
  double m_MSTAR_input, m_MSTAR_misurato;
  double m_cx_input, m_cx_misurato;
  double m_cx_mis_MA, m_cx_mis_TE, m_cx_mis_MSTAR, m_cx_mis_TA, m_cx_mis_TC;

  double err_MA, err_TE, err_MSTAR, err_TA, err_TC;


};

#endif

