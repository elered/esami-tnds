#include "calorimetro.h"
#include <cmath>

EsperimentoCal::EsperimentoCal() :
  m_rgen(1),
  m_TA_input(16.1),
  m_TC_input(90.6),
  m_TE_input(17.2),
  m_MSTAR_input(25),
  m_MA_input(150),
  m_Ca(1),
  m_mC(27.737),
  err_MA(2),
  err_TE(0.2),
  err_MSTAR(5),
  err_TA(0.2),
  err_TC(0.4)
{

  m_cx_input =((m_MA_input+m_MSTAR_input)*m_Ca*(m_TE_input-m_TA_input))/(m_mC*(m_TC_input-m_TE_input));

  	
}  