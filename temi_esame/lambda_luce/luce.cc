#include "luce.h"
#include <cmath>

EsperimentoLuce::EsperimentoLuce() :
  m_rgen(1),
  m_lambda(589),
  m_d(20000),
  m_sigmat(0.001),
  m_th0_input(0)
{

  m_thm1_input = m_th0_input+asin(m_lambda/m_d);
  m_thm2_input = m_th0_input+asin(2*m_lambda/m_d);
  m_thm3_input = m_th0_input+asin(3*m_lambda/m_d);
  m_thm4_input = m_th0_input+asin(4*m_lambda/m_d);
  m_thm5_input = m_th0_input+asin(5*m_lambda/m_d);
  	
}  