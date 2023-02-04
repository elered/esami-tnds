#include "visco.h"
#include <cmath>

EsperimentoVisco::EsperimentoVisco() :
  m_rgen(1),
  m_lambda1(579.1E-9),
  m_lambda2(404.7E-9),
  m_alpha(60.*M_PI/180.),
  m_sigmat(0.3E-3),
  m_A_input(2.7),
  m_B_input(60000E-18)
{
  // calcolo degli indici di rifrazione
  m_n1_input = sqrt(m_A_input + m_B_input / (m_lambda1*m_lambda1) );
  m_n2_input = sqrt(m_A_input + m_B_input / (m_lambda2*m_lambda2) );

  // theta0 e' arbitrario, scelgo M_PI/2.

  m_th0_input = M_PI/2.;

  // determino theta1 e theta2

  m_dm1_input = 2.*asin( m_n1_input * sin (0.5 * m_alpha) ) - m_alpha ;
  m_th1_input = m_th0_input + m_dm1_input ;
  m_dm2_input = 2.*asin( m_n2_input * sin (0.5 * m_alpha) ) - m_alpha ;
  m_th2_input = m_th0_input + m_dm2_input ;
  	
}  