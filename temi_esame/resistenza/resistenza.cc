#include "resistenza.h"
#include <cmath>

Resistenza::Resistenza(double r) :
  m_rgen(1),
  m_R_input(r),
	m_X_input(500),
  m_errR(5E-03*m_R_input),
  m_errl1(0.002),
  m_L(1)
{

	m_l1_input = m_L/((m_X_input/m_R_input)+1);
	m_l2_input = m_L-m_l1_input;

}  