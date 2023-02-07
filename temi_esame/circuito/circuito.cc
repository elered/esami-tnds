#include "circuito.h"
#include <cmath>

Circuito::Circuito(double r) :
  m_rgen(1),
  m_R1_input(90),
	m_R2_input(10),
  m_R3_input(50),
  m_errR(r)
{

	m_rin_input = (m_R1_input*m_R2_input+m_R2_input*m_R3_input+m_R3_input*m_R1_input)/(m_R1_input+m_R2_input);

}  