#include "visco.h"
#include <cmath>

Viscosita::Viscosita() :
    m_rgen(1),
    m_g(9.81),
    m_rho(2700),
    m_rho0(1250),
    m_eta(0.83),
    m_R1_input(0.01),
    m_R2_input(0.005),
    m_s_input(0.60),
    m_sigmat(0.01),
    m_sigmas(0.001),
    m_sigmaR(0.0001)
    {
    
    //Determino il tempo impiegato a scendere
    m_t1_input = (9*m_eta*m_s_input)/(2*m_R1_input*m_R1_input*m_g*(m_rho-m_rho0));

    m_t2_input = (9*m_eta*m_s_input)/(2*m_R2_input*m_R2_input*m_g*(m_rho-m_rho0));
      
    }  
  