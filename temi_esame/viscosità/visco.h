#ifndef _Viscosita_h_
#define _Viscosita_h_

#include "randomgen.h"

class Viscosita {

 public :

  Viscosita();
  ~Viscosita(){;};
  
  //Sporco dati input con statistica (t,s,R)
  
  void Esegui(){
    m_s_misurato = m_rgen.Gaus(m_s_input,m_sigmas);
    m_t1_misurato = m_rgen.Gaus(m_t1_input,m_sigmat);
    m_t2_misurato = m_rgen.Gaus(m_t2_input,m_sigmat);
    m_R1_misurato = m_rgen.Gaus(m_R1_input,m_sigmaR);
    m_R2_misurato = m_rgen.Gaus(m_R2_input,m_sigmaR);
  }

  //Calcolo nuovi valori con dati sporcati

  void Analizza(){

    //Eta1 calcolato con errori su t,s,R
    m_eta1_misurato = (2*m_R1_misurato*m_R1_misurato*m_g*m_t1_misurato*(m_rho-m_rho0))/(9*m_s_misurato);

    //Eta1 calcolato con errori su t
    m_eta1A_misurato = (2*m_R1_input*m_R1_input*m_g*m_t1_misurato*(m_rho-m_rho0))/(9*m_s_input);

    //Eta1 calcolato con errori su s
    m_eta1B_misurato = (2*m_R1_input*m_R1_input*m_g*m_t1_input*(m_rho-m_rho0))/(9*m_s_misurato);

    //Eta1 calcolato con errori su R
    m_eta1C_misurato = (2*m_R1_misurato*m_R1_misurato*m_g*m_t1_input*(m_rho-m_rho0))/(9*m_s_input);
    

    //Eta2 calcolato con errori su t,s,R
    m_eta2_misurato = (2*m_R2_misurato*m_R2_misurato*m_g*m_t2_misurato*(m_rho-m_rho0))/(9*m_s_misurato);

    //Eta2 calcolato con errori su t
    m_eta2A_misurato = (2*m_R2_input*m_R2_input*m_g*m_t2_misurato*(m_rho-m_rho0))/(9*m_s_input);

    //Eta2 calcolato con errori su s
    m_eta2B_misurato = (2*m_R2_input*m_R2_input*m_g*m_t2_input*(m_rho-m_rho0))/(9*m_s_misurato);

    //Eta1 calcolato con errori su R
    m_eta2C_misurato = (2*m_R2_misurato*m_R2_misurato*m_g*m_t2_input*(m_rho-m_rho0))/(9*m_s_input);

    }
    
  //Get misurati 

  double geteta1mis() {return m_eta1_misurato;};
  double geteta2mis() {return m_eta2_misurato;};
  double geteta1Amis() {return m_eta1A_misurato;};
  double geteta2Amis() {return m_eta2A_misurato;};
  double geteta1Bmis() {return m_eta1B_misurato;};
  double geteta2Bmis() {return m_eta2B_misurato;};
  double geteta1Cmis() {return m_eta1C_misurato;};
  double geteta2Cmis() {return m_eta2C_misurato;};
  double getR1mis() {return m_R1_misurato;};
  double getR2mis() {return m_R2_misurato;};
  double getsmis() {return m_s_misurato;};
  double gett1mis() {return m_t1_misurato;};   
  double gett2mis() {return m_t2_misurato;};   
  
 private:
                                        
  RandomGen m_rgen; //Generatore di numeri casuali 
  
  //Parametri dell'apparato sperimentale  
  double m_sigmat, m_sigmas, m_sigmaR, m_g, m_eta, m_rho, m_rho0;

  //valori delle quantita' misurabili :              
	//input: valori assunti come ipotesi nella simulazione 
	//misurato: valore dopo la simulazione di misura                                                     
	double m_eta1_misurato, m_eta1A_misurato, m_eta1B_misurato, m_eta1C_misurato;
  double m_eta2_misurato, m_eta2A_misurato, m_eta2B_misurato, m_eta2C_misurato;
  double m_R1_input, m_R1_misurato;
  double m_R2_input, m_R2_misurato;
  double m_s_input, m_s_misurato;
  double m_t1_input, m_t1_misurato;
  double m_t2_input, m_t2_misurato;

};

#endif

