/*
 **********************************************************
 
 				Power Grid Simulator
 		(Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: parameter.h,v $
 *    Authors: Xin Huang 
 *    Functions: parameter definition
 *    $Date: 2020/03/17 $
 *
 */


#ifndef PARAMETER_H_
#define PARAMETER_H_

/******************************************************************************************************************
* t_nuc = tau*exp(E_v/k/T)*kT/omega/B*exp(-eta-e*Z*rho*L/(4kT)*j)ln((e*Z*rho*L/(4*omega)*j)/(sigma_T+psi+(e*Z*rho*L/(4*omega)*j)-sigma_n))

*** all convert to capitals
*******************************************************************************************************************/
#define H_M12 (3e-7)
#define H_M34 (6e-7)
#define H_M56 (1.2e-6)
#define W_M12 (3e-7)
#define W_M34 (6e-7)
#define W_M56 (1.2e-6)
#define E (1.6e-19)
#define L_UNIT (1e-6) // L_UNIT 1e-9 for ibmpg3
#define G_E_V (0.674*E)
#define G_E_D (0.65*E) //EV: 0.687
#define G_E_A (0.86*E)
#define G_K   (1.38e-23) // unit 'J/K'
#define G_T (373)
#define G_PI (3.14159265358)
#define G_R_CU  (1.28e-10) // copper atomic radius, unit 'm'
#define G_OMEGA  (1.66e-29)  // copper atomic volume
#define G_B  (1e11)  // unit 'Pa'
//#define G_E  (1.602e-19) // 'e', unit 'C'
#define G_Z  10 // eZ: effective charge of the migrating atoms
#define G_THERM_STRESS  (4e8)  // unit 'Pa' , 400MPa is critical value of thermal stress
#define G_CRITICAL_STRESS  (5e8) // unit 'Pa', 600MPs
#define G_RESIS_CU (3e-8) // unit 'omega*m'
#define G_GRAIN_RADIUS (1e-7) // unit 'm' for 100nm copper layer
#define G_X (E*G_Z*G_RESIS_CU/G_OMEGA)
#define G_RSTRESS (G_B*G_R_OVER_DELTA*exp(-1*G_E_V/G_K/G_NONSTRESS_TEMP)/9)
#define I_SCALE (0.4) // the amplitude of current source, special: 0.3
/***********************************************************************
* tau = (l^2/D0)*exp(E_D/k/T)
***********************************************************************/
#define G_D0 (7.56e-5) //unit 'm^2/s'

/***********************************************************************
* eta = f*omega*sigma_T/k/T
***********************************************************************/
#define G_F  (0.6)
#define G_MODIFY  (0.1)

/***********************************************************************
* psi = (B/3)*(R/delta)*exp(-E_v/k/T_ZS)
***********************************************************************/
#define G_NONSTRESS_TEMP  (623) // unit 'K'
#define G_R_OVER_DELTA  (1e3)
#define G_RESIDU_STRESS (G_B/3*G_R_OVER_DELTA*)

/***********************************************************************
* deltaR = rho_Ta*V*deltaTime/(3*segWidth*h_Ta) -
* rho_Cu*V*deltaTime/(segWidth*segHeight)
***********************************************************************/
#define G_RESIS_TA  (1.31e-7) // unit 'omega*m'
#define G_HEIGHT_Cu (9e-8)  // 90nm
#define G_HEIGHT_TA (2e-8) // barrier thickness 20nm

/***********************************************************************
* voltage drop related parameters
***********************************************************************/
#define G_VDD (1.8) // unit 'V'
#define G_VOLDROP_PERC (0.1) // for ground network, FAIL_PERC = delV/VDD

/***********************************************************************
* blech effect related parameters
***********************************************************************/
#define G_JL_CRIC 2e5 // calculated from (OMEGA*SIGMA_N/(E*Z*RESIS_CU)), identical to experimental results: (2000~10000 A/cm), unit 'A/m'

/***********************************************************************
* EM simulation time step
***********************************************************************/
#define TSTEP 5e6 // second

#endif /* PARAMETER_H_ */
