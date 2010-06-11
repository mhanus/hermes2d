#ifndef __HERMES2D_SEA_BREEZE_PARAMS_H
#define __HERMES2D_SEA_BREEZE_PARAMS_H

/* This file is generated by equations.py

   Just modify the python script if you want to change the quantities. Note
   that all the quantities in this file are related, so you *cannot* just
   change R or g (say) without also recalculating p, T and rho.
*/

// Everything is in SI units
// initial presure p in Pascals
#define p_z(z) {{ params.p_z }}
// initial temperature T in Kelvin
#define T_z(z) {{ params.T_z }}
// initial gas density
#define rho_z(z) {{ params.rho_z }}

// other physical constants
#define R {{ params.R }}           // Gas constant [J/(kg*K)]
#define g {{ params.g }}           // gravitational acceleration [m/s^2]
#define c_v {{ params.c_v }}       // specific heat capacity [J/(kg*K)]

// main characteristic constants
#define l_r {{ params.l_r }} // units: m
#define u_r {{ params.u_r }} // units: m/s
#define rho_r {{ params.rho_r }} // units: kg/m^3
// other characteristic constants
#define t_r (l_r/u_r)  // time
#define p_r (rho_r*u_r*u_r) // presure
#define E_r (rho_r*u_r*u_r) // energy
#define g_r (l_r/(t_r*t_r)) // gravitational constant

// other constants
#define kappa (1 + R/c_v)

#endif