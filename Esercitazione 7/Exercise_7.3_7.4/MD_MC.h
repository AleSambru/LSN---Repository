/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid__
#define __fluid__

//Random numbers
#include "random.h"

#include <string>
using namespace std ;
int seed[4];
Random rnd;

// input parameter of the operation of the program
int code_mode, stop, mode;
int n_step_eq ;

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, it, ie, iw, ip;
double walker[m_props];
double vtail, ptail;

// g(r) distribution parameters
const int nbins = 100;
double bin_size, sd;

const int wd = 12 ;

// averages
double blk_av[m_props], blk_norm, accepted, attempted , bin_av[nbins];
double glob_av[m_props], glob_av2[m_props], glob_bin[nbins], glob_bin2[nbins];

double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp, stima_g[nbins];
double err_pot, err_press, err_kin, err_etot, err_temp, err_gdir, err_g[nbins];
double bins[nbins], r[nbins], DeltaV[nbins], g_norm[nbins];

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;

// simulation
int iNVET, nstep, nblk, restart;
double delta;


//pigreco
const double pi=3.1415927;
const double kB=8.31/(6.022*pow(10, 23));



//functions

void Program_SetUp (void);
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);
void Input(string input_file_name) ;
void SetUp(int mode) ;
void Input(string input_file_name) ;
void Averages(int iblk, int mode) ;
void Print_gr(int ) ;
#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
