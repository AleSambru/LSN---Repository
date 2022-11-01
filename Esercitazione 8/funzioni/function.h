/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *Laboratorio di Simulazione Numerica
 * Studentessa: Sambruna Alessia
 * n matricola: 942397
 * Esercitazione 8.1 : Determination of the ground state of the ground state of a quantum system using a variational method.
 * |==> Headerfile function.h
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __function_h__
#define __function_h__

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "function.h"
#include "random.h"

using namespace std;

// parameters
//double hbar = 1 ;
//double m = 1 ;

//statistics
template<typename T> double error(vector <T> AV, vector <T> AV2, int n);
vector <double> Accumulate(vector<double> &v, int n_blk, int L );
void Sum_Prog(vector <double> &v, string filename, int n_blk, int L);

// Hamiltonian
double Psi(double x, double mu, double sigma);
double derivative(double x, double mu, double sigma);
double V(double x);
double K(double x, double mu, double sigma );
double Eloc(double x, double mu, double sigma);




#endif