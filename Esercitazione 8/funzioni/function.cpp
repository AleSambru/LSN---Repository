/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *Laboratorio di Simulazione Numerica
 * Studentessa: Sambruna Alessia
 * n matricola: 942397
 * Esercitazione 8.1 : Determination of the ground state of the ground state of a quantum system using a variational method.
 * |==> function.cpp
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
#include <cmath>
#include <vector>
#include "function.h"
#include "random.h"
#include <fstream>

using namespace std;


double Psi(double x, double mu, double sigma) {
    // Funzione che calcola la funzione d'onda di test
    double a1 = (x - mu) / (sqrt(2.) * sigma);
    double a2 = (x + mu) / (sqrt(2.) * sigma);
    double psi = exp(-a1 * a1) + exp(-a2 * a2);
    return psi;
};

// functions - Hamilonian
double derivative(double x, double mu, double sigma) {
    double alpha = pow((x+mu)/sigma,2);
    double beta = pow((x-mu)/sigma,2);
    return -0.5*1./sigma/sigma * ( alpha*exp(-alpha/2) + beta*exp(-beta/2) - exp(-alpha/2) - exp(-beta/2) );
};

double V(double x) {
    return pow(x, 4) - 5*x*x/2;
};

double K(double x, double mu, double sigma ) {
    return derivative(x, mu, sigma)/*hbar*hbar*//(2*/*m*/Psi(x, mu, sigma));
};

double Eloc(double x, double mu, double sigma) {
    return V(x) + K(x, mu, sigma);
};

// Integral
// mu, sigma -> parameters
// N -> number of point to use to sample the integral




// Statistics
template<typename T> double error(vector <T> AV, vector <T> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt((AV2[n] - pow(AV[n], 2) )/n);
    };
};


vector <double> Accumulate(vector<double> &v, int n_blk, int L ) {

    vector <double> v_acc ;
    cout << "Accumulation\n" ;
    for (int i = 0; i < n_blk; ++i) {
        double sum = 0.;
        for (int j = 0; j < L; ++j) {
            int k = j + i * L;
            sum += v.at(k);
        }
        sum /= L;
        v_acc.push_back(sum);
        cout << "-";
    }
    cout << "End Accumulation\n" ;
    return v_acc ;
};

void Sum_Prog(vector <double> &v, string filename, int n_blk, int L) {
    cout << "Creazione Files \n" ;
    vector <double> sum_prog, sum_prog2 ;
    ofstream output  (filename, ios::app) ;
    for (int i = 0; i < n_blk; i++) {
        double sum = 0. , sum2 = 0.;


        for (int j = 0; j < i + 1; j++) {
            sum += v.at(j);
            sum2+= v[j]*v[j] ;
        };

        // update vector
        sum_prog.push_back(sum / (i + 1));
        sum_prog2.push_back(sum2 / (i + 1));
        // output
        output  << i << " " << v.at(i) << " " <<  sum_prog.at(i)  << " " << error(sum_prog,sum_prog2, i ) << endl;
    }
    output.close() ;
};






// implementation classes


