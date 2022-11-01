/*******************************************************************************
Laboratorio di Simulazione Numerica
Studentessa: Sambruna Alessia 
n matricola: 942397
Esercitazione 1
Esercizio 3 : 
	Simulate the Buffon’s experiment:  A needle of length $L$ is thrown at random onto a horizontal plane ruled with straight lines a 
	distance d> L apart. The probability P that the needle will intersect one of these lines is: P = 2L/\pi d. 
	This could be used to evaluate pi from throws of the needle: if the needle is thrown down N_thr times and is observed to land on a
	line N_hit of those times, we can make an estimate of pi from
	pi = 2L/Pd = Nthr -> inf 2LNthr/(Nhit d)
*********************************************************************************/
#include <iostream>
#include <fstream>
#include "random.h"
#include "function.h"
#include<cmath>
#include<cstdlib>
#include<vector>

using namespace std;

//********** MAIN  ************
int main() {
    //Variables for random methods-> contained in random.h
    Random rnd;
    //Definition of variables
    double l = 1.;
    double d = 10.;
    // Probability P = 2L/\pi d = N_hit/N_throws
    //Calcolare tramite esperimento N_hit/Nthrows
    int N_th = 10000; // numero di lanci
    int M = 10000;
    vector<double> pi_values;
    for (int j = 0; j < M; j++)// ripeteremo la stima di pi, quindi ripeteremo l'esperimento M volte.
    {
        //--Experiment=> throws
        int N_hit = 0; // numero di hit, numero di volte che stecchino colpisce le d
        for (int i = 0; i < N_th; i++) {
            bool result = experiment(rnd, l, d);
            if (result == true){
                N_hit++;
            }
        }

        double P_measured; // probabilità misurata
        P_measured = double(N_hit) / double(N_th);
        double pi = 2 * l / (P_measured * d);
        pi_values.push_back(pi);
        cout << "________\nExperiment n° " << j << "; pi value = " << pi_values[j] <<"\r"<< flush;
    }

// mean with the block method
    int n_blk= 100;
    vector<double> ave = Accumulate(n_blk, pi_values) ;
    vector<double> err;

//calc error	
    vector<double> sum_prog = Sum_Prog(n_blk, ave);
    vector<double> su2_prog = Sum2_Prog(n_blk, ave);
    for (int i = 0; i < n_blk; i++) {
        err.push_back(calc_error(sum_prog, su2_prog, i));
    }

    // output
    ofstream out("Ex3_1_pi_values.dat", ios::app);
    if (!out) cerr << "Errore file output\n";
    for ( int i = 0; i<n_blk;i++) {
        cout << "Value " << i << " = " << sum_prog.at(i) << " Error = " << err.at(i)<<endl;
        out << i *l << " " << sum_prog.at(i) << " " << err[i] << endl;
    }

// print results
    cout << "------------Final Result--------------------\n" ;
    cout << "Pi_ best = " << Average_calculator(sum_prog) << endl;

// end
    rnd.SaveSeed();
    return 0;
}













