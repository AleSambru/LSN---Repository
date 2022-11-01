/**************************************/
/*Laboratorio di Simulazione Numerica */
/*Studentessa: Sambruna Alessia 	  */
/*n matricola: 942397                 */
/**************************************/
/********************************************************************************
Esercitazione 3 :
>> 1. By sampling directly the final asset price S(T) for a GBM(r,\sigma^2)
>> 2. By sampling the discretized GBM(r,\sigma^2) path of the asset price dividing $[0,T]$ in $10^2$ time intervals: $S(0)$, $S(t_1)$, $S(t_2)$, ... up to $S(t_{100}=T=1)$
Make four pictures for the estimation of the European call-option prices, $C[S(0),0]$ (direct and discretized),
and put-option prices, P[S(0),0] (direct and discretized), with their uncertainties with a large number of asset prices at time t=T,
 say M> 10^4, as a function of the number of blocks, N As usual, in computing statistical uncertainties, use data blocking.
*********************************************************************************/
#include <iostream>
#include <fstream>
#include "function.h"
#include "random.h"
#include<cmath>

using namespace std;

// block mean parameters
int M = 10000;
int n_blk = 100;
int L = int(M / n_blk);

// problem variables
double s0 = 100; // asset price at t = 0
double T = 1; // delivery time
int n_time = 100 ;
double K = 100; // strike price
double r = 0.1; // risk free interest rate
double sigma = 0.25; // volatility
int N = 10e3;
//Put - Call variables - direct and indirect and for the block mean
double C_N, P_N, C_N_step, P_N_step ;
vector<double> C_N_direct, P_N_direct, C_N_indirect, P_N_indirect;
vector<double> C_N_prog_direct, P_N_prog_direct, C_N_prog_indirect, P_N_prog_indirect;
vector<double> C_N2_prog_direct, P_N2_prog_direct, C_N2_prog_indirect, P_N2_prog_indirect;
vector<double> C_N_dir_ave, P_N_dir_ave, C_N_indir_ave, P_N_indir_ave ;

// functions
void Put_Call_direct(Random &);
void Put_Call_indirect(Random &rnd , int n);
void Accumulate() ;
void Sum_Prog();
double N_calc(double d);

//main
int main() {
    //Variables for random methods-> contained in random.h
    Random rnd;
    cout << "The followiong code uses a MonteCarlo simulation tu calculate the Put option and the Call option of the Vanilla options.\n"
            "_________________________________________________________________________________________________________________________\n";

    for (int i = 0; i < M; ++i) {
        cout << "n of iteration = " << i << "/" << M << "\r" << flush ;
        Put_Call_direct(rnd);
        Put_Call_indirect(rnd,n_time ) ;
    }
    cout << "Accumulation\n=================\n" ;
    Accumulate() ;
    Sum_Prog() ;


    rnd.SaveSeed();
    return 0;
}

// function definition
void Put_Call_direct(Random &rnd) {

    double Zi, Si, Ci, Pi;
    double sum_ci = 0.;
    double sum_pi = 0.;


    for (int i = 0; i < N; i++) {
        Zi = rnd.Gauss(0, 1); // 0 : mean; 1: sigma
        Si = s0 * exp((r - 0.5 * pow(sigma, 2)) * T + sigma * Zi * sqrt(T));
        Ci = exp(-r * T) * fmax(0, (Si - K));
        Pi = exp(-r*T) * fmax(0, K-Si);
        //Pi = + Ci - (Si - K );
        sum_ci += Ci;
        sum_pi += Pi;
    }

    C_N = sum_ci / N;
    P_N = sum_pi / N;

    C_N_direct.push_back(C_N);
    P_N_direct.push_back(P_N);

};


double N_calc(double d){
    return 0.5*(1 + erf(d/M_SQRT2));
}

//
void Put_Call_indirect(Random &rnd , int n) {

    double Zi, Si, Ci, Pi;
    double sum_ci = 0.;
    double sum_pi = 0.;

    // n = number of time steps
    double dt = double(T) / n;// t step
    double ti = 0.;
    for (int i = 0; i < N; i++) {
        double t = dt * i + ti;

        Si = s0 ;//* exp((r - 0.5 * sigma*sigma)) * t + sigma * Zi * sqrt(t));
        for (int j = 0; j < n; j++) {
            Zi = rnd.Gauss(0, 1);
            Si = Si * exp((r - 0.5 * sigma*sigma) * dt + sigma * Zi * sqrt(dt));
        }
        Ci = exp(-r * T) * fmax(0, (Si - K));
        Pi = exp(-r*T) * fmax(0, K-Si);
        //Pi =  Ci -Si + K * exp(-1*r*(dt)) ;
        sum_ci += Ci;
        sum_pi += Pi;
    }

    C_N_step = sum_ci / N;
    P_N_step = sum_pi / N;

    C_N_indirect.push_back(C_N_step) ;
    P_N_indirect.push_back(P_N_step) ;

} ;

void Accumulate() {
    //cout << "Accumulation\n" ;
    //cout << "-" << "\r" << flush ;
    for (int i = 0; i < n_blk; ++i) {
        double sum_C_direct = 0, sum_P_direct = 0., sum_C_indirect = 0, sum_P_indirect = 0.;

        for (int j = 0; j < L; ++j) {
            int k = j + i * L;
            sum_C_direct += C_N_direct.at(k);
            sum_P_direct += P_N_direct.at(k);
            sum_C_indirect += C_N_indirect.at(k);
            sum_P_indirect += P_N_indirect.at(k);
        }
        sum_C_direct = sum_C_direct / L;
        sum_P_direct = sum_P_direct / L;
        sum_C_indirect = sum_C_indirect / L;
        sum_P_indirect = sum_P_indirect / L;

        C_N_dir_ave.push_back(sum_C_direct);
        P_N_dir_ave.push_back(sum_P_direct);
        C_N_indir_ave.push_back(sum_C_indirect);
        P_N_indir_ave.push_back(sum_P_indirect);
        cout << "==============================================================================\n";
        cout << "Direct: Cn = " << C_N_dir_ave.at(i) << " ; Pn = " << P_N_dir_ave.at(i) << endl;
        cout << "Indirect: Cn = " << C_N_indir_ave.at(i) << " ; Pn = " << P_N_indir_ave.at(i) << endl;
        cout << "==============================================================================\n";
    }
    //cout << "End Accumulation\n" ;
};

void Sum_Prog() {
    cout << "Creazione Files \n" ;
    ofstream out_dir ("Put_Call_direct_prog.dat", ios::app) ;
    ofstream out_indir ("Put_Call_indirect_prog.dat", ios::app) ;
    for (int i = 0; i < n_blk; i++) {
        double sum_C_direct = 0., sum_P_direct = 0.;
        double sum2_C_direct = 0., sum2_P_direct = 0.;
        double sum_C_indirect = 0., sum_P_indirect = 0.;
        double sum2_C_indirect = 0., sum2_P_indirect = 0.;

        for (int j = 0; j < i + 1; j++) {
            sum_C_direct += C_N_dir_ave[j];
            sum_P_direct += P_N_dir_ave[j];
            sum2_C_direct+= C_N_dir_ave[j]*C_N_dir_ave[j] ;
            sum2_P_direct+= P_N_dir_ave[j]*P_N_dir_ave[j] ;

            sum_C_indirect += C_N_indir_ave[j];
            sum_P_indirect += P_N_indir_ave[j];
            sum2_C_indirect+= C_N_indir_ave[j]*C_N_indir_ave[j] ;
            sum2_P_indirect+= P_N_indir_ave[j]*P_N_indir_ave[j] ;
        };

        // update vector
        C_N_prog_direct.push_back(sum_C_direct / (i + 1));
        P_N_prog_direct.push_back(sum_P_direct / (i + 1));
        C_N2_prog_direct.push_back(sum2_C_direct / (i + 1));
        P_N2_prog_direct.push_back(sum2_P_direct / (i + 1));

        C_N_prog_indirect.push_back(sum_C_indirect / (i + 1));
        P_N_prog_indirect.push_back(sum_P_indirect / (i + 1));
        C_N2_prog_indirect.push_back(sum2_C_indirect / (i + 1));
        P_N2_prog_indirect.push_back(sum2_P_indirect / (i + 1));
        // output
        out_dir << i*L << " " << P_N_prog_direct.at(i) << " " << error(P_N_prog_direct,P_N2_prog_direct, i ) << " " << C_N_prog_direct.at(i) << " " << error(C_N_prog_direct,C_N2_prog_direct, i ) << endl;
        out_indir << i*L <<  " " << P_N_prog_indirect.at(i) << " " << error(P_N_prog_indirect,P_N2_prog_indirect, i ) <<" " << C_N_prog_indirect.at(i) << " " << error(C_N_prog_indirect,C_N2_prog_indirect, i ) << endl;
    }
    out_indir.close() ;
    out_dir.close() ;
};

