/*******************************************************************************
Laboratorio di Simulazione Numerica
Studentessa: Sambruna Alessia
n matricola: 942397
Esercitazione 8.1 : Determination of the ground state of the ground state of a
 quantum system using a variational method.
*********************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;
#define _USE_MATH_DEFINES

ofstream out;
ofstream out_par;

double accepted_ann, accepted_met, attempted_ann, attempted_met;
double muold, sigmaold;
double m = 1, hbar = 1;
double beta, delta_ann;
int eq_step;
double mu = 0.5;
double sigma = 0.5;

#define PREC 0.0015
#define delta 2

//funzioni
double Psi_test(double x);

double V(double x);

double derivative(double x);

double K(double x);

double Eloc(double x);

double Metropolis(Random &rnd, double x);

double Estrazione(Random &rnd, double &x);

double Integral(Random &rnd, int N);

double update_param(Random &rnd, double beta);

template<typename T>
double error(vector <T> AV, vector <T> AV2, int n);

vector<double> Accumulate(vector<double> &v, int n_blk, int L);

void Sum_Prog(vector<double> &v, string filename, int n_blk, int L);

double Integral(Random &rnd, int M, int L, string filename);

using namespace std;

//=====================
// Main
//=====================
int main(int argc, char *argv[]) {

    Random rnd;
    double H;
    int L = 15; // Numero di punti estratti in ogni blocco (a temperatura costante)
    vector<double> mus, sigmas, betas;
    double betamin = 1, dbeta = 2.5;
    double err = 0.5;
    int k = 0, iter = 0;
    double error = PREC, beta;
    double aveH = 0, ave2H = 0;
    double ave_mu, ave2_mu, err_mu;//vettori per blocchi annealing
    double ave_sigma, ave2_sigma, err_sigma;//vettori per blocchi

    // Part 1 and 2 : Show a picture of H (with statistical uncertainties) as a function of the SA steps of the algorithm
    out.open("Part_1_H_SA_steps.dat");
    out_par.open("Part_2_mu_sigma.dat");
    while (err >= PREC) {
        iter++;
        betas.push_back(betamin + k * dbeta); // update temperature
        beta = betas.at(k);
        delta_ann = 1. / sqrt(beta);
        double sum_H = 0, sum2_H = 0, sum_mu = 0, sum_sigma = 0;

        //Ciclo nel blocco a temperatura costante
        for (int j = 0; j < L; j++) {
            H = update_param(rnd, beta);
            sum_H += H;
            sum2_H += H * H;
            sum_mu += mu;
            sum_sigma += sigma;
        }

        // Deviazione standard sulla stima di H  a temperatura costante.
        err = sqrt(sum2_H / L - sum_H * sum_H / L / L);

        //stime progressive
        sum_H = sum_H / L;
        mus.push_back(sum_mu / L);
        sigmas.push_back(sum_sigma / L);

        err_mu = sqrt(mus.at(k) / L - mus.at(k) * mus.at(k) / L / L);
        err_sigma = sqrt(sigmas.at(k) / L - sigmas.at(k) * sigmas.at(k) / L / L);
        /*ave_mu+= mus.at(k) ;
        ave2_mu+= (mus.at(k))*(mus.at(k)) ;

        ave_sigma+= sigmas.at(k) ;
        ave2_sigma+= (sigmas.at(k))*(sigmas.at(k)) ;
        */
        aveH += sum_H;
        ave2H += sum_H * sum_H;
        error = k ? sqrt((ave2H - pow(ave2H, 2)) / k) : PREC;

        cout << "Ave H  = " << aveH << "Ave H ^2 = " << ave2H << " ; error 1  = " << err << endl;
        // Stime mu e sigma

        H = aveH / (k + 1);
        // Formato: #blocco, media blocco, errore blocco
        out << k + 1 << "\t" << H << "\t" << error << "\t" << sum_H << "\t" << err << endl;
        out_par << mus.at(k) << "\t" << sigmas.at(k) << endl;
        cout << endl << "BLOCK =  " << k << ", T = " << 1. / beta << ", passo = " << delta_ann << endl;
        cout << " mu = " << mu << ", sigma = " << sigma << endl;
        cout << " Metropolis Acceptance Rate: " << (double) accepted_met / attempted_met << endl;

        k++;
        cout << "==============================================================" << endl;

        if (iter == 1e4) { // max number of iteration
            cout << "max n of iteration reacher: unable to find result with precision " << PREC << endl;
            break;
        }
    }

    cout << endl << "mu = " << mu << endl << "sigma = " << sigma << endl;
    out.close();
    out_par.close();

    out.open("SA_Results.dat", ios::app);
    out << mu << "\t" << err_mu << "\t" << sigma << "\t" << err_sigma << "\t" << H << "\t" << err << endl;
    out.close();
//Part 3: show a picture of the estimation of H and its statistical uncertainty as a function of the number
// of blocks/MC steps for the set of parameters which minimize H
    mu = mus.at(mus.size() - 1);
    sigma = sigmas.at(sigmas.size() - 1);

    int M = 100000;
    L = 100;
    /*int nblk = int(M / L);
    int N = 70000;

    eq_step = 100;
    double x = 1, sum; // Punto del campionamento con metropolis

    // Equilibrazione
    for (int i = 0; i < eq_step; i++) {
        Estrazione(rnd, x);
    }
    vector<double> H_values, H_av_values;

    //upload data on vector
     */
    cout << "Part 3 : Saving values for H for fixed values of mu and sigma \n";
    /*for (int i = 0; i < M; i++) {
        cout << "values = " << i << "/" << M << "\r" << flush;
        H_values.push_back(Integral(rnd, N));
    }
    H_av_values = Accumulate(H_values, nblk, L);
    Sum_Prog(H_av_values, "Part_3_H_nblk.dat", nblk, L);
*/

    double Hfix = Integral(rnd,M, L, "Part_3_H_nblk.dat");
//PUNTO 4: preparazione dell'istogramma per mostrare la densità di probabilità campionata
    out.open("Part_4_psi_modulus_squared.dat");
    mu = mus.at(mus.size() - 1);
    sigma = sigmas.at(sigmas.size() - 1);
    int max = 50000;
    double x = 1;
    cout << "Sampling |psi(x)|^2\n ";
    for (int i = 0; i < max; i++) {
        cout << "i = " << i << "/" << max << "\r" << flush;
        x = Metropolis(rnd, x);
        out << Metropolis(rnd, x) << endl;
    }

    out.close();

    rnd.SaveSeed();
    return 0;
}


/**********************************************************************************************
IMPLEMENTATION FUNCTIONS
***********************************************************************************************/

double Psi_test(double x) {
    // Funzione che calcola la funzione d'onda di test
    double a1 = (x - mu) / (sqrt(2.) * sigma);
    double a2 = (x + mu) / (sqrt(2.) * sigma);
    double psi = exp(-a1 * a1) + exp(-a2 * a2);
    return psi;
};


double V(double x) {
    return pow(x, 4) - 5 * x * x / 2;
};

double derivative(double x) {
    double alpha = pow((x + mu) / sigma, 2);
    double beta = pow((x - mu) / sigma, 2);
    return -1. / sigma / sigma * (alpha * exp(-alpha / 2) + beta * exp(-beta / 2) - exp(-alpha / 2) - exp(-beta / 2));
};

double K(double x) {
    return derivative(x) * hbar * hbar / (2 * m * Psi_test(x));
};

double Eloc(double x) {
    return V(x) + K(x);
};

double Metropolis(Random &rnd, double x) {
    // Campionamento della densità di probabilità tramite Metropolis
    attempted_met++;
    double xnew = x + rnd.Rannyu(-delta, delta);
    double psi_new2 = Psi_test(xnew) * Psi_test(xnew);
    double psi2 = Psi_test(x) * Psi_test(x);
    double p = fmin(1., psi_new2 / psi2); // Probabilità con cui accettare la mossa del Metropolis
    double thrw = rnd.Rannyu(); // sempre minore di 1 
    if (thrw < p) {
        accepted_met++;
        return xnew;
    } else return x;
}

double Estrazione(Random &rnd, double &x) {
    x = Metropolis(rnd, x);
    return Eloc(x);
};


double Integral(Random &rnd, int N) {
    eq_step = 500;
    double x = 1; // Punto del campionamento con metropolis
    double inc; //incertezza: errore progressivo

    // Equilibrazione
    for (int i = 0; i < eq_step; i++) {
        Estrazione(rnd, x);
    }
    double sum = 0;
    for (int k = 0; k < N; k++) {
        sum += Estrazione(rnd, x);
    }
    return sum / N;
};

double Integral(Random &rnd, int M, int L, string filename) {
    ofstream H_out;
    H_out.open("Part_3_H_nblk.dat", ios::app);
    int N = int(M/L);

    eq_step = 500;
    double x = 1, sum; // Punto del campionamento con metropolis
    double inc; //incertezza: errore progressivo

    // Equilibrazione
    for(int i=0; i<eq_step; i++)
        Estrazione(rnd, x);

    double avg=0;
    for(int k = 0; k<N; k++) {

        sum = 0;
        for(int j = 0; j<L; j++)
            sum += Estrazione(rnd, x);

        sum = sum/L;
        avg += sum;

        double avg2;
        avg2 += sum*sum;
        inc = k ? sqrt((avg2/(k+1)-pow(avg/(k+1),2))/k) : 0;
        H_out << k+1 << "\t" << avg/(k+1) << "\t" << inc << endl;
    }
    H_out.close();
    return avg/N;
};

double update_param(Random &rnd, double beta) {
    attempted_ann++;

    muold = mu;
    sigmaold = sigma;
    double H = Integral(rnd, 70000);

    mu += rnd.Rannyu(-delta_ann, delta_ann);
    sigma += rnd.Rannyu(-delta_ann, delta_ann);
    double Hnew = Integral(rnd, 70000);

    double P = (Hnew - H > 0) ? exp(-beta * (Hnew - H)) : 1;
    double y = rnd.Rannyu();

    if (Hnew - H > 0) {

        if (y < P) {
            accepted_ann++;
            H = Hnew;
        } else {
            mu = muold;
            sigma = sigmaold;
        }
    }
    return H;
};


template<typename T>
double error(vector <T> AV, vector <T> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt((AV2[n] - pow(AV[n], 2)) / n);
    };
};


vector<double> Accumulate(vector<double> &v, int n_blk, int L) {

    vector<double> v_acc;
    cout << "Accumulation\n";
    for (int i = 0; i < n_blk; ++i) {
        cout << "block " << i << "/" << n_blk << "\r" << flush;
        double sum = 0.;
        for (int j = 0; j < L; ++j) {
            int k = j + i * L;
            sum += v.at(k);
        }
        sum /= L;
        v_acc.push_back(sum);
        cout << "-";
    }
    cout << "End Accumulation\n";
    return v_acc;
};

void Sum_Prog(vector<double> &v, string filename, int n_blk, int L) {
    cout << "Creazione Files \n";
    vector<double> sum_prog, sum_prog2;
    ofstream output(filename, ios::app);
    for (int i = 0; i < n_blk; i++) {
        cout << "block " << i << "/" << n_blk << "\r" << flush;
        double sum = 0., sum2 = 0.;


        for (int j = 0; j < i + 1; j++) {
            sum += v.at(j);
            sum2 += v[j] * v[j];
        };

        // update vector
        sum_prog.push_back(sum / (i + 1));
        sum_prog2.push_back(sum2 / (i + 1));
        // output
        output << i << "\t" << v.at(i) << "\t" << sum_prog.at(i) << "\t" << error(sum_prog, sum_prog2, i) << endl;
    }
    cout << "output file created \n";
    output.close();
};

