/*******************************************************************************
/* Laboratorio di Simulazione Numerica*/
/* Studentessa: Sambruna Alessia 	  */
/* n matricola: 942397                */
/**************************************/
/* function.h: File containing the function used in exercise 1.1. 
***************************************************************************************************************************************/
#ifndef __functions_h__
#define __functions_h__

#include <iostream>
#include <fstream>
#include <string>
#include<cmath>
#include<cstdlib>
#include<vector>
#include "functions.h"
#include "random.h"

using namespace std;

template<typename T>
double error(vector <T> AV, vector <T> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(((AV2[n] - pow(AV[n], 2)) / n));
    };
};
// functions for block mean

vector<double> Accumulate(vector<double> &v, int nblk, int M) {
    if (M != v.size()) cerr << "v size must be " << M << endl;
    int L = int(M / nblk);
    vector<double> ave;
    for (int i = 0; i < nblk; i++) {
        cout << "blk n =\t" << i << "/" << nblk << "\r" << flush ;
        double sum = 0;
// calculate the average on the measurements on the block
        for (int j = 0; j < L; j++) {
            int k = j + i * L; // taking into account that
            sum = sum + v.at(k);
        }
        ave.push_back(sum / L); // each block has a dimention of L
        //av2.push_back(pow(ave[i], 2));
    };
    return ave;
};

void Progressive_Mean(vector<double> ave, int nblk, int L, string filename) {
    ofstream out;
    out.open(filename, ios::app);
    if (nblk != ave.size()) cerr << "input vector size must be " << nblk << endl;
    vector<double> av2;
    vector<double> sum_prog;
    vector<double> su2_prog;
    vector<double> err_prog;
    for (int i = 0; i < nblk; i++) {
        av2.push_back(pow(ave[i], 2));
    };

    for (int i = 0; i < nblk; i++) {
        cout << "blk n =\t" << i << "/" << nblk << "\r" << flush ;
        double sum = 0;
        double sum2 = 0;
        for (int j = 0; j < i + 1; j++) {
            sum += ave[j];
            sum2 += av2[j];
        };

        sum_prog.push_back(sum / (i + 1));
        su2_prog.push_back(sum2 / (1 + i));
        err_prog.push_back(error(sum_prog, su2_prog, i));

        out << i * L << " " << sum_prog[i] << " " << err_prog[i] << endl;
    };

    out.close();
};

vector<double> Chi_Quad_Calculator(Random &rnd, int Mtot, bool print) {
    ofstream hist_out;
    hist_out.open("Chi_quad_distr.dat");
    vector<double> chi_quads;
    for (int j = 0; j < Mtot; ++j) {

        double chi_quad = 0;
        int M = 100;
        int n_tot = 10000;
// estrmi del jesimo intervallo da analizzare
        double a = 0;
        double b = 1;
        int ni = 0;

        vector<int> n_throws; // vector of the "measured" throws

// per ogni intervallo [(m-1)/M, m/M]
        for (int m = 1; m <= M; m++) {
            ni = 0;
            int n_exp = n_tot / M;
            a = double(m - 1) / M;
            b = double(m) / M;

            // per ciascun intervallo si fa un test di 10^4 tiri
            for (int i = 0; i < n_tot; i++) {
                double temp = rnd.Rannyu();
                if (temp < b and temp >= a) ni++;
            }

            n_throws.push_back(ni);
            chi_quad += pow(ni - n_exp, 2) / n_exp;
        }
        chi_quads.push_back(chi_quad);
        cout << "Chi = " << chi_quad << "\t=> " << j << "/" << Mtot << "\r" << flush ;
        if (print == 0) {
            hist_out << chi_quad << endl;
        }
    }

    hist_out.close();
    return chi_quads;
};


//function_random_generator_setup -> DA CANCELLARE INCLUSA IN COSTRUTTORE DI rANDOM


//Function for statistical uncertainty estimation:
/*Description of the punction: the function estimate the uncertainty between 

Inputs:
AV: vector containing the cumulative averages 
AV2: vector containig the cumulative averages^2
n: number of block 
*/




// =====================================================================
template<typename T>
double CalcolaMedia(const vector <T> &v) {
    T accumulo = 0;
    for (int k = 0; k < v.size(); k++) {
        accumulo += v[k];
    }

    return accumulo / v.size();

};

template<typename T>
double CalcolaMediana(vector <T> v) {

    sort(v.begin(), v.end());  // Use the STL sort

    T mediana; // creo variabile che voglio ritornare
    int n = v.size() / 2; // posizione a metà
    if (v.size() % 2 == 0)//caso pari
    {
        mediana = (v[n - 1] + v[n]) / 2;
    } else //caso dispari
    {
        mediana = v[n];
    }

    return mediana;
};

template<typename T>
double CalcolaVarianza(const vector <T> &v) {

    T sumquad = 0;// inizializzo la somma dei quadrati
    for (int i = 0; i < v.size(); i++) {
        sumquad += pow((v[i] - CalcolaMedia(v)), 2);
    }

    return sumquad / v.size();

};


#endif

