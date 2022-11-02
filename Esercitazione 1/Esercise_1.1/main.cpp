/*******************************************************************************/
/* Laboratorio di Simulazione Numerica*/
/* Studentessa: Sambruna Alessia 	  */
/* n matricola: 942397                */
/**************************************/
/* Esercitazione 1: Test the given Pseudo-Random Number generator by:

>> Part 1:	estimating <r> = ∫ r dr = 1/2. 
			Make a picture of the estimation of <r> and its uncertainty as a function of the number of blocks N (notebook)
>> Part 2:	estimating σ^2= ∫ (r-1/2)^2 dr = 1/12.
			Make a picture of the estimation of σ^2 and its uncertainty as a function of the number of blocks, N (notebook)
>> Part 3: 	Divide [0,1] into M identical sub-intervals and implement the χ^2 test. The number of expected events observed in each sub-
			interval after n throws, according to a uniform distribution, is np = n * 1/M= n/M. Fix M=10^2 and use for n the first 10^4
			pseudo-random numbers, then the successive 10^4 pseudo-random numbers, and so on ... 100 times. In this case the chi-square
			statistic is: χ^2 = Σ_{i=1}^M ( n_i - n/M )^2/(n/M)
			We should expect on average that $(n_i - n/M)^2 \simeq n/M$ and thus $\chi^2 \simeq 100
***************************************************************************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "functions.h"
#include<cmath>
#include<cstdlib>
#include<vector>

using namespace std;


int main(int argc, char *argv[]) {
    //Variables for random methods-> contained in random.h
    Random rnd;
//Exercices 1.1
    cout << "===================\nPart 1 \n===================\n ";
// variables 
    int M = 100000;     // Total number of throws
    int N = 100;        // Number of blocks
    int L = int(M / N); // Number of throws in each block, please use for M a multiple of N

//averages vectors: vectors of N elements, one for each block
    vector<double> ave;
    vector<double> av2;
// progressive sums vectors: vectors of N elements, one for each block
    vector<double> sum_prog;
    vector<double> su2_prog;

    vector<double> err_prog;

    vector<double> r;  //vector with all the measurements
    vector<double> x; // vector with integer values from 0 to N blocks: number of throws at each progressive block

//filling of the vector with the value of measurements
    for (int i = 0; i < M; i++)
        r.push_back(rnd.Rannyu());

    for (int i = 0; i < N; i++)
        x.push_back(i * L);

// calculation of the average on each block
    for (int i = 0; i < N; i++) {
        double sum = 0;
// calculate the average on the measurements on the block 
        for (int j = 0; j < L; j++) {
            int k = j + i * L; // taking into account that
            sum = sum + r[k];
        }
        ave.push_back(sum / L); // each block has a dimention of L
        av2.push_back(pow(ave[i], 2));
    };
// Cumulative square average
//for each block I have a measurement -> use these measurement to calculate the mean
    for (int i = 0; i < N; i++) {
        double sum = 0;
        double sum2 = 0;
        for (int j = 0; j < i + 1; j++) {
            sum += ave[j];
            sum2 += av2[j];
        };

        sum_prog.push_back(sum / (i + 1));
        su2_prog.push_back(sum2 / (1 + i));
        err_prog.push_back(error(sum_prog, su2_prog, i));
    };


// Printing results : Print a file number of throws ave[j]-0.5 err_prog[j]
    ofstream out_1("Exercise_1.1_values.dat");
    if (!out_1) cerr << "Errore file output\n";
    for (int i = 0; i < N; i++) {
        out_1 << i * L << " " << (sum_prog[i] - 0.5) << " " << err_prog[i] << endl;
    }
    out_1.close();
// end exercise 1.1
    cout << "===================\nPart 2 \n===================\n ";

    // exercise 1.2:
//averages vectors: vectors of N elements, one for each block
    vector<double> std;
    vector<double> std2;

// progressive sums vectors: vectors of N elements, one for each block
    vector<double> prog_std;
    vector<double> prog_std2;

// calculation 
    for (int i = 0; i < N; i++) {
        double sum = 0;
        // calculate the average on the measurements on the block
        for (int j = 0; j < L; j++) {
            int k = j + i * L;
            sum += pow(r[k] - 0.5, 2);
        }
        std.push_back(sum / L);
        std2.push_back(pow(std[i], 2));
    }

    for (int i = 0; i < N; i++) {
        double sum = 0;
        double sum2 = 0;
        for (int j = 0; j < i + 1; j++) {
            sum += std[j];
            sum2 += std2[j];
        };

        prog_std.push_back(sum / (i + 1));
        prog_std2.push_back(sum2 / (i + 1));
    };

    vector<double> err_prog_std;
    for (int i = 0; i < N; i++)
        err_prog_std.push_back(error(prog_std, prog_std2, i));


// Print a file number of throws ave[j]-1/12 err_prog[j]
    ofstream out_2("Exercise_1.2_values.dat");
    if (!out_2)
        cerr << "Errore file output\n";

    for (int i = 0; i < N; i++) {
        out_2 << i * L << " " << prog_std[i] - double(1 / 12) << " " << err_prog_std[i] << endl;
    }

    out_2.close();


//Esercizio 1.3 : Divide  [0,1]  into  𝑀  identical sub-intervals and implement the  𝜒2  test.

    cout << "Part 3 - Chi Quad values \n";
    int Ntot = 10000;
    int nblk = 100;
    N = Ntot / nblk;
    cout << "Extimation of Chi quad values \n";
    vector<double> chi_quads = Chi_Quad_Calculator(rnd, Ntot, 0);
    
   	Single_Chi_Quad_Calculator(rnd, 0); 
    cout << "\n##########################\n ## Average calculator ## \n\n";
    vector<double> ave_chi_quads = Accumulate(chi_quads, nblk, Ntot);

    cout << "\n##########################\n ## Prog ave calculator ## \n\n";
    Progressive_Mean(ave_chi_quads, nblk, L, "Ex3_chi_quad_values.dat");

	cout << "Part 3 - Single Chi Quad values \n";
	
	
	

// end
    rnd.SaveSeed();
    return 0;
}
