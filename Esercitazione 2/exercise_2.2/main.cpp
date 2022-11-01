/*******************************************************************************
Laboratorio di Simulazione Numerica
Studentessa: Sambruna Alessia 
n matricola: 942397
Esercitazione 2: Esercizio 2.2: 3D Random Walks (RW) on a cubic lattice and in the continuum
	>> C++ Script : 
	_______________________________________________________________________________________________
	Repeat many times (M = 10^4 ) the simulation of a random walk in 3D always starting at the origin
	- on a cubic lattice with lattice constant 𝑎=1: at each discrete time the walker makes a forward or backward step of length equal to
	𝑎 in one of the 3 principal directions of the lattice: 𝑥, 𝑦 or 𝑧. 
	- in the continuum; at each discrete time the walker makes a step of length equal to 𝑎(=1) along a random direction obtained by
	sampling uniformly the solid angle: 𝜃∈[0,𝜋] and 𝜙∈[0,2𝜋]

	>> Jupyter Notebook "Inserire nome del jupiter collegato allo script "
	_____________________________________________________________________________
	
	- Show a picture of sqrt ⟨|𝑟⃗ 𝑁|2⟩𝑅𝑊 for both RWs, with their statistical uncertainties, as function of the step
	𝑖∈[0,102].
*********************************************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include "functions.h"

using namespace std ; 

/*********** MAIN **********************/ 

int main (){
	// variables setup
	int M = 10000 ; // number of repetition
	int N = 100 ; // numbers of blocks
	int L = int(M/N) ; // dimention of each step

    double a = 1. ;
	

	int steps_tot  = 1000 ;
    int n_steps = 10 ;
	// output file
    string filename_disc = "Discrete_RW_distance_nstep.dat" ;
    string filename_cont = "Continuous_RW_distance_nstep.dat" ;
    Results_output(filename_disc, n_steps, steps_tot, 0, M, N, a) ;
    Results_output(filename_cont, n_steps, steps_tot, 1, M, N, a) ;

    // Print example RW
    Print_RW(0, 100) ;
    Print_RW(1 , 100) ;

    /*
    ofstream out_disc (filename_disc, ios::app) ;
    ofstream out_cont (filename_cont, ios::app) ;

    while (n_steps < steps_tot) {
        Position P (0., 0., 0. ) ;// creation of a variable of type position, containing the actual position of the RW: discrete case
        Position Q (0., 0., 0. ) ;
        vector<double> discrete_distances;
        vector<double> continuous_distances;

        for (int i = 0; i < M; ++i) {//distanza al quadrato
            discrete_distances.push_back(P.get_RW_distance(n_steps, 1, a));
            continuous_distances.push_back(Q.get_RW_distance(n_steps, 0, a));
        }
        vector<double> average_disc = Accumulate(N, discrete_distances);
        vector<double> ave_prog_disc = Sum_Prog(N, average_disc, M);
        vector<double> ave2_prog_disc = Sum2_Prog(N, average_disc);

        double error_disc = calc_error(ave_prog_disc, ave2_prog_disc, N-1) ;

        vector<double> average_cont = Accumulate(N, continuous_distances);
        vector<double> ave_prog_cont = Sum_Prog(N, average_cont, M);
        vector<double> ave2_prog_cont = Sum2_Prog(N, average_cont);

        double error_cont = calc_error(ave_prog_cont, ave2_prog_cont, N-1) ;

        out_disc << n_steps << " " << ave_prog_disc.back() << " " << error_disc << endl ;
        cout << "n steps = " << n_steps << " ; Average = " << ave_prog_disc.back() << "; Error = " << error_disc << endl ;
        out_cont << n_steps << " " << ave_prog_cont.back() << " " << error_cont << endl ;
        cout << "n steps = " << n_steps << " ; Average = " << ave_prog_disc.back() << "; Error = " << error_cont << endl ;

        n_steps +=10 ;
    }
    out_cont.close() ;
    out_disc.close() ;


*/
   return 0;

}
	
	


	
