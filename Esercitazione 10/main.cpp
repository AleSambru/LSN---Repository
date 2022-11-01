/**************************************/
/*Laboratorio di Simulazione Numerica */
/*Studentessa: Sambruna Alessia 	  */
/*n matricola: 942397                 */
/**************************************/
/********************************************************************************
Parallelizzazione di MPI: 
- module load mpi/mpich-x86_64
- make 
- mpiexec -np 6 ./Esercizio_10.exe 6
// suddivido risultati in cartelle -> corrispondenti al caso in cui si sia lavorato con 2, 4, 6, 8, 10 ranks 
*********************************************************************************/
#include "mpi.h"
#include <iostream>
#include <fstream>
#include "function.h"
#include "random.h"
#include <string>
#include <cstdlib>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

    // Parameters
    Random rnd;
    int n_gen = 100 ;
    int n_it = 15000 ;
    int exchange_step = 5 ; // exchanging info frequency

    // Parameters for the use of MPI
    //vector <int> cores = {0, 1,  2,  3,  4,  5}; // core labels : defininig the cores exchanging info
    int size, rank;
    //int n_ranks = 6 ;

    // MPI Setup
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);// number of ranks whithin a specific communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);// which process am I ?
    MPI_Status stat1, stat2, stat3;
    //MPI_Request req ;
    //int itag=1; int itag2=2; int itag3 = 3;  // tags for the different messages that are sent

    // Code Setup - Setting the number of chores
    // note to the developer : the important thing is that the number of chores is even
    if (argc != 2) {
        cerr << "This program is inteded to do a parallelization using MPI working on 6 cores \n";
        cerr << argv[0] << " <mode>\n" << argv[1];
        return -1;
    }
    
    int n_ranks = atoi(argv[1]) ;
    if(n_ranks%2!=0)cerr<< "Error! n of ranks must be even!"<< endl ;    
    
    int n_msg = int(double (n_ranks)/2) ;  
    vector <int> cores (n_ranks) ;
    for(int i = 0 ; i < n_ranks; i++){
        cores[i] = i;
        } 
    int itag[n_msg] ; 
    for(int i = 0 ; i < n_msg; i++){
    itag[i] = i+1; 
    //cout << "itag["<< i << "] = " << itag[i] << endl ; 
    }
    //MPI_Status stat[n_msg] ; 
    
    //___________________
    // Start of The Code
    //___________________

    // creation of the population read from file
    Population population (rnd, "American_capitals.dat", n_gen);
    int n_cities = population.Get_n_cities();
    //  creation of output files
    fstream out;
    string output_filename = "American_capitals_Results_RankUsed_" + to_string(n_ranks) + ".dat" ;
    out.open(output_filename, ios::app);
    
    fstream out_best_path;
    string out_best_path_filename = "BestPath_" + to_string(n_ranks) + "_RanksUsed_" + to_string(rank) + ".dat" ;
    out_best_path.open(out_best_path_filename, ios::app);

    // Actual Parallelization
    /* Every exchange_step the cores exchange their best path, the best path is exchanged using the appo
     * variable "best_path" and then the new path randomly substitute one of the paths of the population*/
     
    for (int i = 0; i < n_it; ++i) {
        cout << "=============================\nIterazione " << i << "/" << n_it << endl << endl;
        population.New_Gen(rnd, n_gen); // n_gen mutations
        // exchange of the best paths between chores
        if (n_it % exchange_step == 0){
            cout << "n it = " << n_it << "\nUse of parallelization\n";
            vector<int> best_path{population.Get_Best_Path()}; //
            random_shuffle(cores.begin(), cores.end());
                //1 , 2 
                if (rank == cores[0]) {
                    //cout << "core " << rank << " sending a message to core 1\n"; 
                    MPI_Send(&best_path[0], n_cities, MPI_INTEGER, cores[1], itag[0], MPI_COMM_WORLD);// sto partendo da j
                } else if (rank==cores[1]){   
                    //cout << "core " << rank << " receiving the message from core 0 " << endl ;             
                    MPI_Recv(&best_path[0], n_cities, MPI_INTEGER, cores[0], itag[0], MPI_COMM_WORLD, &stat1);                
                }
                // 2, 3
                else if (rank == cores[2]) {
                    MPI_Send(&best_path[0], n_cities, MPI_INTEGER, cores[3], itag[1], MPI_COMM_WORLD);// sto partendo da j
                }else if(rank==cores[3]){
                    MPI_Recv(&best_path[0], n_cities, MPI_INTEGER, cores[2], itag[1], MPI_COMM_WORLD, &stat2);                 
                }                
                // 3, 4 
                else if (rank == cores[4]) {
                    MPI_Send(&best_path[0], n_cities, MPI_INTEGER, cores[5], itag[2], MPI_COMM_WORLD);// sto partendo da j
                }else if(rank ==cores[5]){
                    MPI_Recv(&best_path[0], n_cities, MPI_INTEGER, cores[4], itag[2], MPI_COMM_WORLD, &stat3);
                }
             
            population.Replace(rnd, best_path) ;
         }

        // output info
        population.Statistics() ;
        cout << "\nD = " << population.Get_Best()<< "\n\n=============================\n";
        out << i << "\t" << population.Get_Best() << "\t" << population.Get_Mean() << "\t" << population.Get_Sigma() << endl ;
       
      }
          
    population.Get_Population().at(0).print_path(out_best_path_filename);
    // end of the program
    out.close() ;
    out_best_path.close(); 
    
    MPI_Finalize() ; //MPI finalization
    return 0;
}





