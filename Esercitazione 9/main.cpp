/**************************************/
/*Laboratorio di Simulazione Numerica */
/*Studentessa: Sambruna Alessia 	  */
/*n matricola: 942397                 */
/**************************************/
/********************************************************************************
Esercitazione 9 :The Traveling Salesman Problem: the salesman must visit one and only
 one time every city and must be back to the first city in the end of the path.
 The challenge of the problem is that the traveling salesman wants to minimize
 the total length of the trip. the problem is represented as a 1D problem where the vector
 represents the order in which the cities are visited
*********************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include "function.h"
#include "random.h"
#include<cmath>
#include<cstdlib>
#include<vector>
#include <algorithm>
#include <random>

using namespace std;

//main
int main(int argc, char **argv) {
    Random rnd;
    random_setup(rnd);
    double R = 1.;
    double l = 1.;


    // Program operating mode
    if (argc != 2) {
        cerr << "Digit:\n0-> to study circle case;\n1-> to study square case \n";
        cerr << argv[0] << " <mode>\n" << argv[1];
        return -1;
    }

    int mode = atoi(argv[1]);
    int n_cities = 34;
    int N_paths = 500;
    int n_iteration = 5000; // number of iteration of the algorithm
    int n_couples = 150;

    Population population(n_cities, rnd, R, mode, N_paths); // original population,

    population.Length_Selection_Operator();// cities are randomly shuffled
    for (int i = 0; i < population.Get_Population().size(); ++i) {
        cout << "path " << i << " ==> D = " << population.Get_Population().at(i).total_distance() << endl ;
    }

    if (mode == 0) {
        population.Get_Population().at(0).print_path("circle_path.dat") ;
    } else {
        population.Get_Population().at(0).print_path("square_path.dat") ;
    }


    fstream out_len;
    if (mode == 0) {
        out_len.open("Circle_data_length.dat", ios::app);
    } else {
        out_len.open("Square_data_length.dat", ios::app);
    }

    for (int i = 0; i < n_iteration; ++i) {
        cout << "=============================\nIterazione "<< i <<"/"<< n_iteration<<  endl << endl;
        population.New_Gen(rnd, n_couples);
        population.Statistics() ;
        cout << "\nD = " << population.Get_Best()<< "\n\n=============================\n";
        out_len << i << "\t" << population.Get_Best() << "\t" << population.Get_Mean() << "\t" << population.Get_Sigma() << endl ;
    }
    out_len.close();

    // Part 2 : print the best path
    population.Length_Selection_Operator();


    if (mode == 0) {
        population.Get_Population().at(0).print_path("Best_Path_circle.dat");
    } else {
        population.Get_Population().at(0).print_path("Best_Path_square.dat");
    }

    cout << "------------------------------------\n";
    cout << "Best Path-> D= " << population.Get_Best() <<endl ;
    population.Get_Population().at(0).print_path();

    //end
    rnd.SaveSeed();

    return 0;
}
/*
    Paths P (20) ;
    P.circle_path_generator(rnd, R) ;
    P.print_path() ;

    for(int i = 0; i < 10 ; i++) {
    P.pair_permutation(rnd) ;

    }

    cout<< endl<< endl;
    P.print_path() ;

    //vector<int> v {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};

    /*for (int i = 0; i < v.size(); ++i) {
        cout << "v  = " << v.at(i) << endl ;
    }


    vector<int> w {v.begin() + 1, v.end()};// vector excluding the starting point
    int n = v.size() -1 ;
    int pos = rnd.Int( 0, n/2);
    int start1 = rnd.Int(0, n );
    int start2;
    do{
        start2 = rnd.Int( 0, n );
    }while(start2 < fmod(start1 + pos, n) );

    vector<int> tmp;
    for(int i = 0; i < 5; i++){
        tmp.push_back( w.at( fmod(start1 + i, n) ));
        w.at( fmod(start1 + i, n) ) = w.at( fmod(start2 + i, n) );
        w.at( fmod(start2 + i, n) ) = tmp.at(i);
    }
    v.erase(v.begin() + 1, v.end());
    v.insert(v.end(), w.begin(), w.end());

    for (int i = 0; i < v.size(); ++i) {
        cout << "v  = " << v.at(i) << endl ;
    }

    cout << endl ;

*/