#ifndef __functions_h__
#define __functions_h__

#include <iostream>
#include <fstream>
#include <string>
#include<cmath>
#include<cstdlib>
#include<vector>
#include "functions.h"

using namespace std; 

//Function for statistical uncertainty estimation:
/*Description of the punction: the function estimate the uncertainty between 

Inputs:
AV: vector containing the cumulative averages 
AV2: vector containig the cumulative averages^2
n: number of block 
*/

template <typename T> double error( vector<T> AV, vector<T> AV2, int n )
{
	if (n==0)
	{
		return 0 ; 
	}
	else{
		return sqrt((AV2[n] - pow(AV[n], 2)/n)); 
	}; 	
}; 


void random_setup(Random &rnd)//function_random_generator_setup
{

	int seed[4];
   	int p1, p2;
   	ifstream Primes("Primes");
   	if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
   	Primes.close();

   	ifstream input("seed.in");
   	string property;
   	if (input.is_open()){
      	while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
};


void Data_Creation_Exp_Cauchy(int N, int M, Random &rnd) 
{
	// creation of variables for names
	string n_times = to_string(N) ;
	string title_exp = "Esponenziale_"+n_times + ".dat" ; 
	string title_cauchy = "Cauchy_"+n_times+".dat"  ; 
	string title_unif = "Unif_"+n_times+".dat"  ; 
	
	// creation of files
	ofstream out_exp(title_exp) ; 
	if (!out_exp) cerr << "Errore file output\n" ;
	
	ofstream out_cauchy(title_cauchy) ; 
	if (!out_cauchy) cerr << "Errore file output\n" ;
	
	ofstream out_unif(title_unif) ; 
	if (!out_unif) cerr << "Errore file output\n" ;
	
	
			for (int i = 0 ; i < M ; i ++ ) 
		{	
			double sum_exp = 0 ; 
			double sum_cauchy = 0 ; 
			double sum_unif = 0 ; 
			for (int j = 0; j < N ; j++ )
			{
				sum_exp += rnd.Exp(1); 
				sum_cauchy += rnd.Cauchy(1, 0); 
				sum_unif += rnd.Rannyu() ; 
			}
		
			out_exp << (sum_exp)/N << endl;
			out_cauchy << (sum_cauchy/N) << endl;
			out_unif<< (sum_unif/N) << endl ; 
		}

	out_exp.close() ; 
	out_cauchy.close() ; 
	out_unif.close() ; 
};


#endif
