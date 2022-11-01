/*******************************************************************************
Laboratorio di Simulazione Numerica
Studentessa: Sambruna Alessia 
n matricola: 942397
Esercitazione 1 :	

	Make 3 pictures with the histograms obtained filling them with  10^4  realizations of  𝑆𝑁=1𝑁∑𝑁𝑖=1𝑥𝑖  (for  𝑁=1,2,10,100 ), being  𝑥𝑖  	a
	random variable sampled throwing a standard dice (fig.1), an exponential dice (fig.2, use  𝜆=1 ) and a Lorentzian dice (fig.3, use  
	𝜇=0  and  Γ=1 )

*********************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include<cmath>
#include<cstdlib>
#include<vector>
using namespace std;

//function definition -> implementation at the end of the page 
	template <typename T> double error( vector<T> AV, vector<T> AV2, int n );
	void random_setup(Random &rnd) ; 
	void Data_Creation_Exp_Cauchy(int N, int M, Random &rnd) ; 
//*___________________________________*
//*______________main_________________*
//*___________________________________*
int main (int argc, char *argv[])
{
	//Variables for random methods-> contained in random.h
   Random rnd;
   random_setup(rnd) ; 

   int M = 10000; //M : number of extraction of a pseudo random number
	
	// N : number of repetition on each number
	Data_Creation_Exp_Cauchy(1,M,rnd) ; //N = 1
	Data_Creation_Exp_Cauchy(2,M,rnd) ; // N = 2
	Data_Creation_Exp_Cauchy(10,M,rnd) ; //N = 10
	Data_Creation_Exp_Cauchy(100,M,rnd) ; // N = 100

	// end
	   rnd.SaveSeed();
	   return 0;
}

//*___________________________________*
//*______function implementation______*
//*___________________________________*


//N : number of times repeating the measure
// M : number of points

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



template <typename T> double error( vector<T> AV, vector<T> AV2, int n )
{
	if (n==0)
	{
		return 0 ; 
	}
	else{
		return sqrt(((AV2[n] - pow(AV[n], 2))/n)); 
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













