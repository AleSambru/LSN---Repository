#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include <iostream>
#include <fstream>
#include <vector>
#include "function.h"
//#include "random.h"
#include<cmath>
using namespace std; 



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


//function error
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




#endif
