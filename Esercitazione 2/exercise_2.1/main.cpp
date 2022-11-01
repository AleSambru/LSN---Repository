/***************************************/
/* Laboratorio di Simulazione Numerica */
/* Studentessa: Sambruna Alessia 	   */	
/* n matricola: 942397				   */
/***************************************/

/* DA migliorare :
 - controllare che grafici siano sensati 
 - Fare una funzione per metodo della media
 - correggere i valori dell'importance sampling non completamente corretti 
 - scremare il codice 
 - Aggiungere considerazioni su python
 - caratterizzare, titolare, curare la grafica dei grafici
*/

/*Esercitazione 2 : exercise 02.1: Compute the following 1D integral via Monte Carlo : I = \int_0^1 \frac{\pi}{2}\cos(\pi x/2) dx = 1$$
	1. sampling a uniform distribution in $[0,1]$
	2. using importance sampling (i.e. sampling a non-uniform probability in $[0,1]$)
	>> Python: 
	Show a picture of both your estimations of $I$</span> and their uncertainties with a large number of *throws* M
	(e.g. $M > 10^4) as a function of the number of blocks, N.
*********************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "function.h"
#include "IntegralMC.h"

using namespace std;

//*______________MAIN_________________*

int main (int argc, char *argv[])
{
// variables 
 	Random rnd;
 	random_setup(rnd) ; 
	
	BaseFunction *f = new Integranda() ;
    BaseFunction *p = new Line(-2, 2) ;
    IntegralMC myintegral ;
	
	vector <double> I_Mean_MC ; 
	vector <double> I_IS_MC ; 

	double xmin = 0 ; 
	double xmax = 1 ; 
	int n_points = 1000 ; 

    /*for (int i = 0; i < 10 ; i ++) {
        cout << "VALORE = " << rnd.Fx(-2, 2) << endl ;
    }*/


//parameters for block method
	int M = 10000; // number of throws: number of points
	int N = 100; // number of blocks 
	int L = M/N ; // number of element for each block

	for (int i = 0 ; i < M; i++)// for each block
	{ 	
		I_Mean_MC.push_back( myintegral.Integral_Mean(rnd, n_points, xmin , xmax , f));
		I_IS_MC.push_back( myintegral.Integral_Imp_Sampl(rnd, n_points, xmin , xmax , f, p));
		cout << "Valore integrale media = " << I_Mean_MC[i] << endl ; 
		cout << "Valore integrale importance sampling= " << I_IS_MC[i] << endl ; 

	}
	

// mean block

vector <double> ave_mean ;//vector containing the mean for each block - mean method
vector <double> ave_IS ;  //vector containing the mean for each block - Importance sampling method

vector <double> ave2_mean ;//vector containing the mean for each block - mean method
vector <double> ave2_IS ;

	for (int i = 0 ; i < N; i++)// for each block
	{
		double sum_mean = 0. ; 
		double sum_IS = 0. ; 
		
		for (int j = 0 ; j < L ; j ++) // for each element of the block
		{
			int k = j + i*L ; 
			sum_mean+= I_Mean_MC[k] ; 
			sum_IS+= I_IS_MC[k]; 			
		}
		
		ave_mean.push_back(sum_mean/L ); 
		ave_IS.push_back(sum_IS/L ); 
		
		ave2_mean.push_back(pow(ave_mean[i], 2));
		ave2_IS.push_back(pow(ave_IS[i], 2));	
	}
	
vector <double> ave_prog_mean ;//progressive average - mean method
vector <double> ave_prog_IS ;

vector <double> ave2_prog_mean ;//progressive average - mean method
vector <double> ave2_prog_IS ;

vector <double> err_prog_IS ;
vector <double> err_prog_mean ;

	for (int i = 0 ; i <N; i ++)
	{
		double sum_mean = 0 ; 
		double sum2_mean = 0 ; 
		
		double sum_IS = 0 ; 
		double sum2_IS = 0 ; 
		
		for ( int j = 0 ; j < i+1; j ++)
		{
			sum_mean += ave_mean[j];  
        	sum2_mean += ave2_mean[j]; 		
        	
        	sum_IS += ave_IS[j];  
        	sum2_IS += ave2_IS[j]; 	
        			
		};
		
		
		ave_prog_mean.push_back(sum_mean/(i+1)) ;
		ave_prog_IS.push_back(sum_IS/(i+1)) ;
		
		ave2_prog_mean.push_back(sum2_mean/(1+i)); 
		ave2_prog_IS.push_back(sum2_IS/(1+i)); 
		
		err_prog_IS.push_back( error(ave_prog_IS,ave2_prog_IS,i)) ; 
		err_prog_mean.push_back( error(ave_prog_mean,ave2_prog_mean,i)) ; 
	}; 

	
// end block method
	
//output files 
	ofstream out_mean ("Integral_MC_mean.dat") ;
	if (!out_mean) cerr << "Errore file output\n" ;

	ofstream out_IS ("Integral_MC_IS.dat") ;
	if (!out_IS) cerr << "Errore file output\n" ;	

	for (int i = 0 ; i <N; i ++)
	{
		out_mean << i*L << " "<< ave_prog_mean[i]<< " " << err_prog_mean[i] <<  endl ; 
		out_IS << i*L << " "<< ave_prog_IS[i]<< " " << err_prog_IS[i] << endl ; 
	}
	
	
	out_mean.close() ; 
	out_IS.close() ; 

	// end
	   rnd.SaveSeed();
	   return 0;
}














