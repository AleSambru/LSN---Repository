#include "IntegralMC.h"
#include "random.h"

using namespace std ;

/*====================================================
 Mean Method: int f(x)dx
 =====================================================
 * - N points are sampled from a uniform distribution from a range [xmin, xmas]
 * - For each point f(x) is calculated
 * - f(x) are summed
 * - I = sum * (xmax - xmin )/N
 */
double IntegralMC::Integral_Mean(Random &rnd, int N, double xmin , double xmax , BaseFunction * f)
{
	double sum = 0 ; 
	for (int i = 0 ; i < N ; i ++)
	{
		//generation of point uniformly distribuited between x_min and max 
		double x = rnd.Rannyu(xmin,xmax);  	
		//evaluation of the function in that point 
		sum += f->Eval(x) ;
		// media
	}
	return sum*(xmax-xmin)/N ; 
};

/*====================================================
 Importance Sampling Method: int f(x)p(x)dx
 =====================================================
 * - N points are sampled from a distribution p(x) from a range [xmin, xmas], where p(x) is similar to the f(x)
 * - For each point f(x) is calculated
 * - f(x) are summed
 * - I = sum * (xmax - xmin )/N
 */

double IntegralMC::Integral_Imp_Sampl(Random &rnd,int N, double xmin , double xmax , BaseFunction * f,BaseFunction * p )
{
	double sum = 0 ; 
	for (int i = 0 ; i < N ; i ++)
	{
		//generazione di punti
		double x = rnd.Fx() ;
		//valuto il valore della funzione in quel punto
        double fi = (f->Eval(x))/(p->Eval(x));
        sum = sum + fi;
		//sum += f->Eval(x) ;
		// media
	}
	return sum*(xmax-xmin)/N ; 
};




 
