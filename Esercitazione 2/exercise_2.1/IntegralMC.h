#ifndef __INTEGRALMC_H__
#define __INTEGRALMC_H__

//#include "IntegralMC.h"
#include <fstream>
//#include "function.h"
#include "random.h"
#include <cmath>
#include <vector>

using namespace std ; 

/*
In this file are defined:
> BaseFunction class
> GaussFunction class
> Integranda class
> IntegralMC class 

*/

class BaseFunction
{
	public:
	virtual double Eval(double x) const = 0;// è una classe virtuale pura
};


//function that gets integrated in the exercice
class Integranda: public BaseFunction
{
	public:
	Integranda(){;};
	~Integranda();
	virtual double Eval(double x) const {return M_PI*0.5*(cos(M_PI*x*0.5) );};
};

class Line : public BaseFunction {

public:

    Line() {m_m = 0; m_q = 0;};
    Line(double m, double q) {m_m = m; m_q = q;};

    ~Line() {;};

    double Get_q() const {return m_m;};
    double Get_m() const {return m_q;};
    void Set_m(double m) {m_m = m;};
    void Set_q(double q) {m_q = q;};
    double Inverse(double x) const {return (1-sqrt(1-x));};

    virtual double Eval(double x) const {return m_m*x + m_q;};

private:
    double m_m, m_q;

};

class GaussFunction: public BaseFunction 
{
public:
	GaussFunction(double mean, double sigma)
	{
		m_mean = mean; 
		m_sigma = sigma ; 
	};
	~GaussFunction();
	virtual double Eval(double x) const 
	{
		return (exp(pow(x - m_mean, 2)/(2*m_sigma*m_sigma) ))/sqrt(m_sigma*m_sigma*M_PI*2);
	}; 
private: 
double m_sigma ; 
double m_mean ;
}; 




//______________________________________
//Classe Integral MC
//____________________________________
class IntegralMC{
public:
	IntegralMC(){;};// costruttore 
	~ IntegralMC() {;} ; 
	double PROVA (double x) {return x  ; } ; 
	double Integral_Mean(Random &rnd, int N, double xmin, double xmax, BaseFunction *f) ;
    double Integral_Imp_Sampl(Random &rnd,int N, double xmin , double xmax , BaseFunction * f,BaseFunction * p ) ;

private:
	//Random m_rand;
	double m_errore ;
};


#endif
