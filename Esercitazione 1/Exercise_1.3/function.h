#ifndef __Function__
#define __Function__

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "function.h"
#include<cmath>
#include<cstdlib>
#include<vector>

using namespace std;

/*
2 - Position class:
class Position {

public:

    // constructors
    Position() {
        m_x = 0;
        m_y = 0;
        m_z = 0;
    };

    Position(double x, double y, double z) {
        m_x = x;
        m_y = y;
        m_z = z;

    };

    // distructor
    ~Position() {};

    // methods
    // return cartesians coordinates
    double getX() const { return m_x; };

    double getY() const { return m_y; };

    double getZ() const { return m_z; };

    void setX(double x) { m_x = x; };

    void setY(double y) { m_y = y; };

    void setZ(double z) { m_z = z; };

    // return spherical coordinates
    double getR() const { return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2)); };

    double getPhi() const { return atan(m_y / m_x); };

    double getTheta() const { return acos(m_z / getR()); };

    double getRho() const { sqrt(pow(m_x, 2) + pow(m_y, 2)); };     // ray of the cilindric coordinates

    // calcuate the distance between the point and another point Position
    double Distance(const Position &) const {
        return sqrt(pow(getX() - m_x, 2) + pow(getY() - m_y, 2) + pow(getZ() - m_z, 2));
    }; // distanza da un altro punto

private:

    double m_x, m_y, m_z;

};
*/

template<typename T>
double calc_error(vector <T> AV, vector <T> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(((AV2[n] - pow(AV[n], 2)) / n));
    };
};


/*Funzione check hit*/
// Considero intervallo 0, d
/*
This Function Perform the experiment of throwing a stick and checks wheter it has touched the lines. The function return the result of the experiment
*/

bool experiment(Random &rnd, double L, double d) {

    double x_A = rnd.Rannyu(0, d); // first point is thrown
    bool result = false;
    //Determino un angolo theta casuale
    double theta, x, y;
    bool x2_y2 = 0;
    while(x2_y2!=1){
        x = rnd.Rannyu();
        y = rnd.Rannyu();
        x2_y2 = x*x + y*y < 1;
    }
    theta = atan2(y, x);
    double cosx = cos(theta) ;
    double O = 0;
    double x_B = x_A + L * cosx;

    //initial checks: the line is touched if one of the points lays exactly on the line
    if (x_A == O or x_A == d or  x_B == 0 or x_B == d) result = true;
    if (x_A < 0 or x_A > d) cerr << "Error in the random generator!\n";
    //A is generated between 0 and d, therefore we check what happens to the B point
    if (x_B < O or x_B > d) {
        result = true;
    } else if (x_B > 0 and x_B < d) {
        result = false;
    } else {
        cerr << "Implementation Error;\n";
    }
    return result;
};


vector<double> Accumulate(int n_blk, vector<double> v) {
    int N = v.size();
    int L = N / n_blk;
    vector<double> average;
    for (int i = 0; i < n_blk; ++i) {
        double sum = 0.;
        for (int j = 0; j < L; ++j) {
            int k = j + i * L;
            sum += v.at(k);
        }
        sum = sum / L;
        average.push_back(sum);
        cout << "pi = " << sum << endl;
    }

    return average;
};

vector<double> Sum_Prog(int n_blk, vector<double> v) {
    int N = v.size();
    int L = N / n_blk;
    vector<double> sum_prog;
    for (int i = 0; i < n_blk; i++) {
        double sum = 0;
        for (int j = 0; j < i + 1; j++) {
            sum += v[j];
        };
        sum_prog.push_back(sum / (i + 1));
    }
    return sum_prog;
};

vector<double> Sum2_Prog(int n_blk, vector<double> v) {
    int N = v.size();
    int L = N / n_blk;
    vector<double> sum_prog;
    for (int i = 0; i < n_blk; i++) {
        double sum = 0;
        for (int j = 0; j < i + 1; j++) {
            sum += pow(v[j], 2);
        };
        sum_prog.push_back(sum / (i + 1));
    }
    return sum_prog;
};

void print(vector<double> &v, vector<double> &err, string file_name, int L) {
    ofstream out(file_name, ios::app);
    if (!out) cerr << "Errore file output\n";
    for (int i = 0; i < v.size(); i++) {
        //cout << "Value " << i << " = " << v.at(i) << " Error = " << err.at(i) << endl;
        out << i << " " << v.at(i) << " " << err[i] << endl;
    }
};

// Block Mean: DA CONTROLLARE 
/*

vector <double> error_bloch (int N, vector <double> ave)
{
	vector <double> error ; 
	vector <double> sum_prog;
	vector <double> su2_prog; 
	for (int i = 0 ; i <N; i ++)
	{
		double sum = 0 ; 
		double sum2 = 0 ; 
		for ( int j = 0 ; j < i+1; j ++)
		{
			sum += ave[j];  
        	sum2 += pow(ave[j], 2); 				
		};
	
	sum_prog.push_back(sum/(i+1)) ;
	su2_prog.push_back(sum2/(1+i)); 
	error.push_back( calc_error(sum_prog,su2_prog,i)) ; 
	}

	return error; 
}
*/

template<typename T>
double Average_calculator(const vector <T> &v) {
    T sum = 0;
    for (int k = 0; k < v.size(); k++) {
        sum += v[k];
    }

    return sum / v.size();

};


#endif


