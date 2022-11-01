#ifndef __functions_h__
#define __functions_h__


#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"


using namespace std;

// This file contains all the function and classes used in the exercise: 
/* INDEX:
 1 - Class Position
 2 -
*/

/*
Position class:
*/

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

    // destructor
    ~Position() {};

    // methods
    // return cartesian coordinates
    double getX() const { return m_x; };

    double getY() const { return m_y; };

    double getZ() const { return m_z; };

    // return spherical coordinates
    double getR() const { return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2)); };

    double getPhi() const { return atan(m_y / m_x); };

    double getTheta() const { return acos(m_z / getR()); };

    double getRho() const { return sqrt(pow(m_x, 2) + pow(m_y, 2)); };     // ray of the chilindric coordinates

    // calcuate the distance between the point and another point Position
    double Distance(const Position &P) const {
        return sqrt(pow(P.getX() - m_x, 2) + pow(P.getY() - m_y, 2) + pow(P.getZ() - m_z, 2));
    }; // distanza da un altro punto

    double Distance2(const Position &P) const {
        return (pow(P.getX() - m_x, 2) + pow(P.getY() - m_y, 2) + pow(P.getZ() - m_z, 2));
    };

    void step_x (double a, bool vs) // True : positive step; False : negative step
    {
        if (vs)
            m_x = m_x + a;
        else
            m_x = m_x - a;
    };

    void step_y (double a, bool vs) // True : positive step; False : negative step
    {
        if (vs)
            m_y = m_y + a;
        else
            m_y = m_y - a;
    };
    void step_z (double a, bool vs) // True : positive step; False : negative step
    {
        if (vs)
            m_z = m_z + a;
        else
            m_z = m_z - a;
    }

    void discrete_step (double a = 1.)
    {

        double v = m_rnd->Int(0, 2) ;  // random number to define the verse "Positive/negative" of the step
        double d = m_rnd->Rannyu(0, 3) ; // random number to define the direction of the step

        if (a<0) cerr << "Error: argument step must be positive\n" ;

        if (v==1)
        {
            a = a ;
        } else if (v==0)
        {
            a = -1*a ;
        } else{
            cerr << "Error in the random generator : check !\n" ;
        }

        if (d <1)
            m_x = m_x + a ;
        else if(d < 2 and d >=1)
            m_y = m_y + a ;
        else if (d < 3 and d >= 2)
            m_z = m_z + a ;
        else cout << "Error in the random generator : check !\n" ;
    };

    void continuous_step(double l = 1.) {
        double theta = m_rnd->Rannyu(-M_PI, M_PI); // definizione dell'angolo coseno
        double phi = m_rnd->Rannyu(0, 2 * M_PI);
        m_x = m_x + l * cos(phi) * sin(theta);
        m_y = m_y + l * sin(phi) * sin(theta);
        m_z = m_z + l * cos(phi);
    };

    void Print () {
        cout << "("<< m_x << "; " << m_y << "; " << m_z << ")" << endl ;
    }

    // return the distance between the starting point and the position
    double get_RW_distance(int n_steps, bool mode, double a = 1.){
        Position P_0 (0., 0., 0. ) ; // creation of a variable of type position, containing the actual position of the RW: discrete case
        m_x = m_y = m_z = 0. ;
        if( mode ){// discrete RW
            for (int i = 0; i < n_steps; ++i) {// numero di passi del RW
                discrete_step(a);
            }
        }else{// continuous RW
            for (int i = 0; i < n_steps; ++i) {// numero di passi del RW
                continuous_step(a);
            }
        }
        return Distance2(P_0) ;
    };

private:

    double m_x, m_y, m_z;
    Random *m_rnd = new Random();

};


// return a vector of N elements, of the mean for each block
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
    }
    return average;
};


//Progressive block mean functions
vector<double> Sum_Prog(int n_blk, vector<double> v, int M) {
    if (v.size()!=n_blk) cerr << "Error! vector size must be " << n_blk << endl ;
    int L = M / n_blk;
    vector<double> sum_prog;
    for (int i = 0; i < n_blk; i++) {
        double sum = 0;
        for (int j = 0; j < i + 1; j++) {
            sum += v[j];
        };
        sum_prog.push_back(sqrt(sum / (i + 1)));
    }
    return (sum_prog);
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
        sum_prog.push_back(sqrt(sum / (i + 1)));
    }
    return (sum_prog);
};

template<typename T>
double calc_error(vector <T> AV, vector <T> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(((AV2[n] - pow(AV[n], 2)) / n));
    };
};

void Results_output(string filename, int n_steps, int steps_tot, bool mode, int M , int N, double a =1.){

    // creation of output file
    ofstream out (filename, ios::app) ;
    if(mode ==0 ){
        cout << "\n___________________\nDiscrete Random Walk\n___________________\n" ;
    } else {
        cout << "\n___________________\nContinuous Random Walk\n___________________\n" ;
    }

    //RW simulation -> calculation of the distance each 10 steps
    while (n_steps < steps_tot) {
        Position P (0., 0., 0. ) ;// creation of a variable of type position, containing the actual position of the RW: discrete case
        vector<double> distances;

        for (int i = 0; i < M; ++i) {//distanza al quadrato
            distances.push_back(P.get_RW_distance(n_steps,mode, a));
        }
        vector<double> average = Accumulate(N, distances);
        vector<double> ave_prog = Sum_Prog(N, average, M);
        vector<double> ave2_prog = Sum2_Prog(N, average);

        double error_disc = calc_error(ave_prog, ave2_prog, N-1) ;

        cout << "n steps = " << n_steps << " ; Average = " << ave_prog.back() << "; Error = " << error_disc << endl ;
        out << n_steps << " " << ave_prog.back() << " " << error_disc << endl ;
        n_steps +=10 ;
    }
    out.close() ;
}

void Print_RW (bool mode, int n_steps ) {

    if (mode) {// discrete output
        Position P (0., 0., 0.) ;
        fstream out ("Discrete_RW.dat", ios::app) ;
        out << P.getX() << " " << P.getY() << " " << P.getZ() << endl ;
        for (int i = 0; i < n_steps; i ++) {
            P.discrete_step();
            out << P.getX() << " " << P.getY() << " " << P.getZ() << endl ;
        }
        out.close();
    }else{// continuous
        Position P (0., 0., 0.) ;
        fstream out ("Continuous_RW.dat", ios::app) ;
        out << P.getX() << " " << P.getY() << " " << P.getZ() << endl ;
        for (int i = 0; i < n_steps; i ++) {
            P.continuous_step();
            out << P.getX() << " " << P.getY() << " " << P.getZ() << endl ;
        }
        out.close();

    }
} ;


#endif