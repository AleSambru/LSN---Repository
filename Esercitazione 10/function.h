
#ifndef __function_h__
#define __function_h__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "function.h"
#include "random.h"
#include <iomanip>
#include <algorithm>    // std::find_if_not
#include <array>
#include <numeric>
#include <iterator>
#include <functional>

using namespace std;
// Definition of Variables




//------------------
//--- Class City ---
//------------------
class City {
public:
    // costruttori
    City() {
        m_x = 0;
        m_y = 0;
    };

   City(double x, double y, int label, string state, string city) {
        m_x = x;
        m_y = y;
        m_label = label;
        m_city = city ;
        m_state = state ;
    };

    // distruttore
    ~City() {};

    // metodi
    // get
    double getX() const { return m_x; };        // Cartesian Coordinates
    double getY() const { return m_y; };
    double getPhi() const { return atan(m_y / m_x); };
    double getRho() const { return sqrt(pow(m_x, 2) + pow(m_y, 2)); };// Cilindric coordinates
    int getLabel() const { return m_label; };

    // set
    void setLabel(int n) { m_label = n; }

    void setCity(City C);

    //distance between the argument and the private member
    double Distance(const City &C) const {
        return sqrt(pow(C.getX() - m_x, 2) + pow(C.getY() - m_y, 2));
    }; // distanza da un altro punto


    //===================
    //operators overload
    //==================
    // operator assign
    City &operator=(const City &C) {
        m_x = C.getX();
        m_y = C.getY();
        m_label = C.getLabel();
        return *this;
    };

    // check if equal
    bool operator==(City &a) {
        bool c;
        if (a.getX() == m_x and a.getY() == m_y and a.getLabel() == m_label) c = true;
        else c = false;
        return c;
    };

    // print
    void print() { cout << m_x << setw(10) << m_y << setw(10) << endl; };

private:
    double m_x, m_y;
    int m_label;
    string m_state, m_city;
};

//=============
// Paths Class
//=============

class Paths // ciascun elemento della path è costituito da un dato percorso, quindi da un vettore di City, cioè di posizioni
{
public:
    // constructors
    Paths() = default;

    Paths(int n) {
        m_n = n; // numbers of cities that make the path -> size of the vector m_path
    };

    Paths(vector <City> p) {// ok
        m_n = p.size();
        for (int i = 0; i < m_n; ++i) {
            m_data.push_back(p.at(i));
            m_current_path.push_back(p.at(i));
            m_path.push_back( i);
        }
    };

    // destructor
    ~Paths() { ; };

    // print methods
    void print_path() ;

    void print_path(int j) ;

    void print_path(string filename) ;
    //distance calculator
    double total_distance();

    // Get methods
    double Get_total_distance() const { return m_distance; };

    vector <City> Get_Current_Path() const { return m_current_path; };

    vector <City> Get_Data() const { return m_data; };

    vector<int> Get_label_path() const { return m_path; };

    int Get_n_cities() const { return m_n; };

    // Set methods
    void Set_City_Order();

    void Set_City_Order(vector<int> ord_path) ;

    void Check_path();

    // Mutation operators
    void mutation (Random &rnd);
    void inversion_operator(Random &rnd);
    void shift(Random &rnd);
    void swap_city_operator(Random &rnd);
    void pair_permutation(Random &rnd);
    void successive_permutation(Random &rnd);
    void pair_swap(Random &rnd) ; 
    // Crossing over
    void Crossingover(Random &rnd, vector<int> partner);

    void operator= (Paths i) {
        this -> m_data = i.Get_Data();
        this -> m_current_path = i.Get_Current_Path();
        this ->m_path = i.Get_label_path() ;
    }

private:
    int m_n; // number of cities
    vector <City> m_data; // "Cities Database"
    vector <City> m_current_path;
    vector<int> m_path; // label vector representing the visiting order of the cities
    double m_distance;
};

// function definition
void remove(vector<int> &v);
bool check_equal(int &value_1);
void Length_Selection_Operator(vector <Paths> &v) ;

//==================
// Population Class
//==================
class Population {
public:

    // costruttore
    Population() = default;

    Population(Random &rnd, vector <City> &path_0, int N, int n_cities);

    Population(Random &rnd, string inputfile, int N_gen);

    void operator=(Population i) {
        this->m_generation = i.Get_Population();
    }

    // print operator
    void Print();

    // Get
    vector <Paths> Get_Population() {return m_generation;};
    double Get_Mean() { return m_mean; };
    double Get_Sigma() { return m_sigma; };
    double Get_Best() { return m_best; };
    int Get_n_cities() {return m_n;} ;
    vector<int> Get_Best_Path();
    void Selection_Operator(Random &rnd) ;
    void Length_Selection_Operator() ;
    void Check_population();
    void Print_Best(string filename);
    void Statistics() ;
    void New_Gen(Random &rnd, int n_couples) ;
    void Replace(Random &rnd, Paths paths);
    void Replace(Random &rnd, vector <int> best_path ) ;
private:
    vector <Paths> m_generation;
    double m_mutability;
    double m_p; // selection operator exponent
    int m_N; // number of individuals : number of different path
    int m_n; // number of cities in each path
    double m_mean;
    double m_sigma;
    double m_best;
};



#endif
