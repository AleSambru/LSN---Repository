//
// function.h implementation file
//
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

// Functions
void remove(vector<int> &v) {
    auto end = v.end();
    for (auto it = v.begin(); it != end; ++it) {
        end = std::remove(it + 1, end, *it);
    }

    v.erase(end, v.end());
};

bool check_equal(int &value_1) {
    bool c;
    if (value_1 == 0) c = true;
    else c = false;
    return c;
};

//============
// City Class
//============

void City::setCity(City C) {
    m_x = C.getX();
    m_y = C.getY();
    m_label = C.getLabel();
};


// ===========
// Class Path
// ===========
// print methods
void Paths::print_path() {
    for (int i = 0; i < m_n; i++) {
        cout << "Città " << i << setw(3) << " label = " << m_current_path[i].getLabel() << "=>( " << setw(3)
             << m_current_path[i].getX() << setprecision(3) << setw(3);
        cout << " ; " << setw(3) << m_current_path[i].getY() << setprecision(3) << setw(3) << " )" << endl;
    }
};

void Paths::print_path(int j) {
    cout << "Vector" << j << endl << "____________________" << endl;
    cout << "Numero di città = " << m_n << endl;
    for (int i = 0; i < m_n; i++) {
        cout << "Città " << i << " ( " << setw(3) << m_data[i].getX() << setprecision(3) << setw(3) << " ; "
             << setw(3);
        cout << m_data[i].getY() << setprecision(3) << setw(3) << " )" << endl;
    }
};

void Paths::print_path(string filename) {
    ofstream out;
    out.open(filename, ios::app);
    for (int i = 0; i < m_n; i++) {
        out << m_current_path[i].getLabel() << setw(10) << m_current_path[i].getX() << setw(10)
            << m_current_path[i].getY() << endl;
    }
    // printing the initial city
    out << m_current_path.at(0).getLabel() << setw(10) << m_current_path.at(0).getX() << setw(10)
        << m_current_path.at(0).getY() << endl;
    out.close();
};

//distance calculator
double Paths::total_distance() {
    double sum = 0.;
    for (int i = 1; i < m_n; i++) {
        sum += m_current_path.at(i - 1).Distance(
                m_current_path.at(i));// calcola la distanza tra il punto precedente e il successivo
    }

    // calculating the distance to come back to begining of the path
    sum += m_current_path.at(m_n - 1).Distance(m_current_path.at(0));
    m_distance = sum; // saving distance
    return sum;
};

// Set methods
void Paths::Set_City_Order() {
    for (int i = 0; i < m_n; i++) {
        // trova la città iesima
        int index = m_path[i];
        m_current_path.at(i).setCity(m_data.at(index));

    }
};

void Paths::Set_City_Order(vector<int> ord_path) {
    for (int i = 0; i < m_n; i++) {
        // trova la città iesima
        int index = ord_path.at(i);
        m_current_path.at(i).setCity(m_data.at(index));
    }
};

void Paths::Check_path() {
    //sposto la prima città all'inizio
    auto it = find(m_path.begin(), m_path.end(), 0);
    vector<int> cut(it, m_path.end());
    m_path.erase(it, m_path.end());
    m_path.insert(m_path.begin(), cut.begin(), cut.end());
    //cout << "Controllo fatto ! \n" ;
};

// Mutation operators
void Paths::mutation(Random &rnd)  {
        int mode = rnd.Int(0, 4); // chosing which operator to perform
        switch (mode) {
                case 0:
                    //cout << "Performing inversion\n";
                    swap_city_operator(rnd);
                    Check_path();
                    break;
                case 1:
                    //cout << "Performing shift\n";
                    shift(rnd);
                    Check_path();
                    break;
                case 2:
                    //cout << "Performing continuous permutation\n";
                    inversion_operator(rnd);
                    Check_path();
                    break;
                case 3:
                    //cout << "Performing pair permutation\n";
                    pair_swap(rnd);
                    Check_path();
                    break;
        }
    };
// swap two cities
    void Paths::swap_city_operator(Random &rnd) {
        // choosing appropriate indexes
        int i = rnd.Int(1, m_n);
        int j = rnd.Int(1, m_n);
        while (i == j) {
            i = rnd.Int(1, m_n);
            j = rnd.Int(1, m_n);
        }
        // swapping city at i and j (m_path)
        int temp = m_path.at(i);
        m_path.at(i) = m_path.at(j);
        m_path.at(j) = temp;
        // upload m_currentpath with the new order set by m_path
        Set_City_Order(); //reset of the current path
    };

    //cities are shifted of n shift steps
    void Paths::shift(Random &rnd) {
        int n_shift = rnd.Int(1, m_n - 1);
        vector<int> w{m_path.begin() + 1, m_path.end()};// vector excluding the starting point
        vector<int> r(w); // copia di r
        for (int j = 0; j < w.size(); j++) {
            int t = (j + n_shift) % (w.size());
            r.at(t) = w.at(j);
        }
        int k = 0;
        for (int i = 1; i < m_path.size(); i++) {
            m_path.at(i) = r.at(k);
            k++;
        }
        Set_City_Order();
    };

    void Paths::inversion_operator(Random &rnd) {
        int pos = rnd.Int(1, m_n / 2);
        vector<int> w{m_path.begin() + 1, m_path.begin() + pos};// vector excluding the starting point
        vector<int> r{m_path.begin() + pos, m_path.end()}; // copy of m_path from the starting point
        m_path.erase(m_path.begin() + 1, m_path.end()); // only the first element is kept
        m_path.insert(m_path.end(), r.begin(), r.end()); // the last part of the original vector is put at the beginning
        m_path.insert(m_path.end(), w.begin(), w.end()); // the first part is put ad the end

        Set_City_Order();
    };

    void Paths::successive_permutation(Random &rnd) {// permutazione di un certo indice con il successivo

        vector<int> w{m_path.begin() + 1, m_path.end()};// vector excluding the starting point
        int i = rnd.Int(0, m_n - 1);
        int j = (i + 1) % (m_n - 1);
        swap(w.at(i), w.at(j));
        m_path.erase(m_path.begin() + 1, m_path.end());
        m_path.insert(m_path.end(), w.begin(), w.end());

        Set_City_Order();
    };

    void Paths::pair_permutation(Random &rnd) {
        vector<int> w{m_path.begin() + 1, m_path.end()};// vector excluding the starting point
        int n = m_n - 1;
        int pos = rnd.Int(0, n / 2);
        int start1 = rnd.Int(0, n);
        int start2;
        do {
            start2 = rnd.Int(0, n);
        } while (start2 < fmod(start1 + pos, n));

        vector<int> tmp;
        for (int i = 0; i < pos; i++) {
            tmp.push_back(w.at(fmod(start1 + i, n)));
            w.at(fmod(start1 + i, n)) = w.at(fmod(start2 + i, n));
            w.at(fmod(start2 + i, n)) = tmp.at(i);
        }
        m_path.erase(m_path.begin() + 1, m_path.end());
        m_path.insert(m_path.end(), w.begin(), w.end());

        Set_City_Order();
    };

    void Paths::pair_swap(Random &rnd) {
        int start = 0;
        int finish = 0;

        while (start == finish) {
            start = rnd.Int(1, m_n);
            finish = rnd.Int(1, m_n);
        }
        if (start > finish)
            swap(start, finish);

        if (finish - start % 2 == 0) { //pari
            while (finish != start) {
                swap(m_path.at(start), m_path.at(finish));
                start += 1;
                finish -= 1;
            }

        } else { //dispari
            int i = 0;
            while (i < (finish - start) / 2) {
                swap(m_path.at(start), m_path.at(finish));
                start += 1;
                finish -= 1;
                i++;
            }
        }
    };


//===Crossing Over===
void Paths::Crossingover(Random &rnd, vector<int> partner) {
    // m_path is one of the parents
    int n = partner.size() - 1;
    int ncut = rnd.Int(1, n);
    // partner and initial vector must have the same size
    if (partner.size() != m_n) cerr << "Error in the n of city elements of paths\n";
    double prob = 0.8; // crossing over probability
    if (rnd.Rannyu() < 1) {// doing the crossing over
        vector<int> son{m_path.begin(), m_path.begin() + ncut};// first part to keep
        son.insert(son.end(), partner.begin(), partner.end());
        m_path.insert(m_path.begin() + ncut, partner.begin(), partner.end());
        remove(m_path);
    }

    Set_City_Order(m_path);
    return;
};


//---------
//Functions
//---------
template<typename T>
void Swap(T &a, T &b) {
    T temp = a;
    a = b;
    b = temp;
};

//====================================================================================================
// Length_Selection_Operator: takes as input a vector of paths, and order the vector according to the
// total distance for each path
//====================================================================================================
void Length_Selection_Operator(vector <Paths> &v) {
    int n;
    int N = v.size();
    // check on the number of cities for each path of the vector
    for (int i = 1; i < v.size(); i++) {
        n = v.at(i - 1).Get_n_cities();
        double n_new = v.at(i).Get_n_cities();
        if (n != n_new) cerr << "All paths ov vector v must have the same number of City elements\n";
    }

    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double d_bef = v.at(i).total_distance();
            double d_aft = v.at(j).total_distance();
            if (d_aft < d_bef) {
                Swap(v.at(i), v.at(j)); // swap
            }
        }
    }
};

//==================================================================================================================
// Selection Operator
// - order the paths vector
// - write a selection operator which simply uses the order in the ordered population with 𝑀
//   individuals, e.g. select the individual 𝑗 with the algorithm: 𝑗=𝑖𝑛𝑡(𝑀×𝑟𝑝)+1 where 𝑟 is a uniform random number 𝑟∈[0,1)
//   and 𝑝 a convenient exponent
//
//==================================================================================================================

//====================
//=Classe Popolazione
//====================
Population::Population(Random &rnd, vector <City> &path_0, int N, int n_cities) {// starting population

    if(path_0.size()!= n_cities) cerr<< "The number of cities of the path must be " << n_cities << endl ;

    m_n = path_0.size();
    m_N = N;
    Paths first_path(path_0);

    m_mutability = 0.15;
    m_p = 2.5;

    //cout << "first path created \n";

    m_generation.push_back(first_path);
    for (int i = 1; i < m_N; i++) {
        first_path.mutation(rnd);
        m_generation.push_back(first_path);
    }
    if (m_N != m_generation.size()) cerr << "Error in the creation of population!" << endl;
};

Population::Population(Random &rnd, string filename, int N_gen) {// starting population
    // temporary variables
    ifstream indata;
    string state, city;
    double city_x, city_y ;
    int counter = 0;
    vector<City> starting_path;

    // open file
    indata.open(filename);
    if(!indata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    // read file
    while(!indata.eof())
    {   indata >> city >> state >> city_x >> city_y ;
        City C(city_x, city_y, counter, state, city);
        starting_path.push_back(C) ;
        //cout << counter << "\t" << city  << "\t\t" <<state <<  "\t\t\t" << city_x << "\t" << city_y<< endl ;
        counter++;
    }
    indata.close() ;
    // definition of parameters
    m_n = starting_path.size(); // number of cities read
    m_N = N_gen; //number of paths of the generation
    Paths first_path(starting_path);

    m_mutability = 0.15;
    m_p = 2.5;
    m_generation.push_back(first_path);//
    for (int i = 1; i < m_N; i++) {
        first_path.mutation(rnd);
        m_generation.push_back(first_path);
    }
    if (m_N != m_generation.size()) cerr << "Error in the creation of population!" << endl;
};

// print operator
void Population::Print() {
    //    cout << "Printing individuals of the population\n--------------------------------------\n\n";
    //  cout <<  m_generation.size()<<"Dimensione del vettore "  << endl ;
    for (int i = 0; i < m_generation.size(); i++) {
        cout << "=================\nIndividuals " << " = " << i << endl;
        cout << "=================" << endl;
        m_generation.at(i).print_path();
    }
};


void Population::Selection_Operator(Random &rnd) {
    Length_Selection_Operator(); // ordering the vector of paths
    double r = rnd.Rannyu(); // mutation parameter
    int j = int(m_N * pow(r, m_p)) + 1; // chosing the path to change
    if (rnd.Rannyu() < m_mutability) {
        m_generation.at(j).mutation(rnd);
    }
    m_generation.at(j).Check_path();
    Length_Selection_Operator();
};

void Population::Length_Selection_Operator() {
    int n;
    int N = m_generation.size();
    // check on the number of cities for each path of the vector
    for (int i = 1; i < m_generation.size(); i++) {
        n = m_generation.at(i - 1).Get_n_cities();
        double n_new = m_generation.at(i).Get_n_cities();
        if (n != n_new) cerr << "All paths ov vector v must have the same number of City elements\n";
    }

    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double d_bef = m_generation.at(i).total_distance();
            double d_aft = m_generation.at(j).total_distance();
            if (d_aft < d_bef) {
                Swap(m_generation.at(i), m_generation.at(j)); // swap
            }
        }
    }
};

void Population::Check_population() {
    int delta = m_N - m_generation.size();
    for (int i = 0; i < m_N; i++) {
        m_generation.at(i).Check_path();
    }

};

void Population::Print_Best(string filename) {
    Length_Selection_Operator();
    m_generation.at(0).print_path(filename);
};

vector<int> Population::Get_Best_Path() {
    Length_Selection_Operator();
    return m_generation.at(0).Get_label_path();
};

void Population::Statistics() {
    Length_Selection_Operator();
    double half = 0;
    double sum = 0;
    double sum2 = 0;

    for (int i = 0; i < m_N; i++) {
        sum += m_generation.at(i).Get_total_distance();
        sum2 += pow(m_generation.at(i).Get_total_distance(), 2);

        if (i == int(m_N / 2))
            half = sum;
    }

    m_mean = half / double(int(m_N / 2));
    m_best = m_generation.at(0).Get_total_distance();
    m_sigma = sqrt(sum2 / double(m_N) - pow(sum / double(m_N), 2));
};

void Population::New_Gen(Random &rnd, int n_couples) {
        for (int i = 0; i < n_couples; i++) {
            Length_Selection_Operator();
            int j = int(ceil(m_N * pow(rnd.Rannyu(),m_p)))-1;
            int k = int(ceil(m_N * pow(rnd.Rannyu(),m_p)))-1;

            cout << "Pairing = " << i << "/" << n_couples << "\r" << flush;
            //int k = (int(m_N * pow(rnd.Rannyu(), m_p)) + 1) % m_N;
            //int j = (int(m_N * pow(rnd.Rannyu(), m_p)) + 1) % m_N;

            vector<int> partner = m_generation.at(k).Get_label_path();
            m_generation.at(j).Crossingover(rnd, partner);
            Length_Selection_Operator();

            j = (int(m_N * pow(rnd.Rannyu(), m_p)) + 1) % m_N;
            double r = rnd.Rannyu();
            if (r < m_mutability) {
                m_generation.at(j).mutation(rnd);
            }
            Replace(rnd, m_generation.at(j));
            Length_Selection_Operator();
        }
        Length_Selection_Operator();
        return;
    };

void Population::Replace(Random &rnd, Paths paths) {
    int n = m_N - 1;
    int x = rnd.Int(int(m_N / 2), n);
    m_generation.at(x) = paths;

};

void Population::Replace(Random &rnd, vector <int> best_path ){
    int n = m_N - 1;
    int x = rnd.Int(int(m_N / 2), n);
    m_generation.at(x).Set_City_Order(best_path) ; // I take the x element and set the new order ot this path
};

