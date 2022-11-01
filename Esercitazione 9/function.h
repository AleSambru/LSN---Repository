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

int wd = 20;

void random_setup(Random &rnd) {

    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
};

// function
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

class City {
public:
    // costruttori
    City() {
        m_x = 0;
        m_y = 0;
    };

    City(double x, double y) {
        m_x = x;
        m_y = y;
        m_label = 0;
    };

    City(double x, double y, int label) {
        m_x = x;
        m_y = y;
        m_label = label;
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

    void setCity(City C) {
        m_x = C.getX();
        m_y = C.getY();
        m_label = C.getLabel();
    }

    //distance between the argument and the private member
    double Distance(const City &C) const {
        return sqrt(pow(C.getX() - m_x, 2) + pow(C.getY() - m_y, 2));
    }; // distanza da un altro punto

    // generators
    void circle_position_gen(Random &rnd, double R = 1., int label = 0)// posizioni random su una circonferenza
    {
        double theta = rnd.Rannyu(0, 2 * M_PI); // polar angle
        m_x = R * cos(theta);
        m_y = R * sin(theta);
        m_label = label;
    };

    void square_position_gen(Random &rnd, double l = 1., int label = 0)// random positions inside a square
    {
        m_x = rnd.Rannyu(-1 * l, l);
        m_y = rnd.Rannyu(-1 * l, l);
        m_label = label;
    };

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
    void print() { cout << m_x << setw(wd) << m_y << setw(wd) << endl; };

private:
    double m_x, m_y;
    int m_label;
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

    Paths(vector <City> p) {
        m_n = p.size();
        for (int i = 0; i < m_n; ++i) {
            m_data.push_back(p.at(i));
            m_path.push_back(i);
        }
    };

    // destructor
    ~Paths() { ; };

    // print methods
    void print_path() {
        for (int i = 0; i < m_n; i++) {
            cout << "Città " << i << setw(3) << " label = " << m_current_path[i].getLabel() << "=>( " << setw(3)
                 << m_current_path[i].getX() << setprecision(3) << setw(3);
            cout << " ; " << setw(3) << m_current_path[i].getY() << setprecision(3) << setw(3) << " )" << endl;
        }
    };

    void print_path(int j) {
        cout << "Vector" << j << endl << "____________________" << endl;
        cout << "Numero di città = " << m_n << endl;
        for (int i = 0; i < m_n; i++) {
            cout << "Città " << i << " ( " << setw(3) << m_data[i].getX() << setprecision(3) << setw(3) << " ; "
                 << setw(3);
            cout << m_data[i].getY() << setprecision(3) << setw(3) << " )" << endl;
        }
    };

    void print_path(string filename) {
        ofstream out;
        out.open(filename, ios::app);
        for (int i = 0; i < m_n; i++) {
            out << m_current_path[i].getLabel() << setw(wd) << m_current_path[i].getX() << setw(wd)
                << m_current_path[i].getY() << endl;
        }
        // printing the initial city
        out << m_current_path.at(0).getLabel() << setw(wd) << m_current_path.at(0).getX() << setw(wd)
            << m_current_path.at(0).getY() << endl;
        out.close();
    };

    // path generators
    void circle_path_generator(Random &rnd, double R) {
        City C; // mlabel posta a zero
        for (int i = 0; i < m_n; i++) {
            C.circle_position_gen(rnd, R, i);
            m_data.push_back(C); // definisco l'elemento iesimo, cioè quanto vale City
            m_current_path.push_back(C);
            m_path.push_back(i);
        }
    };

    void square_path_generator(Random &rnd, double l) {
        City C;
        for (int i = 0; i < m_n; i++) {
            C.square_position_gen(rnd, l, i);
            m_data.push_back(C); // definisco l'elemento iesimo, cioè quanto vale City
            m_current_path.push_back(C);
            m_path.push_back(i);
        }
    };

    //distance calculator
    double total_distance() {
        double sum = 0.;
        for (int i = 1; i < m_n; i++) {
            sum += m_current_path.at(i - 1).Distance(
                    m_current_path.at(i));// calcola la distanza tra il punto precedente e il successivo
        }

        // calculating the distance to come back to beginning of the path
        sum += m_current_path.at(m_n - 1).Distance(m_current_path.at(0));
        m_distance = sum; // saving distance
        return sum;
    };

    // Get methods
    double Get_total_distance() const { return m_distance; };

    vector <City> Get_Current_Path() const { return m_current_path; };

    vector <City> Get_Data() const { return m_data; };

    vector<int> Get_label_path() const { return m_path; };

    int Get_n_cities() const { return m_n; };

    // Set methods
    void Set_City_Order() {
        for (int i = 0; i < m_n; i++) {
            // trova la città iesima
            int index = m_path[i];
            m_current_path.at(i).setCity(m_data.at(index));

        }
    };

    void Set_City_Order(vector<int> ord_path) {
        for (int i = 0; i < m_n; i++) {
            // trova la città iesima
            int index = ord_path.at(i);
            m_current_path.at(i).setCity(m_data.at(index));
        }
    };

    void Check_path() {//move the first city at the begining
        auto it = find(m_path.begin(), m_path.end(), 0);
        vector<int> cut(it, m_path.end());
        m_path.erase(it, m_path.end());
        m_path.insert(m_path.begin(), cut.begin(), cut.end());
        Set_City_Order();
    };
    //====================
    //=Mutation operators=
    //====================

    void mutation(Random &rnd) {
        int mode = rnd.Int(0, 5); // chosing which operator to perform
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
                case 5:
                    //cout << "Performing pair permutation\n";
                    pair_permutation(rnd);
                    Check_path();
                    break;
        }
    };

    // swap two cities
    void swap_city_operator(Random &rnd) {
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
    void shift(Random &rnd) {
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

    void inversion_operator(Random &rnd) {
        int pos = rnd.Int(1, m_n / 2);
        vector<int> w{m_path.begin() + 1, m_path.begin() + pos};// vector excluding the starting point
        vector<int> r{m_path.begin() + pos, m_path.end()}; // copy of m_path from the starting point
        m_path.erase(m_path.begin() + 1, m_path.end()); // only the first element is kept
        m_path.insert(m_path.end(), r.begin(), r.end()); // the last part of the original vector is put at the beginning
        m_path.insert(m_path.end(), w.begin(), w.end()); // the first part is put ad the end

        Set_City_Order();
    };

    void successive_permutation(Random &rnd) {// permutazione di un certo indice con il successivo

        vector<int> w{m_path.begin() + 1, m_path.end()};// vector excluding the starting point
        int i = rnd.Int(0, m_n - 1);
        int j = (i + 1) % (m_n - 1);
        swap(w.at(i), w.at(j));
        m_path.erase(m_path.begin() + 1, m_path.end());
        m_path.insert(m_path.end(), w.begin(), w.end());

        Set_City_Order();
    };

    void pair_permutation(Random &rnd) {
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

    void pair_swap(Random &rnd) {
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

    //===================
    //===Crossing Over===
    //===================

    void Crossingover(Random &rnd, vector<int> partner) {
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


    void operator=(Paths i) {
        this->m_data = i.Get_Data();
        this->m_current_path = i.Get_Current_Path();
        this->m_path = i.Get_label_path();
    }

private:
    int m_n; // number of cities
    vector <City> m_data; // "Cities Database"
    vector <City> m_current_path;
    vector<int> m_path; // label vector representing the visiting order of the cities
    double m_distance;
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

//===================
//=Classe Popolazione
//===================
class Population {
public:

    // costruttore
    Population() = default;

    Population(int n_cities, Random &rnd, double R, int &mode, int N) {// starting population

        m_N = N;
        Paths first_path(n_cities);
        if (mode == 0) {
            first_path.circle_path_generator(rnd, R);
        } else {
            first_path.square_path_generator(rnd, R);
        }

        m_mutability = 0.15;
        m_p = 1.5;
        m_n = n_cities;

        m_generation.push_back(first_path);
        for (int i = 0; i < m_N - 1; i++) {
            first_path.mutation(rnd);
            m_generation.push_back(first_path);
        }
        if (m_N != m_generation.size()) cerr << "Error in the creation of population!" << endl;
        Length_Selection_Operator();
    };

    void operator=(Population i) {
        this->m_generation = i.Get_Population();
    }

    // print operator
    void Print() {
        //    cout << "Printing individuals of the population\n--------------------------------------\n\n";
        //  cout <<  m_generation.size()<<"Dimensione del vettore "  << endl ;
        for (int i = 0; i < m_generation.size(); i++) {
            cout << "=================\nIndividuals " << " = " << i << endl;
            cout << "=================" << endl;
            m_generation.at(i).print_path();
        }
    };

    // Get
    vector <Paths> Get_Population() {
        return m_generation;
    };

    double Get_Mean() { return m_mean; };

    double Get_Sigma() { return m_sigma; };

    double Get_Best() {
        Length_Selection_Operator();
        return m_best;
    };

    /*void Selection_Operator(Random &rnd) {
        Length_Selection_Operator(); // ordering the vector of paths
        double r = rnd.Rannyu(); // mutation parameter
        int j = int(m_N * pow(r, m_p)) + 1; // chosing the path to change
        if (rnd.Rannyu() < m_mutability) {
            m_generation.at(j).mutation(rnd);
        }
        m_generation.at(j).Check_path();
        Length_Selection_Operator();
    };*/

    // order the paths of the popolation from the shortest to the longest
    void Length_Selection_Operator() {
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

    void Check_population() {
        int delta = m_N - m_generation.size();
        for (int i = 0; i < m_N; i++) {
            m_generation.at(i).Check_path();
        }

    };

    void Print_Best(string filename) {
        Length_Selection_Operator();
        m_generation.at(0).print_path(filename);
    };

    void Statistics() {
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

    // generate a new individuals to substitute to the generation with one of the worst
    void New_Gen(Random &rnd, int n_couples) {
        for (int i = 0; i < n_couples; i++) {
            Length_Selection_Operator();
            int j = int(ceil(m_N * pow(rnd.Rannyu(),m_p)))-1;
            int k = int(ceil(m_N * pow(rnd.Rannyu(),m_p)))-1;

            //cout << "Pairing = " << i << "/" << n_couples << "\r" << flush;
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

    void Replace(Random &rnd, Paths paths) {
        int n = m_N - 1;
        int x = rnd.Int(int(m_N / 2), n);
        m_generation.at(x) = paths;

    };


private:
    vector <Paths> m_generation;
    double m_mutability;
    double m_p; // selection operator exponent
    int m_N; // number of individuals : number of different path
    int m_n; // number of cities in each path
    double m_mean; // statistics parameters
    double m_sigma;
    double m_best; // distance of the best path
};



/*
void Selection_Operator(Random &rnd, vector<Paths> &p){
    Length_Selection_Operator (p) ; // ordering the vector of paths
    double r = rnd.Rannyu();
    double p = 0.4 ;
    int M = p.size() ;
    int j  = int (M*r*p) +1 ;

    Paths new_path (p.at(j)) ;
    double m_mutabilità = 0.7 ;

    if (rnd.Rannyu() <m_mutabilità ){

        switch (rnd.Int(1,4)){
            case 1: figlio.Mutazione1();
                break;
            case 2: figlio.Mutazione2();
                break;
            case 3: figlio.Mutazione3();
                break;
            case 4: figlio.Shift();
                break;
        }
    }
   // figlio.Vincoli();
    //figlio.Fit(m_p);
    //Sostituisci(figlio);
}

Order();


}

*/


// Problem Operator

// permutation among 𝑚 contiguous cities (except for the first city) with other (different!) 𝑚 contiguous cities
// (𝑚<𝑁/2), e.g. [1,2,3,4,5]→[1,4,5,2,3]


// selected parents : a, b
// n_cut : indice a cui fare il cut del vettore
/*
void Crossingover (Random& rnd , double p , Paths &a, Paths &b, int n_cut ) {

    if (a.Get_Path().size() != b.Get_Path().size())  cerr << "Error in the n of city elements of paths\n" ;
    int size  = a.Get_Path().size() ;
    double r = rnd.Rannyu() ;
   // cout << r << endl  ;
    if (r < p ){// crossingover
        // copio tutto

        vector <City> a_tale, b_tale;// code dei percorsi
        vector <City> a_new, b_new ;

        for (int i = 0; i < n_cut ; ++i) { // up to n_cut the same path is kept
            a_new.push_back(a.Get_Path().at(i));
            b_new.push_back(b.Get_Path().at(i));
        }

        // definisco tutta la tale
        for (int i = 0; i < size ; ++i) {
            a_tale.push_back(b.Get_Path().at(i));
            b_tale.push_back(a.Get_Path().at(i));
        }

        // elimino i doppioni se trovo che nella tale sono stati presi elementi che erano già stati presi nella prima parte
        for(int i =  n_cut; i <size; i++){
            for(int j = 0; j < size; j++){
                //Verifico se b.get_city(j) è presente in missing_a
                //  Se non è presente, j++
                //  Se è presente: parents.at(0).set_city(i, b.get_city(j));
                auto e = find(a_new.begin(), a_new.end(), b.Get_Path(j));
                if ( e != missing_a.end() ){
                    parents.at(0).set_city(i, b.get_city(j));
                    missing_a.erase(e);
                    break;
                }
            }
        }

        Crossingover(Random& rnd , double p , Paths &a, Paths &b) {

            int n_cut = rnd.Int(1,m_lDNA-2);
            vector<int> patrimonio (m_DNA.begin(), m_DNA.begin()+index);
            int j=0;

            while (index < m_lDNA){
                if ( find(patrimonio.begin(), patrimonio.end(), parent2.GetGene(j)) == patrimonio.end() ) {
                    patrimonio.push_back (parent2.GetGene(j));
                    index++;
                }
                j++;
            }
            Individuo figlio (m_rnd, m_lDNA, patrimonio);
            return figlio;
        };


        for (int i = 0; i < n_cut; i++){
            for (int j = 0; j < a_tale.size(); j++){
                if(b_new[i].getLabel() !=b_tale[j].getLabel()){
                    b_new.push_back(b_tale.at(j)) ;
                }
            }
        }

        for (int i = n_cut; i < size  ; ++i) {

            a.Set_City(i, a_new.at(i)) ;
            b.Set_City(i, b_new.at(i)) ;

        }


    }
};

 */

/*


 */
// genetic mutation operator
// swap

/*
void shift(int n_shift) {

    int size = m_n;
    int app_size = size - 1;
    vector <City> w(app_size);
    int k = 0;
    for (int i = 1; i < size; i++) {
        w.at(k) = (m_path.at(i));
        k++;
    }
    vector <City> r(app_size);

    for (int j = 0; j < app_size; j++) {
        int t = (j + n_shift) % (app_size);
        // cout << "t = " << t << endl ;
        r.at(t) = w.at(j);
    }

    k = 0;
    for (int i = 1; i < size; i++) {
        m_path.at(i) = r.at(k);
        k++;
    }
};

void mutation_shift(Random rnd, double p, int n_shift) {
    double r = rnd.Rannyu();
    cout << r << endl;
    if (r < p) {

        shift(n_shift);
    }
};
*/


// general
template<typename T>
double error(vector <T> AV, vector <T> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt((AV2[n] - pow(AV[n], 2) / n));
    };
};



//At this point you can add a crossover operator (that fulfils the bonds and that you will call with probability 𝑝(𝑖)𝑐>50%
//given a selected mother and father, e.g.
//[1,2,3,4,5][1,5,4,3,2]
//
//    cut their paths at the same position:
//    [1,2,3|4,5][1,5,4|3,2]
//
//conserve the first part of the paths:
//[1,2,3|𝑋,𝑌][1,5,4|𝑋,𝑌]
//complete the paths with the missing cities adding them in the order in which they appear in the consort:
//[1,2,3|5,4][1,5,4|2,3]
//
//Here you are: you have two sons that you can add to the new generation ... and so on!


//inversion of the order in which they appear in the path of 𝑚
//cities (except for the first city and 𝑚≤𝑁), e.g. [1,2,3,4,5]→[1,4,3,2,5] for the inversion of the cities from 2 to 4.
//


/*
vector <City> shift(Paths P, int n_shift, int size) {
    int app_size = size - 1;
    vector <City> v = P.Get_Path();
    vector <City> w(size - 1);
    int k = 0;
    // filling w vector
    for (int i = 1; i < size; i++) {
        w.at(k) = (v.at(i));
        k++;
    }
    vector <City> r(app_size);

    for (int j = 0; j < app_size; j++) {
        int t = (j + n_shift) % (app_size);
        // cout << "t = " << t << endl ;
        r.at(t) = w.at(j);
    }

    k = 0;
    for (int i = 1; i < size; i++) {
        v.at(i) = r.at(k);
        k++;
    }

    return v;
};*/


// check function
// le città devono essere visitate una sola volta
// l'ultima città coincide con la prima
// se le città da visitare sono 34, il vettore path avrà 35 elementi


//===================
//Mutation operators=
//===================
// swap operator: take as input the vector to swap, the positions of the vector values to swap

/*template<typename T>
void inversion_operator(vector <T> &v, int i, int j) {
    if (i != 0 and j != 0) {
        swap(v.at(i), v.at(j));
    } else cerr << "Starting points can't be changed\n\n";
}
*/




// shift
/*
template<typename T>
void shift(int n_shift, vector <T> &v) {
    vector <T> w{v.begin() + 1, v.end()};// vector excluding the starting point
    vector <T> r(w); // copia di r
    for (int j = 0; j < w.size(); j++) {
        int t = (j + n_shift) % (w.size());
        r.at(t) = w.at(j);
    }
    int k = 0;
    for (int i = 1; i < v.size(); i++) {
        v.at(i) = r.at(k);
        k++;
    }
};

template<typename T>
void cont_perm(int pos, vector <T> &v) {
    vector <T> w{v.begin() + 1, v.begin() + pos};// vector excluding the starting point
    vector <T> r{v.begin() + pos, v.end()}; // copia di r
    v.erase(v.begin() + 1, v.end());
    v.insert(v.end(), r.begin(), r.end());
    v.insert(v.end(), w.begin(), w.end());
};
*/



/*
template<typename T> void pair_perm_operator(vector <T> &v, int i){
    if(i!=0) {
        swap(v.at(i), v.at(i+1));
    } else cerr << "Starting points can't be changed\n\n";
}
*/



// =====================================================================
template<typename T>
double CalcolaMedia(const vector <T> &v) {
    T accumulo = 0;
    for (int k = 0; k < v.size(); k++) {
        accumulo += v[k];
    }

    return accumulo / v.size();

};

template<typename T>
double CalcolaMediana(vector <T> v) {

    sort(v.begin(), v.end());  // Use the STL sort

    T mediana; // creo variabile che voglio ritornare
    int n = v.size() / 2; // posizione a metà
    if (v.size() % 2 == 0)//caso pari
    {
        mediana = (v[n - 1] + v[n]) / 2;
    } else //caso dispari
    {
        mediana = v[n];
    }

    return mediana;
};

template<typename T>
double CalcolaVarianza(const vector <T> &v) {

    T sumquad = 0;// inizializzo la somma dei quadrati
    for (int i = 0; i < v.size(); i++) {
        sumquad += pow((v[i] - CalcolaMedia(v)), 2);
    }

    return sumquad / v.size();

};


#endif

