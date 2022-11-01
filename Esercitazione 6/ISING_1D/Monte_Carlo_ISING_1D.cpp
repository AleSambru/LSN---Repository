/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"
//#include <format>
using namespace std;

//________________________________________________________MAIN_____________________________________________________
int main(int argc, char **argv) {
    //Inizialization
    Input();
    
    // Equilibration
    /*if (argc != 2) {
        cerr << "Digitare:\n0-> per valutare l'equilibrazione del sistema;\n1-> eseguire la simulazione;\n"
             << argv[0] << " <mode>\n" << argv[1];
        return -1;
    }
    int mode = atoi(argv[1]); //numero indicante la modalità di funzionamento
*/
    //  if (mode==0){


    //}
        metro = 0;// definition of method: 1:metropolis; 0:gibbs
        if (metro == 0) cout << "Gibbs Sampling Method\n";
        else if (metro == 1) cout << "Metropolis Method\n";
        else cerr << "Error in the definition of metro parameter!\n ";
        cout << "______________________________________________________________\n";
        cout << "Equilibration\n________________________\n\n";

        // Equilibration
        // Rappresentiamo l'andamento di diverse quantità in funzione del numero di step di equilibrazione per valutare il numero di step necessari nell'algoritmo
        
        // senza blocchi
        
        int n_eq_step = 10000; 
        
        vector <double> T  {1.,2., 1.5 } ; 
        vector <double> hs {2., 0.} ;
        
        // a blocchi 
        
            for (int i = 0; i < T.size(); i++) {// 3 diverse t iniziali
            	for (int j = 0; j < hs.size(); j++){ // 2 different h condition
            		
            		temp = T.at(i) ; 
            		beta = 1./temp ; 
            		h = hs.at(j) ;  		
         
            		cout << "=================\nInitial Condition\n" ; 
            		cout << "T = " << temp << "; h = " << h << endl ;            		
            		cout << "=================\n" << endl;
					
					ofstream ene_eq; 
					string file ; 
					if (metro == 0) file = "output.gibbs.eq.ene.h=" + to_string(h) + "t=" + to_string(temp) + ".dat" ;
        			else if (metro == 1) file = "output.metro.eq.ene.h=" + to_string(h) + "t=" + to_string(temp) + "";
            		
            		ene_eq.open(file, ios::app);
            		for (int iblk = 1; iblk <= nblk; iblk++) //Simulation
            		{
				        cout << "n block = " << iblk << "/" << nblk << "\r" << flush;
				        Reset(iblk);   //Reset block averages
				        for (int istep = 1; istep <= nstep; ++istep) {
				            //cout << "n step " << istep << "/" << nstep << "\r" << flush ;
				            Move(metro);
				            Measure();
				            Accumulate(); //Update block averages
				        }
				        cout << "\n\n";
				        Averages(iblk, temp);  //Print results for current block
				        ene_eq << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double) iblk << setw(wd) << err_u << endl;
            		}
            		         	 	
           	 	ene_eq.close() ; 
        	}
        }
        // avendo definito una nuova temperatura


        
        
        /*for (int iblk = 1; iblk <= nblk; iblk++) //Simulation
            {
                cout << "n block = " << iblk << "/" << nblk << "\r" << flush;
                Reset(iblk);   //Reset block averages
                for (int istep = 1; istep <= nstep; ++istep) {
                    //cout << "n step " << istep << "/" << nstep << "\r" << flush ;
                    Move(metro);
                    Measure();
                    Accumulate(); //Update block averages
                }
                cout << "\n\n";
                Averages(iblk, temp);  //Print results for current block
            }
            //ConfFinal(); //Write final configuration
        }*/
    


//CREAZIONE EFFETTIVA DEI GRAFICI DA INSERIRE POI IN UNA FUNZIONE/ DA RISCRIVERE IN MODO PIÙ CARINO

    // creation and opening files

/*
    ofstream Ene_t, Heat_t, Mag_t, Chi_t; // plot delle quantità in funzione della temperatura
    Ene_t.open("output.ene(t)", ios::app);
    Heat_t.open("output.heat(t)", ios::app);
    Mag_t.open("output.mag(t)", ios::app);
    Chi_t.open("output.chi(t)", ios::app);

    // filling vector
    for (int i = 0; i <= n_temp; i++) {
        cout << "----------- t = " << T[i] << "--------------" << endl;

        // avendo definito una nuova temperatura
        for (int iblk = 1; iblk <= nblk; iblk++) //Simulation
        {

            Reset(iblk);
            for (int istep = 1; istep <= nstep; istep++) {
                Move(metro);
                Measure();
                Accumulate();
            }
            Averages(iblk);
        }
        Ene_t << setw(wd) << T[i] << setw(wd) << glob_av[iu] / (double) nblk << setw(wd) << err_u << endl;
    }

    Ene_t.close();
*/

// FINE DELLE COSE DA RISCRIVERE IN MODO PIÙ CARINO

    return 0;
}

//__________________________________________________end of main ________________________________________________________
//__________________________________________________Functions___________________________________________________________
void Input(void) {
    ifstream ReadInput;

    cout << "Classic 1D Ising model\n________________________" << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed, p1, p2);
    input.close();

//Read input informations
    ReadInput.open("input.dat");

    ReadInput >> temp;
    beta = 1.0 / temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;

    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl << endl;

    ReadInput >> metro; // if=1 Metropolis else Gibbs

    ReadInput >> nblk;

    ReadInput >> nstep;

    ReadInput >> t_in; // lettura telle temperature iniziali e finali del range
    // entro cui studiare l'andoamento delle quantità

    ReadInput >> t_fin;

    ReadInput >> n_temp;

    delta_t = double((t_fin - t_in)) / n_temp;


    // filling tempertature vector
    for (int i = 0; i <= n_temp; i++) {
        T.push_back(t_in + delta_t * i);
    }


    if (metro == 1)
        cout << "The program perform Metropolis moves" << endl;// metro: indicates the modes in which we use the code
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();

//Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility

    n_props = 4; //Number of observables

//initial configuration: random initial configuraiton
    for (int i = 0; i < nspin; ++i) {
        if (rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
    }

//Evaluate energy etc. of the initial configuration
    Measure();

//Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu] / (double) nspin << endl;
}

// Makes a move
// metro = 0 -> Gibbs Sampling
// metro = 1 -> Metropolis
void Move(int metro) {
    int o;
    double p, energy_old, energy_new, sm, delta_E;
    double p_up, p_down;

    for (int i = 0; i < nspin; ++i) {
        //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
        o = (int) (rnd.Rannyu() * nspin);// prendo uno spin a caso

        if (metro == 1) //Metropolis
        {
            //Metodo di Metropolis:
            //Il metodo consiste nel partire da una configurazione casuale:
            // si definisce casualmente un punto del reticolo e si flippa lo spin in
            // in quella posizione: se si ha una diminuzione di energia teniamo
            // DA RISCRIVERE MEGLIO LA SPIEGAZIONE
            sm = s[o];
            energy_old = Boltzmann(sm, o);
            sm = -1 * sm;
            energy_new = Boltzmann(sm, o);
            delta_E = energy_new - energy_old;
            p = fmin(1, exp(-beta * delta_E));// probabilità di accettare la mossa
            if (rnd.Rannyu() < p) {
                s[o] = -s[o];
                accepted++;
            }
            attempted++;
        } else //Gibbs sampling
        {
            /*energy_up = Boltzmann(1, o);
            energy_down = Boltzmann(-1, o);
            //the move is always accepted
            p = 1. / (1. + exp(beta * (energy_up - energy_down)));
            if (rnd.Rannyu() < p) {
                s[o] = +1;
            } else {
                s[o] = -1;
            }*/
            p_up = exp(-beta * Boltzmann(1,o));
            p_down = exp(-beta * Boltzmann(-1,o));
            p = p_up / (p_up + p_down);

            if(rnd.Rannyu() < p)
                s[o] = +1;
            else
                s[o] = -1;

        }
    }

} // check


double Boltzmann(int sm, int ip) {
    double ene = -J * sm * (s[Pbc(ip - 1)] + s[Pbc(ip + 1)]) - h * sm;
    return ene;
}// check

// misura di tutte le quantità
void Measure() {
    //int bin;
    double u = 0.0, m = 0.0, u2 = 0.0, m2 = 0.0;

//cycle over spins
    for (int i = 0; i < nspin; ++i) {

        u += -J * s[i] * s[Pbc(i + 1)] - 0.5 * h * (s[i] + s[Pbc(i + 1)]);
        //u2 += pow((-J * s[i] * s[Pbc(i + 1)] - 0.5 * h * (s[i] + s[Pbc(i + 1)])), 2);
        m += s[i];
        //m2 += pow(s[i], 2);
    }

    // ###### --> check the formula for calculation of this quantities #####
    walker[iu] = u; // internal energy
    walker[ic] = u * u; //pow(beta, 2) * ((u2) - pow((u), 2)); //specific heat
    walker[im] = m; // magnetization
    walker[ix] = m * m; //beta * (m2 - pow(m, 2)); //suscieptivity

    //##########################################àà

} // check


void Reset(int iblk) //Reset block averages // check
{

    if (iblk == 1) {
        for (int i = 0; i < n_props; ++i) {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for (int i = 0; i < n_props; ++i) {
        blk_av[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}


void Accumulate(void) //Update block averages // check
{

    for (int i = 0; i < n_props; ++i) {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
   /* ofstream Ene, Heat, Mag, Chi;

    // openinig files
    Ene.open("output.ene.0", ios::app);
    Heat.open("output.heat.0", ios::app);
    Mag.open("output.mag.0", ios::app);
    Chi.open("output.chi.0", ios::app);
*/
    // calculating values
    // internal energy
    stima_u = (blk_av[iu] / blk_norm) / (double) nspin;
    glob_av[iu] += stima_u;
    glob_av2[iu] += stima_u * stima_u;
    err_u = Error(glob_av[iu], glob_av2[iu], iblk);

    //specific heat
    stima_c =
            pow(beta, 2) * (blk_av[ic] / blk_norm - pow(blk_av[iu] / blk_norm, 2)) / (double) nspin; // Calore specifico
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c * stima_c;
    err_c = Error(glob_av[ic], glob_av2[ic], iblk);

    // magnetization
    stima_m = (blk_av[im] / blk_norm) / (double) nspin;
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m * stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);

    // suscieptivity
    stima_x = beta * blk_av[ix] / blk_norm / (double) nspin;
    glob_av[ix] += stima_g;
    glob_av2[ix] += stima_g * stima_g;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);

    // printing results on files
    /*Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double) iblk << setw(wd) << err_u
        << endl;
    Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double) iblk << setw(wd) << err_c
         << endl;
    Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im] / (double) iblk << setw(wd) << err_m
        << endl;
    Chi << setw(wd) << iblk << setw(wd) << stima_g << setw(wd) << glob_av[ix] / (double) iblk << setw(wd) << err_g
        << endl;

    // closing files
    Ene.close();
    Heat.close();
    Mag.close();
    Chi.close();
    */
}


void Averages(int iblk, double temperature) //Print results for current block
{
    // assignining file names
    string str_temp = to_string(temperature);
    string str_metro = metro ? "metro" : "gibbs";
    //string str_h = format("{:.2f}", h);
    string str_h = to_string(h) ;
    ofstream Ene, Heat, Mag, Chi;

    // openinig files
    Ene.open("output.ene." + str_metro + "." + str_h + "(t=" + str_temp + ")", ios::app);
    Heat.open("output.heat." + str_metro + "." + str_h + "(t=" + str_temp + ")", ios::app);
    Mag.open("output.mag." + str_metro + "." + str_h + "(t=" + str_temp + ")", ios::app);
    Chi.open("output.chi." + str_metro + "." + str_h + "(t=" + str_temp + ")", ios::app);
    // calculating values
    // internal energy
    stima_u = (blk_av[iu] / blk_norm) / (double) nspin;
    glob_av[iu] += stima_u;
    glob_av2[iu] += stima_u * stima_u;
    err_u = Error(glob_av[iu], glob_av2[iu], iblk);

    //specific heat
    stima_c =
            pow(beta, 2) * (blk_av[ic] / blk_norm - pow(blk_av[iu] / blk_norm, 2)) / (double) nspin; // Calore specifico
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c * stima_c;
    err_c = Error(glob_av[ic], glob_av2[ic], iblk);

    // magnetization
    stima_m = (blk_av[im] / blk_norm) / (double) nspin;
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m * stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);

    // suscieptivity
    stima_x = beta * blk_av[ix] / blk_norm / (double) nspin;
    glob_av[ix] += stima_g;
    glob_av2[ix] += stima_g * stima_g;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);

    // printing results on files
    Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double) iblk << setw(wd) << err_u
        << endl;
    Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double) iblk << setw(wd) << err_c
         << endl;
    Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im] / (double) iblk << setw(wd) << err_m
        << endl;
    Chi << setw(wd) << iblk << setw(wd) << stima_g << setw(wd) << glob_av[ix] / (double) iblk << setw(wd) << err_g
        << endl;

    // closing files
    Ene.close();
    Heat.close();
    Mag.close();
    Chi.close();
}


void ConfFinal(void) {
    ofstream WriteConf;

    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");
    for (int i = 0; i < nspin; ++i) {
        WriteConf << s[i] << endl;
    }
    WriteConf.close();

    rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if (i >= nspin) i = i - nspin;
    else if (i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk) {
    if (iblk == 1) return 0.0;
    else return sqrt((sum2 / (double) iblk - pow(sum / (double) iblk, 2)) / (double) (iblk - 1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
