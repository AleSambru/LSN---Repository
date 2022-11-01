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
#include <string>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;

int main() {

    // definition of mode and phase (=solid, liquid, gas)
    Program_SetUp();
    switch (code_mode) {// defining equilibration stap for the AR
        case 0:
            cout << "\n######################################################\n"
                    "Setting parameters for a proper acceptance rate = 50%\n"
                    "######################################################\n\n";
            cout << "Select n step = ";
            cin >> n_step_eq;
            cout << "\nMode -> 1 = solid, 2 = liquid, 3 = gas \n";
            int prel_mode;
            cin >> prel_mode;
            cout << "\n-----------------------------------------------\n"
                    "Printing out Epot values - n step =" << n_step_eq
                 << "\n-----------------------------------------------\n" << endl;
            Instant_Epot_values(prel_mode, n_step_eq);
            cout << "Acceptance rate  = " << accepted / attempted << endl << endl;
            cout << "###################################################################\n";
            break;
        case 1:// n steps -> output data
            n_step_eq = 50000;
            n_steps = 500000; // nsteps to calculate the correlation
            cout << "\n#####################\nAcceptance Rate has been set\n#####################\n\n";
            cout << "\n      Printing out Epot values - n step =" << n_step_eq << "      " << endl;
            cout << "\nSelect Mode -> 1 = solid, 2 = liquid, 3 = gas \n";
            cin >> mode;
            Instant_Epot_values(mode, n_step_eq);
            break;
    }

// end
    rnd.SaveSeed();
    return 0;
}

//================================================================================================
// function definition
void Program_SetUp() {
    cout << "\n#######################################################################\n";
    cout << "\n******************************************************\n";
    cout << "This program perform a MC simulation of a classic Lennard-Jones fluid\n";
    cout << "In the different state : solid, liquid, gas";
    cout << "\n******************************************************\n";
    cout << "\n#######################################################################\n";
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl;
    cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;
    cout << "Data are taken form input files\ninput.in -> gas\ninput.liquid-> liquid\ninput.solid -> solid\n";
    cout << "\n-----------------------------------------------------------------------\n";
    cout << "Program Operation:\n";

    cout << "The first part of the program is used to print out instantaneous values of potential energy per particle\n"
            "MC steps  = $5×10^5$\nthe Metropolis acceptance ratio has to be fixed approximately to 50%\n";


    cout << "Digit:\n0-> to do a preliminary study in order to fix the acceptance rate;"
            "\n1-> to print out M = $5×10^5$ instantaneous values of potential energy;";
    cin >> code_mode;
}

// define the input file
void SetUp(int mode) {
    switch (mode) {
        case 1:
            cout << "\n-----------\nSolid Case\n-----------" << endl;
            Input("input.solid");
            break;
        case 2:
            cout << "\n-----------\nLiquid Case\n-----------" << endl;
            Input("input.liquid");
            break;
        case 3:
            cout << "\n-----------\nGas Case\n-----------" << endl;
            Input("input.in");
            break;
    }
}

// setup variable taking info from a specific file
void Input(string input_file_name) {

    ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

// Random Generator Setup
    //Read seed for random numbers
    int p1, p2;
    Primes.open("Primes");
    Primes >> p1 >> p2;
    Primes.close();


    //Read input informations
    ReadInput.open(input_file_name);

    ReadInput >> iNVET;
    ReadInput >> restart;

    if (iNVET == 0)
        cout << "----------------------------\nPerforming Molecolar Dynamic\n----------------------------\n";
    else if (iNVET == 1)
        cout << "-------------------------------\nPerforming Montecarlo simulation\n-------------------------------\n";;
    cout << "Parameters : \n\n";
    if (restart) Seed.open("seed.out");
    else Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed, p1, p2);
    Seed.close();

    ReadInput >> temp;
    beta = 1.0 / temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double) npart / rho;
    box = pow(vol, 1.0 / 3.0);
    cout << "Volume of the simulation box = " << vol << endl;
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;

    ReadInput >> delta;

    cout << "delta = " << delta << endl;

    ReadInput >> nblk;

    ReadInput >> nstep;

    //cout << "The program perform Metropolis moves with uniform translations" << endl;
    //cout << "Moves parameter = " << delta << endl;
    //cout << "Number of blocks = " << nblk << endl;
    //cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();

//Prepare arrays for measurements
    iv = 0; //Potential energy
    it = 1; //Temperature
    ik = 2; //Kinetic energy
    ie = 3; //Total energy
    ip = 4; // pressure
    n_props = 5; //Number of observables

//Read initial configuration
    cout << "Read initial configuration" << endl << endl;
    if (restart) {
        ReadConf.open("config.out");
        ReadVelocity.open("velocity.out");
        for (int i = 0; i < npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
    } else {
        ReadConf.open("config.in");
        cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < npart; ++i) {
            vx[i] = rnd.Gauss(0., sqrt(temp));
            vy[i] = rnd.Gauss(0., sqrt(temp));
            vz[i] = rnd.Gauss(0., sqrt(temp));
            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        for (int idim = 0; idim < 3; ++idim) sumv[idim] /= (double) npart;
        double sumv2 = 0.0, fs;
        for (int i = 0; i < npart; ++i) {
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];
            sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        }
        sumv2 /= (double) npart;
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        cout << "velocity scale factor: " << fs << endl << endl;
        for (int i = 0; i < npart; ++i) {
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;
        }
    }

    for (int i = 0; i < npart; ++i) {
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = Pbc(x[i] * box);
        y[i] = Pbc(y[i] * box);
        z[i] = Pbc(z[i] * box);
    }
    ReadConf.close();

    for (int i = 0; i < npart; ++i) {
        if (iNVET) {
            xold[i] = x[i];
            yold[i] = y[i];
            zold[i] = z[i];
        } else {
            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    }

    // tale correction
    // aggiunta della correzione per r> rcut
    // correzione sul potenziale
    vtail = 8 * M_PI * rho / 3 * (pow(3, -1) * pow(rcut, -9) - pow(rcut, -3)); // g(r) == 1
    ptail = 32 * M_PI * rho * (pow(6, -1) * pow(rcut, -9) - pow(6, -1) * pow(rcut, -3));// g(r) == 1
    // g(r) != 1 -> importance sampling
    cout << "Tail corrections:\n";
    cout << "<V> tail      = " << vtail << endl;
    cout << "<p> tail   = " << ptail << endl;

//Evaluate properties of the initial configuration
    Measure();

//Print initial values for measured properties
    cout << "Initial potential energy = " << walker[iv] / (double) npart << endl;
    cout << "Initial pressure         = " << walker[ip] << endl;

    return;
}

// Move the Algorithm -> MC/MD
void Move() {
    int o;
    double p, energy_old, energy_new;
    double xnew, ynew, znew;

    if (iNVET) // Monte Carlo (NVT) move
    {
        for (int i = 0; i < npart; ++i) {
            //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
            o = (int) (rnd.Rannyu() * npart);

            //Old
            energy_old = Boltzmann(x[o], y[o], z[o], o);

            //New
            x[o] = Pbc(x[o] + delta * (rnd.Rannyu() - 0.5));
            y[o] = Pbc(y[o] + delta * (rnd.Rannyu() - 0.5));
            z[o] = Pbc(z[o] + delta * (rnd.Rannyu() - 0.5));

            energy_new = Boltzmann(x[o], y[o], z[o], o);

            //Metropolis test
            p = exp(beta * (energy_old - energy_new));
            if (p >= rnd.Rannyu()) {
                //Update
                xold[o] = x[o];
                yold[o] = y[o];
                zold[o] = z[o];
                accepted = accepted + 1.0;
            } else {
                x[o] = xold[o];
                y[o] = yold[o];
                z[o] = zold[o];
            }
            attempted = attempted + 1.0;
        }
    } else // Molecular Dynamics (NVE) move
    {
        double fx[m_part], fy[m_part], fz[m_part];

        for (int i = 0; i < npart; ++i) { //Force acting on particle i
            fx[i] = Force(i, 0);
            fy[i] = Force(i, 1);
            fz[i] = Force(i, 2);
        }

        for (int i = 0; i < npart; ++i) { //Verlet integration scheme

            xnew = Pbc(2.0 * x[i] - xold[i] + fx[i] * pow(delta, 2));
            ynew = Pbc(2.0 * y[i] - yold[i] + fy[i] * pow(delta, 2));
            znew = Pbc(2.0 * z[i] - zold[i] + fz[i] * pow(delta, 2));

            vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
            vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
            vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);

            xold[i] = x[i];
            yold[i] = y[i];
            zold[i] = z[i];

            x[i] = xnew;
            y[i] = ynew;
            z[i] = znew;

            accepted = accepted + 1.0;
            attempted = attempted + 1.0;
        }
    }
    return;
}

// calculate boltzmann energy
double Boltzmann(double xx, double yy, double zz, int ip) {
    double ene = 0.0;
    double dx, dy, dz, dr;

    for (int i = 0; i < npart; ++i) {
        if (i != ip) {
// distance ip-i in pbc
            dx = Pbc(xx - x[i]);
            dy = Pbc(yy - y[i]);
            dz = Pbc(zz - z[i]);

            dr = dx * dx + dy * dy + dz * dz;
            dr = sqrt(dr);

            if (dr < rcut) {
                ene += 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
            }
        }
    }

    return 4.0 * ene;
}

double Force(int ip, int idir) { //Compute forces as -Grad_ip V(r)
    double f = 0.0;
    double dvec[3], dr;

    for (int i = 0; i < npart; ++i) {
        if (i != ip) {
            dvec[0] = Pbc(x[ip] - x[i]);  // distance ip-i in pbc
            dvec[1] = Pbc(y[ip] - y[i]);
            dvec[2] = Pbc(z[ip] - z[i]);

            dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
            dr = sqrt(dr);

            if (dr < rcut) {
                f += dvec[idir] * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8)); // -Grad_ip V(r)
            }
        }
    }

    return f;
}

void Measure() {//Properties measurement
    //cout << "Start measuring\n" ;
    double v = 0.0, kin = 0.0, p = 0.0;
    double vij, pij;
    double dx, dy, dz, dr;

//cycle over pairs of particles
    //cout << "cycle over pairs of particles\n" ;
    for (int i = 0; i < npart - 1; ++i) {
        for (int j = i + 1; j < npart; ++j) {
// distance i-j in pbc
            dx = Pbc(x[i] - x[j]);
            dy = Pbc(y[i] - y[j]);
            dz = Pbc(z[i] - z[j]);


            dr = dx * dx + dy * dy + dz * dz;
            dr = sqrt(dr);

            if (dr < rcut) {
                vij = 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
                pij = 1.0 / pow(dr, 12) - 0.5 / pow(dr, 6); // accumulatore per p
                v += vij;
                p += pij;
            }


            //int index = floor(dr/bin_size); //Rounds x downward, returning the largest integral value that is not greater than x.
            //if (dr < box/2){
            //    bins[index] += 2;
            //}
        }
    }

    //cout << "Calculating kinetic energy\n" ;
    for (int i = 0; i < npart; ++i) kin += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

    walker[iv] = 4.0 * v; // Potential energy
    walker[ik] = kin; // Kinetic energy
    walker[it] = (2.0 / 3.0) * kin / (double) npart; // Temperature
    walker[ie] = 4.0 * v + kin;  // Total energy;
    walker[ip] = rho * kB * walker[it] + 16.0 * p / (vol * double(npart)); // pressure

    return;
}


void Reset(int iblk) //Reset block averages
{

    if (iblk == 1) {
        for (int i = 0; i < n_props; ++i) {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for (int i = 0; i < nbins; i++) {
        glob_bin[i] = 0;
        glob_bin2[i] = 0;
    }

    for (int i = 0; i < n_props; ++i) {
        blk_av[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}


void Accumulate(void) //Update block averages
{

    for (int i = 0; i < n_props; ++i) {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;

    //Accumulo i bins
    for (int i = 0; i < nbins; i++) {
        bin_av[i] += bins[i];
    }
}

void Averages(int iblk) //Print results for current block
{
    ofstream Epot, Pres;
    const int wd = 12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted / attempted << endl << endl;

    Epot.open("output_epot.dat", ios::app);
    Pres.open("output_pres.dat", ios::app);

    stima_pot = (blk_av[iv]) / blk_norm / (double) npart + vtail; //Potential energy-> add tale correction
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot * stima_pot;
    err_pot = Error(glob_av[iv], glob_av2[iv], iblk);

    stima_pres = (blk_av[ip]) / blk_norm + ptail; // tale correction
    glob_av[ip] += stima_pres;
    glob_av2[ip] += stima_pres * stima_pres;
    err_press = Error(glob_av[ip], glob_av2[ip], iblk);

    Epot.close();
    Pres.close();

}

void Instant_Epot_values(int mode, int nstep_eq) {

    SetUp(mode);
    ofstream pre_pot, pre_pres, pot_corr, pres_corr;

    switch (mode) {
        case 1:// solid
            pre_pot.open("pre_pot_solid.dat", ios::app);
            pre_pres.open("pre_pres_solid.dat", ios::app);

            pot_corr.open("pot_solid_corr.dat", ios::app);
            pres_corr.open("pres_solid_corr.dat", ios::app);
            break;
        case 2:// liquid
            pre_pot.open("pre_pot_liquid.dat", ios::app);
            pre_pres.open("pre_pres_liquid.dat", ios::app);

            pot_corr.open("pot_liquid_corr.dat", ios::app);
            pres_corr.open("pres_liquid_corr.dat", ios::app);
            break;
        case 3:// gas
            pre_pot.open("pre_pot_gas.dat", ios::app);
            pre_pres.open("pre_pres_gas.dat", ios::app);

            pot_corr.open("pot_gas_corr.dat", ios::app);
            pres_corr.open("pres_gas_corr.dat", ios::app);
            break;
    }

    cout << "---------------------\nPotential Energy values\n---------------------\n";
    // Testing Equilibration
    nstep_eq = 50000;
    cout << "Equilibration\n=====================\n";
    for (int i = 0; i < nstep_eq; i++) //Equilibrazione
    {
        Move();
        Measure();

        stima_pot = walker[iv] / npart + vtail; //Energia
        stima_pres = walker[ip] + ptail;
        pre_pot << setw(wd) << i << setw(wd) << stima_pot << endl;
        pre_pres << setw(wd) << i << setw(wd) << stima_pres << endl;

        if (i % 100 == 0) {
            cout << "Step " << i << " - acceptance rate = " << double(accepted) / attempted << " ==>  E pot = "
                 << stima_pot << " ==> " << int(100 * (double(i) / nstep_eq)) << " % \r    " << flush;
        }
    }
    //cout << "Final Acceptance rate = " << double (accepted)/attempted << endl ;

    pre_pres.close();
    pre_pot.close();

    if (code_mode == 1) {
        cout << "Printing Instantaneous values\n==============================\n";
        for (int i = 0; i < n_steps; i++) //Equilibrazione
        {

            Move();
            Measure();

            stima_pot = walker[iv] / npart + vtail; //Energia
            stima_pres = walker[ip] + ptail;
            pot_corr << i << setw(wd) << stima_pot << endl;
            pres_corr << i << setw(wd) << stima_pres << endl;

            if (i % 100 == 0) {
                cout << "Step " << i << " - acceptance rate = " << double(accepted) / attempted << " ==>  E pot = "
                     << stima_pot << " ==> " << int(100 * (double(i) / nstep_eq)) << " % \r    " << flush;
            }

        }
    }
};


void ConfFinal(void) {
    ofstream WriteConf, WriteVelocity, WriteSeed;

    cout << "Print final configuration to file config.out" << endl << endl;
    WriteConf.open("config.out");
    WriteVelocity.open("velocity.out");
    for (int i = 0; i < npart; ++i) {
        WriteConf << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
        WriteVelocity << vx[i] << "   " << vy[i] << "   " << vz[i] << endl;
    }
    WriteConf.close();
    WriteVelocity.close();

    rnd.SaveSeed();
}

void ConfXYZ(int nconf) { //Write configuration in .xyz format
    ofstream WriteXYZ;

    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for (int i = 0; i < npart; ++i) {
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " << Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
    }
    WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r / box);
}

double Error(double sum, double sum2, int iblk) {
    return sqrt(fabs(sum2 / (double) iblk - pow(sum / (double) iblk, 2)) / (double) iblk);
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
