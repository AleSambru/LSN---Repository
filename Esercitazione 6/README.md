### Output Data: 
# Exercise 6
make Monte_Carlo_ISING.exe -> compile Monte_Carlo_ISING.exe
./Monte_Carlo_ISING.exe -> execute 

# input files
- input.dat

### output data 
in ISING_1D are containend data relative to the equilibration, in particulair we consider the following outputs: 
- ISING_1D/output.metro.eq.ene.h=0.000000t=1.000000
- ISING_1D/output.metro.eq.ene.h=0.000000t=1.500000
- ISING_1D/output.metro.eq.ene.h=0.000000t=2.000000
- ISING_1D/output.metro.eq.ene.h=2.000000t=1.000000
- ISING_1D/output.metro.eq.ene.h=2.000000t=1.500000
- ISING_1D/output.metro.eq.ene.h=2.000000t=2.000000

- ISING_1D/output.gibbs.eq.ene.h=0.000000t=1.000000
- ISING_1D/output.gibbs.eq.ene.h=0.000000t=1.500000
- ISING_1D/output.gibbs.eq.ene.h=0.000000t=2.000000
- ISING_1D/output.gibbs.eq.ene.h=2.000000t=1.000000
- ISING_1D/output.gibbs.eq.ene.h=2.000000t=1.500000
- ISING_1D/output.gibbs.eq.ene.h=2.000000t=2.000000
Related to the different starting condition of the system in order to define the appropriate n of equilibrium steps. 
In this folder other output data are created but we are not going to consider them in the python code 

In the main folders the following data are used 
output.chi.0.020000gibbs.dat
output.chi.0.020000metro.dat
output.ene.0.020000gibbs.dat
output.ene.0.020000metro.dat
output.heat.0.020000gibbs.dat
output.heat.0.020000metro.dat
output.mag.0.020000gibbs.dat
output.mag.0.020000metro.dat