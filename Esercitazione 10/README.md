# Exercise 10
To correctly compile the program, MPI must be correctly installed. These programs have been executed remotly connecting to the Uni Lab

The following guidelines must be followed.
- module load mpi/mpich-x86_64
- make 
- mpiexec -np n ./Esercizio_10.exe n
 
where n is the number of cores used 

output data: 
- In Results folder :
    - American_Capitals_Best_Path_1_core.dat
    - American_capitals_Results_1_core.dat
    - 2_Rank_Used
    - 4_Ranks_Used
    - 6_Ranks_Used
    - 8_Ranks_Used
    - 10_Ranks_Used
In each folder it's containend
- American_capitals_Results_n_RankUsed_m.dat
- BestPath_n_RanksUsed_m.dat
 where n correspond to the numbers of Ranks used , while m correspond to the specific ranks whose results are
