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

using namespace std;


//Function for statistical uncertainty estimation:
/*Description of the punction: the function estimate the uncertainty between 

Inputs:
AV: vector containing the cumulative averages 
AV2: vector containig the cumulative averages^2
n: number of block 
*/

template<typename T>
double error(vector <T> AV, vector <T> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt((AV2[n] - pow(AV[n], 2) )/n);
    };
};




#endif

