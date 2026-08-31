


#include <iostream>
#include "funcs.h"
#include <cmath>

using namespace std;
void compute_normalized_arrays(const double* U, const double* L, double* U_normal, double* L_normal, int length,double *max_mean,double *min_mean, double *max_std, double *min_std) {
    for (int i = 0; i < length; i++) {

        
        double numerator_U = U[i] - min_mean[i];

        
        if (numerator_U >= 0) {
            U_normal[i] = numerator_U / min_std[i];
        }
        
        else {
            U_normal[i] = numerator_U / max_std[i];
        }

        
        double numerator_L = L[i] - max_mean[i];

        
        if (numerator_L <= 0) {
            L_normal[i] = numerator_L / min_std[i];
        }
        
        else {
            L_normal[i] = numerator_L / max_std[i];
        }
    }
}
