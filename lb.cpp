


#include <iostream>
#include <string.h>
#include "funcs.h"


using namespace std;

void s_to_g_l(const double* L, int L_size, double t_l, double* g_l, int g_size) {
    s_to_g_l(L, L_size, t_l, g_l, g_size, nullptr, 0.0);
}

void s_to_g_u(const double* U, int U_size, double t_u, double* g_u, int g_size) {
    s_to_g_u(U, U_size, t_u, g_u, g_size, nullptr, 0.0);
}










void s_to_g_l(const double* L, int L_size, double t_l, double* g_l, int g_size, const char* tag, double C) {
    if (L == nullptr || g_l == nullptr || L_size != g_size) {
        return; 
    }
    for (int i = 0; i < L_size; i++) {
        g_l[i] = (L[i] >= t_l) ? 1.0 : 0.0;
    }
}










void s_to_g_u(const double* U, int U_size, double t_u, double* g_u, int g_size, const char* tag, double C) {
    if (U == nullptr || g_u == nullptr || U_size != g_size) {
        return; 
    }
    for (int i = 0; i < U_size; i++) {
        g_u[i] = (U[i] <= t_u) ? 1.0 : 0.0;
    }
}
