//
// Created by LZT on 2025-07-27.
//
#include <iostream>
#include "funcs.h"
#include "string.h"
#include <algorithm>

using namespace std;

void q_to_h_l(const double* q, int q_size, double t_l, double C, double* h_l, int h_size, const char* tag) {
    if (q == nullptr || h_l == nullptr || q_size != h_size) {
        return;
    }
    if(strcmp("msm",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if (q[i] < t_l) {
                h_l[i] = std::min(C, t_l - q[i]);
            } else {
                
                h_l[i] = 0.0;
            }
        }
    }else if(strcmp("dtw",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if(q[i]<t_l) {
                h_l[i] = (t_l-q[i])*(t_l-q[i]);
            }else {
                h_l[i] = 0.0;
            }

        }
    }
    else if(strcmp("erp",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if (q[i] < t_l) {
                if ((q[i]-t_l)*(q[i]-t_l) > (q[i]-C)*(q[i]-C)) {
                    h_l[i]=(q[i]-C)*(q[i]-C);
                }
                else {
                    h_l[i] = (q[i]-t_l)*(q[i]-t_l);
                }
            }
            else {
                h_l[i] = 0.0;
            }
        }

    }
    else if(strcmp("edr",tag)==0||strcmp("lcss",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if (q[i] < t_l-C) {
                h_l[i] = 1;
            }
            else {
                h_l[i] = 0.0;
            }
        }
    }
    else if(strcmp("swale",tag)==0) {
        q_to_h_l(q, q_size, t_l, C, h_l, h_size, tag,
                 SWALE_MISMATCH_COST, SWALE_MATCH_COST);
    }
}
void q_to_h_l(const double* q, int q_size, double t_l, double C, double* h_l, int h_size, const char* tag,double p,double r) {
    if (q == nullptr || h_l == nullptr || q_size != h_size) {
        return;
    }
    if(strcmp("swale",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if (q[i] < t_l-C) {
                h_l[i] = p;
            }
            else if(q[i]<t_l) {
                h_l[i] = r;
            }
            else {
                h_l[i] = 0.0;
            }
        }
    }
}

void q_to_h_u(const double* q, int q_size, double t_u, double C, double* h_u, int h_size,const char* tag) {

    if (q == nullptr || h_u == nullptr || q_size != h_size) {
        return;
    }
    if(strcmp("msm",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if (q[i] > t_u) {
                h_u[i] = std::min(C, q[i] - t_u);
            } else {
                h_u[i] = 0.0;
            }
        }
    }else if(strcmp("dtw",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if(q[i]>t_u) {
                h_u[i]=(q[i]-t_u)*(q[i]-t_u);
            }
            else {
                h_u[i] = 0.0;
            }
        }
    }else if(strcmp("edr",tag)==0||strcmp("lcss",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if (q[i] > t_u+C) {
                h_u[i] = 1;
            }
            else {
                h_u[i] = 0.0;
            }
        }
    }else if(strcmp("erp",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if(q[i]>t_u) {
                if((q[i]-t_u)*(q[i]-t_u)>(q[i]-C)*(q[i]-C)) {
                    h_u[i]=(q[i]-C)*(q[i]-C);
                }
                else {
                    h_u[i] = (q[i]-t_u)*(q[i]-t_u);
                }
            }
            else {
                h_u[i] = 0.0;
            }
        }
    }else if(strcmp("swale",tag)==0) {
        q_to_h_u(q, q_size, t_u, C, h_u, h_size, tag,
                 SWALE_MISMATCH_COST, SWALE_MATCH_COST);
    }

}
void q_to_h_u(const double* q, int q_size, double t_u, double C, double* h_u, int h_size,const char* tag,double p,double r) {
    if (q == nullptr || h_u == nullptr || q_size != h_size) {
        return;
    }
    if(strcmp("swale",tag)==0) {
        for (int i = 0; i < q_size; i++) {
            if (q[i] > t_u+C) {
                h_u[i] = p;
            }
            else if(q[i]>t_u) {
                h_u[i] = r;
            }
            else {
                h_u[i] = 0.0;
            }
        }
    }
}
