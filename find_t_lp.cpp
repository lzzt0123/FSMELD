


#include <cstdio>
#include <stdlib.h>
#include "funcs.h"
#include <algorithm>
#include <queue>
#include <cmath>
#include <string.h>
#include <cstring>
#include <execution>

using namespace std;












struct HeapItem {
    double value;
    int index;

    bool operator<(const HeapItem& other) const {
        return value > other.value;
    }
};


struct maxHeapItem {
    double value;
    int index;

    bool operator<(const maxHeapItem& other) const {
        return value < other.value;
}

};

struct UniqueItem {
    double value;
    int count;
    bool operator<(const UniqueItem& other) const {
        
        return value < other.value;
    }
};


int create_unique_struct_array(double *T, UniqueItem *output, int m,int n) {
    if (n <= 0) return 0;
    int j = 0;
    
    for (int i = 0; i < n; i++) {
        if (i == 0 || T[m-1+i] != T[m+i-2]) {
            output[j].value = T[m-1+i];
            output[j].count = 1;
            j++;
        } else {
            output[j-1].count++;
        }
    }
    return j;
}
void get_sorted_struct_array(UniqueItem *T_sort, int n, int isupper) {

    if (isupper == 1) {
        
        std::sort(T_sort, T_sort + n, [](const UniqueItem& a, const UniqueItem& b) {
            return a.value < b.value;
        });
    } else {
        
        std::sort(T_sort, T_sort + n, [](const UniqueItem& a, const UniqueItem& b) {
            return a.value > b.value;
        });
    }
}


void get_sorted_array(double* T_sort, int n, int isupper) {

    if (isupper == 1) {
        
        std::sort(T_sort, T_sort + n);
    } else {
        
        
        std::sort(T_sort, T_sort + n, std::greater<double>());
    }
}
int create_unique_array(double *T, double *T_remove, int *counts, int n) {
    
    if (n <= 0) {
        return 0;
    }

    int j = 0;

    for (int i = 0; i < n; i++) {
        
        if (i == 0 || T[i] != T[i-1]) {
            T_remove[j] = T[i]; 
            counts[j] = 1;      
            j++;                
        }
        else {
            
            
            counts[j-1]++;
        }
    }

    return j;
}
void sort_parallel_arrays(double* values, int* counts, int n) {
    if (n <= 0) return;

    int* p = new int[n];

    for (int i = 0; i < n; ++i) {
        p[i] = i;
    }

    std::sort(p, p + n, [values](int a, int b) {
        return values[a] < values[b];
    });

    double* temp_values = new double[n];
    int* temp_counts = new int[n];

    for (int i = 0; i < n; ++i) {
        
        int original_index = p[i];
        temp_values[i] = values[original_index];
        temp_counts[i] = counts[original_index];
    }

    
    for (int i = 0; i < n; ++i) {
        values[i] = temp_values[i];
        counts[i] = temp_counts[i];
    }

    delete[] p;
    delete[] temp_values;
    delete[] temp_counts;
}
double D(double q_i,double C,const char *tag) {
    if(strcmp(tag,"msm")==0){
        return C;
    }else if(strcmp(tag,"erp")==0) {
        return (q_i-C)*(q_i-C);
    }else if (strcmp(tag,"edr")==0||strcmp(tag,"lcss")==0) {
        return 1;
    }else if(strcmp(tag,"dtw")==0) {
        return (q_i-C)*(q_i-C);
    }
    return 0.0;
}

double t_u(double q_i,double C,const char *tag) {
    if(strcmp(tag,"msm")==0){
        return q_i-C;
    }else if(strcmp(tag,"erp")==0) {
        return q_i-abs(q_i-C);
    }else if (strcmp(tag,"edr")==0||strcmp(tag,"lcss")==0) {
        return q_i-C;
    }
    return q_i;
}
double t_l(double q_i,double C,const char *tag) {
    if(strcmp(tag,"msm")==0){
        return q_i+C;
    }else if(strcmp(tag,"erp")==0) {
        return q_i+abs(q_i-C);
    }else if (strcmp(tag,"edr")==0||strcmp(tag,"lcss")==0) {
        return q_i+C;
    }
    return q_i;
}
double f_p(double *m_p, int p, const char *tag, double t) {
    if (p == 0) {
        if(strcmp(tag,"lcss")==0||strcmp(tag,"edr")==0) {
            return 0;
        }
        return m_p[0];
    }else if (p == 1) {
        return m_p[1] - m_p[0] * t;
    }else if (p == 2) {
        return m_p[2]-m_p[1]*t*2+m_p[0]*pow(t,2);
    }
    return 0.0;
}
int findMaxIndex(const double *arr, int size) {
    if (size <= 0) return -1;  

    int maxIndex = 0;
    for (int i = 0; i < size; ++i) {
        if (arr[i] > arr[maxIndex]) {
            maxIndex = i;
        }
    }
    return maxIndex;
}

double find_tl_lp_dtw(double *L, int l, double *q, int m, int p, double *LB, double *m_p, const char* tag, const double* q_sorted) {
    double tl = 0.0, h_L = 0.0;
    int g_L = 0;
    double t_0;
    int n_remove;

    UniqueItem* L_remove = new UniqueItem[l - 2 * m + 2];

    memset(m_p, 0, sizeof(double) * (p + 1));
    memset(LB, 0, sizeof(double) * (l - 2 * m + 2));

    n_remove = create_unique_struct_array(L, L_remove, m, l - 2 * m + 2);

    
    get_sorted_struct_array(L_remove, n_remove, 0);

    g_L = 0;
    t_0 = L_remove[0].value;

    
    double* q_copy_owned = nullptr;
    const double* q_copy = q_sorted;
    if (q_copy == nullptr) {
        q_copy_owned = new double[m];
        std::memcpy(q_copy_owned, q, m * sizeof(double));
        std::sort(q_copy_owned, q_copy_owned + m);
        q_copy = q_copy_owned;
    }

    int idx1 = m - 1;

    
    
    
    while (idx1 >= 0 && q_copy[idx1] >= t_0) {
        idx1--;
    }

    
    
    for (int i = 0; i <= idx1; i++) {
        for (int a = 0; a < p + 1; a++) {
            m_p[a] += std::pow(q_copy[i], a);
        }
    }

    
    for (int i = 0; i < n_remove; i++) {
        double t = L_remove[i].value;
        g_L = g_L + L_remove[i].count;

        
        
        while (idx1 >= 0) {
            double val = q_copy[idx1];

            
            if (val >= t) {
                for (int a = 0; a < p + 1; a++) {
                    m_p[a] -= std::pow(val, a);
                }
                idx1--; 
            } else {
                
                break;
            }
        }

        
        h_L = std::abs(f_p(m_p, p, tag, t));
        LB[i] = h_L * g_L;
    }

    int index = findMaxIndex(LB, n_remove);
    tl = L_remove[index].value;

    
    delete[] q_copy_owned;
    delete[] L_remove;

    
    return tl;
}
double find_tl_lp(double *L,int l,double *q,int m,int p,double C_,double *LB,double *m_p,const char* tag,const double* q_sorted){

    double tl=0.0,h_L=0.0,C;
    int g_L;
    double t_0;

    
    int n_remove;

    memset(m_p, 0, sizeof(double) * (p + 1));
    memset(LB, 0, sizeof(double) * (l-2*m+2));

    UniqueItem* L_remove = new UniqueItem[l-2*m+2];

    n_remove=create_unique_struct_array(L,L_remove,m,l-2*m+2);

    get_sorted_struct_array(L_remove,n_remove,0);


    g_L=0;
    C=0;
    t_0=L_remove[0].value;

    
    
    

    int idx1 = m - 1;
    
    double* q_copy_owned = nullptr;
    const double* q_copy = q_sorted;
    if (q_copy == nullptr) {
        q_copy_owned = new double[m];
        std::memcpy(q_copy_owned, q, m * sizeof(double));
        std::sort(q_copy_owned, q_copy_owned + m);
        q_copy = q_copy_owned;
    }

    while (idx1 >= 0 && q_copy[idx1] >= t_0) {
        idx1--;
    }

    
    int idx2 = idx1;

    
    while (idx2 >= 0) {
        double val = t_l(q_copy[idx2],C_,tag);

        if (val >= t_0) {
            
            for (int a = 0; a < p + 1; a++) {
                m_p[a] += std::pow(q_copy[idx2], a);
            }
            idx2--;
        } else {
            
            
            break;
        }
    }

    
    for (int j = 0; j <= idx2; j++) {
        C += D(q_copy[j], C_,tag);
    }

    for (int i = 0; i < n_remove; i++) {
        double t = L_remove[i].value;
        g_L = g_L + L_remove[i].count;

        
        
        
        while (idx2 >= 0) {
            double val = t_l(q_copy[idx2],C_,tag);
            
            
            if (val >= t) {
                
                C = C - D(q_copy[idx2], C_,tag);

                
                for (int a = 0; a < p + 1; a++) {
                    m_p[a] += std::pow(q_copy[idx2], a);
                }

                idx2--;
            } else {
                
                break;
            }
        }

        
        
        
        
        while (idx1 > idx2) {
            double val = q_copy[idx1];
            if (val > t) {
                
                for (int a = 0; a < p + 1; a++) {
                    m_p[a] -= std::pow(val, a);
                }
                
                idx1--;
            } else {
                
                break;
            }
        }

        
        h_L = C + std::abs(f_p(m_p, p, tag,t));
        LB[i] = h_L * g_L;
    }

    int index = findMaxIndex(LB,n_remove);
    tl=L_remove[index].value;
    
    delete[] q_copy_owned;

    delete[] L_remove;
    
    return tl;
}
double find_tl_lp_erp(double *L, int l, double *q, int m, int p, double C_, double *LB, double *m_p) {
    double tl = 0.0, h_L = 0.0, C = 0.0;
    int g_L = 0;
    double t_0;
    int n_remove;

    UniqueItem* L_remove = new UniqueItem[l - 2 * m + 2];
    memset(m_p, 0, sizeof(double) * (p + 1));
    memset(LB, 0, sizeof(double) * (l - 2 * m + 2));

    n_remove = create_unique_struct_array(L, L_remove, m, l - 2 * m + 2);
    
    get_sorted_struct_array(L_remove, n_remove, 0);
    t_0 = L_remove[0].value;

    std::priority_queue<maxHeapItem> H2;
    std::priority_queue<maxHeapItem> H1;

    
    for (int i = 0; i < m; i++) {
        double tau = t_l(q[i], C_, "erp");

        if (q[i] >= t_0) {
            
        } else if (tau >= t_0) {
            
            H1.push({q[i], i});
            for (int a = 0; a < p + 1; a++) {
                m_p[a] += std::pow(q[i], a);
            }
        } else {
            
            H2.push({tau, i});
            C += D(q[i], C_, "erp");
        }
    }

    
    for (int i = 0; i < n_remove; i++) {
        double t = L_remove[i].value;
        g_L = g_L + L_remove[i].count;

        
        while (!H2.empty() && H2.top().value >= t) {
            int orig_idx = H2.top().index;
            double orig_q = q[orig_idx];
            H2.pop();

            C -= D(orig_q, C_, "erp");
            H1.push({orig_q, orig_idx});

            for (int a = 0; a < p + 1; a++) {
                m_p[a] += std::pow(orig_q, a);
            }
        }

        
        while (!H1.empty() && H1.top().value >= t) {
            int orig_idx = H1.top().index;
            double orig_q = q[orig_idx];
            H1.pop();

            for (int a = 0; a < p + 1; a++) {
                m_p[a] -= std::pow(orig_q, a);
            }
        }

        h_L = C + std::abs(f_p(m_p, p, "erp", t));
        LB[i] = h_L * g_L;
    }

    int index = findMaxIndex(LB, n_remove);
    tl = L_remove[index].value;

    delete[] L_remove;
    return tl;
}
double find_tu_lp_erp(double *U, int l, double *q, int m, int p, double C_, double *LB, double *m_p) {
    double tu = 0.0, h_U = 0.0, C = 0.0;
    int g_U = 0;
    double t_0;
    int n_remove;

    UniqueItem* U_remove = new UniqueItem[l - 2 * m + 2];
    memset(m_p, 0, sizeof(double) * (p + 1));
    memset(LB, 0, sizeof(double) * (l - 2 * m + 2));

    n_remove = create_unique_struct_array(U, U_remove, m, l - 2 * m + 2);
    get_sorted_struct_array(U_remove, n_remove, 1);
    t_0 = U_remove[0].value;

    
    std::priority_queue<HeapItem> H2; 
    std::priority_queue<HeapItem> H1; 

    
    for (int i = 0; i < m; i++) {
        double tau = t_u(q[i], C_, "erp");

        if (q[i] <= t_0) {
            
        } else if (tau <= t_0) {
            
            H1.push({q[i], i});
            for (int a = 0; a < p + 1; a++) {
                m_p[a] += std::pow(q[i], a);
            }
        } else {
            
            H2.push({tau, i});
            C += D(q[i], C_, "erp");
        }
    }

    
    for (int i = 0; i < n_remove; i++) {
        double t = U_remove[i].value;
        g_U = g_U + U_remove[i].count;

        
        while (!H2.empty() && H2.top().value <= t) {
            int orig_idx = H2.top().index; 
            double orig_q = q[orig_idx];   
            H2.pop();

            C -= D(orig_q, C_, "erp");
            H1.push({orig_q, orig_idx});   

            for (int a = 0; a < p + 1; a++) {
                m_p[a] += std::pow(orig_q, a);
            }
        }

        
        while (!H1.empty() && H1.top().value <= t) {
            int orig_idx = H1.top().index;
            double orig_q = q[orig_idx];
            H1.pop();

            for (int a = 0; a < p + 1; a++) {
                m_p[a] -= std::pow(orig_q, a);
            }
        }

        
        h_U = C + f_p(m_p, p, "erp", t);
        LB[i] = h_U * g_U;
    }

    int index = findMaxIndex(LB, n_remove);
    tu = U_remove[index].value;

    delete[] U_remove;
    return tu;
}
double find_tu_lp_dtw(double *U, int l, double *q, int m, int p, double *LB, double *m_p, const char* tag, const double* q_sorted) {
    double tu = 0.0, h_U = 0.0;
    int g_U = 0;
    double t_0;
    int n_remove;

    UniqueItem* U_remove = new UniqueItem[l - 2 * m + 2];

    memset(m_p, 0, sizeof(double) * (p + 1));
    memset(LB, 0, sizeof(double) * (l - 2 * m + 2));

    n_remove = create_unique_struct_array(U, U_remove, m, l - 2 * m + 2);
    get_sorted_struct_array(U_remove, n_remove, 1);

    g_U = 0;
    t_0 = U_remove[0].value;

    
    double* q_copy_owned = nullptr;
    const double* q_copy = q_sorted;
    if (q_copy == nullptr) {
        q_copy_owned = new double[m];
        std::memcpy(q_copy_owned, q, m * sizeof(double));
        std::sort(q_copy_owned, q_copy_owned + m);
        q_copy = q_copy_owned;
    }

    int idx1 = 0;

    
    while (idx1 < m && q_copy[idx1] <= t_0) {
        idx1++;
    }

    
    
    for (int i = idx1; i < m; i++) {
        for (int a = 0; a < p + 1; a++) {
            m_p[a] += std::pow(q_copy[i], a);
        }
    }

    
    for (int i = 0; i < n_remove; i++) {
        double t = U_remove[i].value;
        g_U = g_U + U_remove[i].count;

        
        
        while (idx1 < m) {
            double val = q_copy[idx1];

            
            if (val < t) {
                
                for (int a = 0; a < p + 1; a++) {
                    m_p[a] -= std::pow(val, a);
                }
                idx1++;
            } else {
                
                break;
            }
        }

        h_U = f_p(m_p, p, tag, t);
        LB[i] = h_U * g_U;
    }

    int index = findMaxIndex(LB, n_remove);
    tu = U_remove[index].value;

    
    delete[] q_copy_owned;
    delete[] U_remove;

    
    return tu;
}

double find_tu_lp(double *U,int l,double *q,int m,int p,double C_,double *LB,double *m_p,const char* tag,const double* q_sorted) {
    double tu=0.0,h_U=0.0,C;
    int g_U;
    double t_0;
    int n_remove;

    UniqueItem* U_remove = new UniqueItem[l-2*m+2];

    memset(m_p, 0, sizeof(double) * (p + 1));
    memset(LB, 0, sizeof(double) * (l-2*m+2));

    n_remove=create_unique_struct_array(U,U_remove,m,l-2*m+2);

    get_sorted_struct_array(U_remove,n_remove,1);

    g_U=0;
    C=0;
    t_0=U_remove[0].value;

    
    double* q_copy_owned = nullptr;
    const double* q_copy = q_sorted;
    if (q_copy == nullptr) {
        q_copy_owned = new double[m];
        std::memcpy(q_copy_owned, q, m * sizeof(double));
        std::sort(q_copy_owned, q_copy_owned + m);
        q_copy = q_copy_owned;
    }


    int idx1 = 0; 
    int idx2 = 0; 

    
    while (idx1 < m && q_copy[idx1] <= t_0) {
        idx1++;
    }

    
    idx2 = idx1;

    
    while (idx2 < m) {
        double val = t_u(q_copy[idx2],C_,tag);
        
        
        if (val <= t_0) {
            for (int a = 0; a < p + 1; a++) {
                m_p[a] += std::pow(q_copy[idx2], a);
            }
            idx2++;
        } else {
            
            break;
        }
    }

    
    
    for (int j = idx2; j < m; j++) {
        C += D(q_copy[j], C_,tag);
    }
    for (int i = 0; i < n_remove; i++) {
        double t = U_remove[i].value;
        g_U = g_U + U_remove[i].count;

        
        while (idx2 < m) {
            double val = t_u(q_copy[idx2],C_,tag);
            
            if (val <= t) {

                C = C - D(q_copy[idx2], C_, tag);
                for (int a = 0; a < p + 1; a++) {
                    m_p[a] += std::pow(q_copy[idx2], a);
                }

                idx2++;
            } else {
                break;
            }
        }

        
        
        
        
        while (idx1 < idx2) {
            double val = q_copy[idx1];
            
            if (val < t) {
                for (int a = 0; a < p + 1; a++) {
                    m_p[a] -= std::pow(val, a);
                }
                idx1++;
            } else {
                break;
            }
        }

        h_U = C + f_p(m_p, p, tag,t);
        LB[i] = h_U * g_U;
    }

    int index = findMaxIndex(LB,n_remove);
    tu=U_remove[index].value;
    
    delete[] q_copy_owned;
    delete[] U_remove;
    
    return tu;
}

