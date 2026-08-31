#include <algorithm>
#include <cmath>
#include <cstring> 
#include "def.h"
#include "funcs.h"
using namespace std;




inline void upper_lemire(double a[], int n, int r, double u[]) {
    int u_index = 0;
    int k = 2 * r + 1;
    int q_u[n];
    int head_u = 0, tail_u = -1;

    for (int i = 0; i < k - 1 - r; i++) {
        while (head_u <= tail_u && a[q_u[tail_u]] <= a[i]) tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r; i < n; i++) {
        while (head_u <= tail_u && a[q_u[tail_u]] <= a[i]) tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k) head_u++;
        u[u_index++] = a[q_u[head_u]];
    }

    for (int i = n; i < n + r; i++) {
        while (q_u[head_u] <= i - k) head_u++;
        u[u_index++] = a[q_u[head_u]];
    }
}
double edr_subcost(double a, double b, double epsilon) {
    return (abs(a - b) <= epsilon) ? 0.0 : 1.0;
}






double bglb_dtw(double* Q, int lenQ, double* T, int lenT, double* L_Q, double* U_Q, double* L_T, double* U_T, int w, double& bsf) {
    double bsf_sq = bsf * bsf;
    double d_boundary = dist(T[0], Q[0]) + dist(T[lenT-1], Q[lenQ-1]);
    double d_f1 = 0.0;

    double d1[lenT];
    memset(d1, 0, sizeof(double) * lenT);
    for (int i = 1; i < lenT - 1; i++) {
        if (T[i] < L_Q[i]) {
            d1[i] = dist(T[i], L_Q[i]);
        } else if (T[i] > U_Q[i]) {
            d1[i] = dist(T[i], U_Q[i]);
        }
        d_f1 += d1[i];
        if (d_boundary + d_f1 > bsf_sq) return MAX_VAL;
    }

    double d_b1 = 0.0;
    double d_f2 = 0.0;
    double d1_U[lenT];
    memset(d1_U, 0, sizeof(double) * lenT);
    upper_lemire(d1, lenT, w, d1_U);

    double d2[lenQ];
    memset(d2, 0, sizeof(double) * lenQ);
    for (int i = 1; i < lenQ - 1; i++) {
        if (Q[i] < L_T[i]) {
            d2[i] = dist(Q[i], L_T[i]);
        } else if (Q[i] > U_T[i]) {
            d2[i] = dist(Q[i], U_T[i]);
        }
        d_b1 += MAX(d2[i] - d1_U[i], 0.0);
        d_f2 += d2[i];
        if (d_boundary + d_f1 + d_b1 > bsf_sq || d_boundary + d_f2 > bsf_sq) return MAX_VAL;
    }

    double d_b2 = 0.0;
    double d2_U[lenQ];
    memset(d2_U, 0, sizeof(double) * lenQ);
    upper_lemire(d2, lenQ, w, d2_U);
    for (int i = 1; i < lenT - 1; i++) {
        d_b2 += MAX(d1[i] - d2_U[i], 0.0);
        if (d_boundary + d_f2 + d_b2 > bsf_sq) return MAX_VAL;
    }

    return sqrt(max(0.0, d_boundary + MAX(d_f1 + d_b1, d_f2 + d_b2)));
}


double bglb_edr(double* Q, int lenQ, double* T, int lenT, double epsilon, double* U_Q, double* L_Q, double* U_T, double* L_T, int w, double& bsf) {
    double d_boundary = 0.0;
    double d_f1 = 0.0;

    int pad_T = lenT + 2 * w;
    double d1[pad_T];
    memset(d1, 0, sizeof(double) * pad_T);
    for (int i = 0; i < lenT; i++) {
        double M = 0, D = 0;
        if (T[i] < L_Q[i]) {
            M = edr_subcost(T[i], L_Q[i], epsilon);
            D = 1;
            d1[i + w] = MIN(M, D);
        } else if (T[i] > U_Q[i]) {
            M = edr_subcost(T[i], U_Q[i], epsilon);
            D = 1;
            d1[i + w] = MIN(M, D);
        }
        d_f1 += d1[i + w];
        if (d_boundary + d_f1 > bsf) return MAX_VAL;
    }

    double d_b1 = 0.0;
    double d_f2 = 0.0;

    int pad_Q = lenQ + 2 * w;
    double d2[pad_Q];
    memset(d2, 0, sizeof(double) * pad_Q);
    for (int i = 0; i < lenQ; i++) {
        double M = 0, D = 0;
        if (Q[i] < L_T[i]) {
            M = edr_subcost(Q[i], L_T[i], epsilon);
            D = 1;
            d2[i + w] = MIN(M, D);
        } else if (Q[i] > U_T[i]) {
            M = edr_subcost(Q[i], U_T[i], epsilon);
            D = 1;
            d2[i + w] = MIN(M, D);
        }
        d_f2 += d2[i + w];
        if (d_boundary + d_f2 > bsf) return MAX_VAL;
    }

    double d1_copy[pad_T];
    memcpy(d1_copy, d1, sizeof(double) * pad_T);
    int next = 0;
    for (int i = 0; i < lenQ; i++) {
        double current_d2 = d2[i + w];
        if (current_d2 > d1_copy[i]) {
            next = MAX(next, i);
            for (int j = next; j <= i + 2 * w; j++) {
                if (current_d2 > d1_copy[j]) {
                    current_d2 -= d1_copy[j];
                    d1_copy[j] = 0;
                    next = j + 1;
                } else {
                    d1_copy[j] -= current_d2;
                    current_d2 = 0;
                    next = j;
                    break;
                }
            }
            d_b1 += current_d2;
            if (d_boundary + d_f1 + d_b1 > bsf) return MAX_VAL;
        }
    }

    double d_b2 = 0.0;
    double d2_copy[pad_Q];
    memcpy(d2_copy, d2, sizeof(double) * pad_Q);
    next = 0;
    for (int i = 0; i < lenT; i++) {
        double current_d1 = d1[i + w];
        if (current_d1 > d2_copy[i]) {
            next = MAX(next, i);
            for (int j = next; j <= i + 2 * w; j++) {
                if (current_d1 > d2_copy[j]) {
                    current_d1 -= d2_copy[j];
                    d2_copy[j] = 0;
                    next = j + 1;
                } else {
                    d2_copy[j] -= current_d1;
                    current_d1 = 0;
                    next = j;
                    break;
                }
            }
            d_b2 += current_d1;
            if (d_boundary + d_f2 + d_b2 > bsf) return MAX_VAL;
        }
    }

    return d_boundary + MAX(d_f1 + d_b1, d_f2 + d_b2);
}


double bglb_lcss(double* Q, int lenQ, double* T, int lenT, double epsilon, double* U_Q, double* L_Q, double* U_T, double* L_T, int w, double& bsf) {
    double min_len = MIN(lenQ, lenT);
    double scaled_bsf = bsf * min_len;
    double edr_dist = bglb_edr(Q, lenQ, T, lenT, epsilon, U_Q, L_Q, U_T, L_T, w, scaled_bsf);

    if (edr_dist == MAX_VAL) return MAX_VAL;
    return edr_dist / min_len;
}


double bglb_erp(double* Q, int lenQ, double* T, int lenT, double g_val, double* U_Q, double* L_Q, double* U_T, double* L_T, int w, double& bsf) {
    double bsf_sq = bsf * bsf;
    double d_boundary = min({
        dist(T[lenT-1], Q[lenQ-1]),
        dist(T[lenT-1], g_val),
        dist(Q[lenQ-1], g_val)
    });

    double d_f1 = 0.0;
    double d1[lenT], M1[lenT], D1[lenT];
    memset(d1, 0, sizeof(double) * lenT);
    memset(M1, 0, sizeof(double) * lenT);
    memset(D1, 0, sizeof(double) * lenT);

    for (int i = 1; i < lenT - 1; i++) {
        if (T[i] < L_Q[i]) {
            M1[i] = dist(T[i], L_Q[i]);
            D1[i] = dist(T[i], g_val);
            d1[i] = MIN(M1[i], D1[i]);
        } else if (T[i] > U_Q[i]) {
            M1[i] = dist(T[i], U_Q[i]);
            D1[i] = dist(T[i], g_val);
            d1[i] = MIN(M1[i], D1[i]);
        }
        d_f1 += d1[i];
        if (d_boundary + d_f1 > bsf_sq) return MAX_VAL;
    }

    double d_b1 = 0.0;
    double d_f2 = 0.0;
    double d1_U[lenT];
    memset(d1_U, 0, sizeof(double) * lenT);
    upper_lemire(d1, lenT, w, d1_U);

    double d2[lenQ], M2[lenQ], D2[lenQ];
    memset(d2, 0, sizeof(double) * lenQ);
    memset(M2, 0, sizeof(double) * lenQ);
    memset(D2, 0, sizeof(double) * lenQ);

    for (int i = 1; i < lenQ - 1; i++) {
        if (Q[i] < L_T[i]) {
            M2[i] = dist(Q[i], L_T[i]);
            D2[i] = dist(Q[i], g_val);
            d2[i] = MIN(M2[i], D2[i]);
        } else if (Q[i] > U_T[i]) {
            M2[i] = dist(Q[i], U_T[i]);
            D2[i] = dist(Q[i], g_val);
            d2[i] = MIN(M2[i], D2[i]);
        }
        d_b1 += MIN(MAX(M2[i] - d1_U[i], 0.0), D2[i]);
        d_f2 += d2[i];
        if (d_boundary + d_f1 + d_b1 > bsf_sq || d_boundary + d_f2 > bsf_sq) return MAX_VAL;
    }

    double d_b2 = 0.0;
    double d2_U[lenQ];
    memset(d2_U, 0, sizeof(double) * lenQ);
    upper_lemire(d2, lenQ, w, d2_U);
    for (int i = 1; i < lenT - 1; i++) {
        d_b2 += MIN(MAX(M1[i] - d2_U[i], 0.0), D1[i]);
        if (d_boundary + d_f2 + d_b2 > bsf_sq) return MAX_VAL;
    }

    return sqrt(max(0.0, d_boundary + MAX(d_f1 + d_b1, d_f2 + d_b2)));
}


double bglb_msm(double* Q, int lenQ, double* T, int lenT, double c, double* U_Q, double* L_Q, double* U_T, double* L_T, int w, double& bsf) {
    double d_boundary;
    if ((T[lenT-2] >= T[lenT-1] && T[lenT-1] >= Q[lenQ-1]) ||
        (T[lenT-2] <= T[lenT-1] && T[lenT-1] <= Q[lenQ-1]) ||
        (Q[lenQ-2] <= Q[lenQ-1] && Q[lenQ-1] <= T[lenT-1]) ||
        (Q[lenQ-2] >= Q[lenQ-1] && Q[lenQ-1] >= T[lenT-1])) {
        d_boundary = abs(T[0] - Q[0]) + MIN(abs(T[lenT-1] - Q[lenQ-1]), c);
    } else {
        d_boundary = abs(T[0] - Q[0]) + min({
            abs(T[lenT-1] - Q[lenQ-1]),
            c + abs(Q[lenQ-1] - Q[lenQ-2]),
            c + abs(T[lenT-1] - T[lenT-2])
        });
    }

    double d_f1 = 0.0;
    double d1[lenT], M1[lenT], D1[lenT];
    memset(d1, 0, sizeof(double) * lenT);
    memset(M1, 0, sizeof(double) * lenT);
    memset(D1, 0, sizeof(double) * lenT);

    for (int i = 1; i < lenT - 1; i++) {
        if (T[i] < L_Q[i]) {
            M1[i] = abs(T[i] - L_Q[i]);
            D1[i] = c;
            d1[i] = MIN(M1[i], D1[i]);
        } else if (T[i] > U_Q[i]) {
            M1[i] = abs(T[i] - U_Q[i]);
            D1[i] = c;
            d1[i] = MIN(M1[i], D1[i]);
        }
        d_f1 += d1[i];
        if (d_boundary + d_f1 > bsf) return MAX_VAL;
    }

    double d_b1 = 0.0;
    double d_f2 = 0.0;
    double d1_U[lenT];
    memset(d1_U, 0, sizeof(double) * lenT);
    upper_lemire(d1, lenT, w, d1_U);

    double d2[lenQ], M2[lenQ], D2[lenQ];
    memset(d2, 0, sizeof(double) * lenQ);
    memset(M2, 0, sizeof(double) * lenQ);
    memset(D2, 0, sizeof(double) * lenQ);

    for (int i = 1; i < lenQ - 1; i++) {
        if (Q[i] < L_T[i]) {
            M2[i] = abs(Q[i] - L_T[i]);
            D2[i] = c;
            d2[i] = MIN(M2[i], D2[i]);
        } else if (Q[i] > U_T[i]) {
            M2[i] = abs(Q[i] - U_T[i]);
            D2[i] = c;
            d2[i] = MIN(M2[i], D2[i]);
        }
        d_b1 += MIN(MAX(M2[i] - d1_U[i], 0.0), D2[i]);
        d_f2 += d2[i];
        if (d_boundary + d_f1 + d_b1 > bsf || d_boundary + d_f2 > bsf) return MAX_VAL;
    }

    double d_b2 = 0.0;
    double d2_U[lenQ];
    memset(d2_U, 0, sizeof(double) * lenQ);
    upper_lemire(d2, lenQ, w, d2_U);
    for (int i = 1; i < lenT - 1; i++) {
        d_b2 += MIN(MAX(M1[i] - d2_U[i], 0.0), D1[i]);
        if (d_boundary + d_f2 + d_b2 > bsf) return MAX_VAL;
    }

    return d_boundary + MAX(d_f1 + d_b1, d_f2 + d_b2);
}

MetricType metric_from_tag(const char* tag) {
    if (strcmp("dtw", tag) == 0) return METRIC_DTW;
    if (strcmp("edr", tag) == 0) return METRIC_EDR;
    if (strcmp("erp", tag) == 0) return METRIC_ERP;
    if (strcmp("lcss", tag) == 0) return METRIC_LCSS;
    if (strcmp("msm", tag) == 0) return METRIC_MSM;
    return METRIC_UNKNOWN;
}

static bool bglb_dtw_forward_pruned(double* Q, int lenQ, double* T, int lenT, double* L_Q, double* U_Q, double bsf) {
    double bsf_sq = bsf * bsf;
    double d_boundary = dist(T[0], Q[0]) + dist(T[lenT - 1], Q[lenQ - 1]);
    double d_f1 = 0.0;
    for (int i = 1; i < lenT - 1; i++) {
        if (T[i] < L_Q[i]) {
            d_f1 += dist(T[i], L_Q[i]);
        } else if (T[i] > U_Q[i]) {
            d_f1 += dist(T[i], U_Q[i]);
        }
        if (d_boundary + d_f1 > bsf_sq) return true;
    }
    return false;
}

static bool bglb_edr_forward_pruned(double* T, int lenT, double epsilon, double* L_Q, double* U_Q, double bsf) {
    double d_f1 = 0.0;
    for (int i = 0; i < lenT; i++) {
        if (T[i] < L_Q[i]) {
            d_f1 += MIN(edr_subcost(T[i], L_Q[i], epsilon), 1.0);
        } else if (T[i] > U_Q[i]) {
            d_f1 += MIN(edr_subcost(T[i], U_Q[i], epsilon), 1.0);
        }
        if (d_f1 > bsf) return true;
    }
    return false;
}

static bool bglb_erp_forward_pruned(double* Q, int lenQ, double* T, int lenT, double g_val, double* L_Q, double* U_Q, double bsf) {
    double bsf_sq = bsf * bsf;
    double d_boundary = min({
        dist(T[lenT - 1], Q[lenQ - 1]),
        dist(T[lenT - 1], g_val),
        dist(Q[lenQ - 1], g_val)
    });
    double d_f1 = 0.0;
    for (int i = 1; i < lenT - 1; i++) {
        if (T[i] < L_Q[i]) {
            d_f1 += MIN(dist(T[i], L_Q[i]), dist(T[i], g_val));
        } else if (T[i] > U_Q[i]) {
            d_f1 += MIN(dist(T[i], U_Q[i]), dist(T[i], g_val));
        }
        if (d_boundary + d_f1 > bsf_sq) return true;
    }
    return false;
}

static bool bglb_msm_forward_pruned(double* Q, int lenQ, double* T, int lenT, double c, double* L_Q, double* U_Q, double bsf) {
    double d_boundary;
    if ((T[lenT - 2] >= T[lenT - 1] && T[lenT - 1] >= Q[lenQ - 1]) ||
        (T[lenT - 2] <= T[lenT - 1] && T[lenT - 1] <= Q[lenQ - 1]) ||
        (Q[lenQ - 2] <= Q[lenQ - 1] && Q[lenQ - 1] <= T[lenT - 1]) ||
        (Q[lenQ - 2] >= Q[lenQ - 1] && Q[lenQ - 1] >= T[lenT - 1])) {
        d_boundary = abs(T[0] - Q[0]) + MIN(abs(T[lenT - 1] - Q[lenQ - 1]), c);
    } else {
        d_boundary = abs(T[0] - Q[0]) + min({
            abs(T[lenT - 1] - Q[lenQ - 1]),
            c + abs(Q[lenQ - 1] - Q[lenQ - 2]),
            c + abs(T[lenT - 1] - T[lenT - 2])
        });
    }

    double d_f1 = 0.0;
    for (int i = 1; i < lenT - 1; i++) {
        if (T[i] < L_Q[i]) {
            d_f1 += MIN(abs(T[i] - L_Q[i]), c);
        } else if (T[i] > U_Q[i]) {
            d_f1 += MIN(abs(T[i] - U_Q[i]), c);
        }
        if (d_boundary + d_f1 > bsf) return true;
    }
    return false;
}

double bglb_elastic_lazy(MetricType metric, double* Q, int lenQ, double* T, int lenT, double arg, double* L_Q, double* U_Q, int w, double& bsf) {
    switch (metric) {
        case METRIC_DTW:
            if (bglb_dtw_forward_pruned(Q, lenQ, T, lenT, L_Q, U_Q, bsf)) return MAX_VAL;
            break;
        case METRIC_EDR:
            if (bglb_edr_forward_pruned(T, lenT, arg, L_Q, U_Q, bsf)) return MAX_VAL;
            break;
        case METRIC_LCSS: {
            double scaled_bsf = bsf * MIN(lenQ, lenT);
            if (bglb_edr_forward_pruned(T, lenT, arg, L_Q, U_Q, scaled_bsf)) return MAX_VAL;
            break;
        }
        case METRIC_ERP:
            if (bglb_erp_forward_pruned(Q, lenQ, T, lenT, arg, L_Q, U_Q, bsf)) return MAX_VAL;
            break;
        case METRIC_MSM:
            if (bglb_msm_forward_pruned(Q, lenQ, T, lenT, arg, L_Q, U_Q, bsf)) return MAX_VAL;
            break;
        default:
            return MAX_VAL;
    }

    double L_T[lenT];
    double U_T[lenT];
    lower_upper_lemire(T, lenT, w, L_T, U_T);

    switch (metric) {
        case METRIC_DTW:
            return bglb_dtw(Q, lenQ, T, lenT, L_Q, U_Q, L_T, U_T, w, bsf);
        case METRIC_EDR:
            return bglb_edr(Q, lenQ, T, lenT, arg, U_Q, L_Q, U_T, L_T, w, bsf);
        case METRIC_LCSS:
            return bglb_lcss(Q, lenQ, T, lenT, arg, U_Q, L_Q, U_T, L_T, w, bsf);
        case METRIC_ERP:
            return bglb_erp(Q, lenQ, T, lenT, arg, U_Q, L_Q, U_T, L_T, w, bsf);
        case METRIC_MSM:
            return bglb_msm(Q, lenQ, T, lenT, arg, U_Q, L_Q, U_T, L_T, w, bsf);
        default:
            return MAX_VAL;
    }
}


