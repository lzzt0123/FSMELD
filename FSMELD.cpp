


#include <stdlib.h>
#include <string.h>
#include "fstream"
#include<iostream>
#include "def.h"
#include "funcs.h"
#include "fftw3.h"
#include <cmath>
#include <algorithm>
#include <cstring>   
#include <cfloat> 
#include <iomanip> 
#include <chrono>

using namespace std;
static constexpr double DIST_TOL = 1e-9;

struct CandidateExclusionState {
    bool active = false;
    long long begin = -1;
    long long end = -1;
    long long next_window_start = 0;
    long long excluded_windows = 0;
};

static thread_local CandidateExclusionState candidate_exclusion;

// Temporarily excludes target windows that overlap a known query segment.
class ScopedCandidateExclusion {
public:
    ScopedCandidateExclusion(long long begin, long long end)
        : previous_(candidate_exclusion) {
        candidate_exclusion = CandidateExclusionState{};
        if (begin >= 0 && end > begin) {
            candidate_exclusion.active = true;
            candidate_exclusion.begin = begin;
            candidate_exclusion.end = end;
        }
    }

    ~ScopedCandidateExclusion() {
        candidate_exclusion = previous_;
    }

    long long visited_windows() const {
        return candidate_exclusion.next_window_start;
    }

    long long excluded_windows() const {
        return candidate_exclusion.excluded_windows;
    }

private:
    CandidateExclusionState previous_;
};

double max_mean = -DBL_MAX, min_mean = DBL_MAX;
double max_std  = -DBL_MAX, min_std  = DBL_MAX;
ofstream csv_file;
#include <queue>   
#include <vector>
#include <utility> 

void trans(double z_l[],double z_u[],double z[],int m,int k_m_l,char *tag) {
    if(strcmp(tag,"erp")==0) {
        for(int i=0,j=m-1;i<k_m_l;i++,j++) {
            z[i]=std::sqrt(std::max(0.0, z_l[j]+z_u[j]));
        }
    }else if (strcmp(tag,"dtw")==0) {
        for(int i=0,j=m-1;i<k_m_l;i++,j++) {
            z[i]=std::sqrt(std::max(0.0, z_l[j]+z_u[j]));
        }
    }else if (strcmp(tag,"msm")==0) {
        for(int i=0,j=m-1;i<k_m_l;i++,j++) {
            z[i]=z_l[j]+z_u[j];
        }
    }else if (strcmp(tag,"edr")==0) {
        for(int i=0,j=m-1;i<k_m_l;i++,j++) {
            z[i]=z_l[j]+z_u[j];
        }
    }else if (strcmp(tag,"lcss")==0) {
        for(int i=0,j=m-1;i<k_m_l;i++,j++) {
            z[i] = std::max(0.0, (z_l[j] + z_u[j]) / m);
        }
    }
}



static void precompute_invalid_mask_fast(const double* t, const double* si_t, int m_k, int m, bool* invalid_mask) {
    const int window_count = m_k - m + 1;
    if (window_count <= 0) return;

    int invalid_count = 0;
    for (int k = 0; k < m; k++) {
        if (!std::isfinite(t[k])) invalid_count++;
    }
    invalid_mask[0] = invalid_count > 0;

    for (int j = 1; j < window_count; j++) {
        if (!std::isfinite(t[j - 1])) invalid_count--;
        if (!std::isfinite(t[j + m - 1])) invalid_count++;
        invalid_mask[j] = invalid_count > 0 || !std::isfinite(si_t[j]) || si_t[j] <= 1e-12;
    }

    invalid_mask[0] = invalid_mask[0] || !std::isfinite(si_t[0]) || si_t[0] <= 1e-12;

    if (candidate_exclusion.active) {
        const long long block_start = candidate_exclusion.next_window_start;
        const long long overlap_begin = std::max(
            block_start, candidate_exclusion.begin);
        const long long overlap_end = std::min(
            block_start + window_count, candidate_exclusion.end);
        for (long long global_index = overlap_begin; global_index < overlap_end; ++global_index) {
            invalid_mask[global_index - block_start] = true;
        }
        if (overlap_end > overlap_begin) {
            candidate_exclusion.excluded_windows += overlap_end - overlap_begin;
        }
        candidate_exclusion.next_window_start += window_count;
    }
}

static double exact_elastic_distance(MetricType metric, double* q, double* t,
                                     int m, int w, double arg) {
    switch (metric) {
        case METRIC_MSM: return msm(q, m, t, m, arg, w);
        case METRIC_ERP: return erp(q, m, t, m, arg, w);
        case METRIC_LCSS: return lcss(q, m, t, m, w, arg);
        case METRIC_EDR: return edr(q, m, t, m, arg, w);
        case METRIC_DTW: return dtw(q, m, t, m, w);
        default: return NAN;
    }
}


static void find_threshold_impl(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w, double arg, int top_k, bool filter_by_initial_epsilon) {
    
    priority_queue<pair<double, long long>> top_k_heap;
    vector<long long> A; 
    MetricType metric = metric_from_tag(tag);

    int i;
    int m_k;
    double d;
    int p;

    if (strcmp(tag, "erp") == 0) p = 2;
    else if (strcmp(tag, "dtw") == 0) p = 2;
    else if (strcmp(tag, "msm") == 0) p = 1;
    else if (strcmp(tag, "edr") == 0) p = 0;
    else if (strcmp(tag, "lcss") == 0) p = 0;

    
    double ex_q = 0, ex2_q = 0, std_val = 0, mean = 0;
    unsigned long long kk1 = 0; 
    unsigned long long exact_count = 0;
    unsigned long long lb_pruned = 0;
    unsigned long long bglb_pruned = 0;

    
    m_k = 1024;
    while (m_k <= 4 * m) {
        m_k = m_k * 2;
    }
    double* t = (double*)malloc(sizeof(double) * m_k);
    double* q = (double*)malloc(sizeof(double) * m);
    double* q_normal = (double*)malloc(sizeof(double) * m);

    
    i = 0;
    FILE* qp = fopen(bi_q, "r");
    while (fscanf(qp, "%lf", &d) != EOF && i < m) {
        ex_q += d;
        ex2_q += d * d;
        q[i] = d;
        i++;
    }
    fclose(qp);

    mean = ex_q / m;
    std_val = ex2_q / m;
    std_val = sqrt(std_val - mean * mean);
    fflush(stdout);

    for (i = 0; i < m; i++)
        q_normal[i] = (q[i] - mean) / std_val;

    double* q_sorted = (double*)malloc(sizeof(double) * m);
    memcpy(q_sorted, q_normal, sizeof(double) * m);
    std::sort(q_sorted, q_sorted + m);

    
    streampos sizeofT;
    string t_name = bi_t;
    ifstream input(t_name.c_str(), ios::binary | ios::in);
    input.seekg(0, ios::end);
    sizeofT = input.tellg();
    input.seekg(0, ios::beg);

    long long k_m_l = m_k - m + 1;
    long long ts_n = (sizeofT / (sizeof(double))) - m + 1;

    std::cout << "ts_n=" << ts_n << std::endl;
    std::cout << "m=" << m << std::endl;
    std::cout << "k_m_l=" << k_m_l << std::endl;
    std::cout << "m_k=" << m_k << std::endl;

    long long forsize = (ts_n) / (k_m_l);
    int flag_of_tail = 0;
    double t_u, t_l;

    size_t tail = ts_n % (k_m_l);
    int siz_t = tail + m - 1;

    if (tail == 0) flag_of_tail = 1;
    if (siz_t <= 2 * m - 2) {
        tail = m;
        siz_t = 2 * m - 1;
    }

    int forsize_divide_100 = forsize / 100;
    int i_ = 0;

    
    double* Z = (double*)malloc(sizeof(double) * m_k);
    fftw_complex* XX = fftw_alloc_complex(m_k / 2 + 1);
    fftw_complex* YY = fftw_alloc_complex(m_k / 2 + 1);
    fftw_complex* ZZ = fftw_alloc_complex(m_k / 2 + 1);
    fftw_plan p_r2c = fftw_plan_dft_r2c_1d(m_k, Z, ZZ, FFTW_ESTIMATE);
    fftw_plan p_c2r = fftw_plan_dft_c2r_1d(m_k, ZZ, Z, FFTW_ESTIMATE);

    double* Z_tail = (double*)malloc(sizeof(double) * siz_t);
    fftw_complex* XX_tail = fftw_alloc_complex(siz_t / 2 + 1);
    fftw_complex* YY_tail = fftw_alloc_complex(siz_t / 2 + 1);
    fftw_complex* ZZ_tail = fftw_alloc_complex(siz_t / 2 + 1);
    fftw_plan p_r2c_tail = fftw_plan_dft_r2c_1d(siz_t, Z_tail, ZZ_tail, FFTW_ESTIMATE);
    fftw_plan p_c2r_tail = fftw_plan_dft_c2r_1d(siz_t, ZZ_tail, Z_tail, FFTW_ESTIMATE);

    double* L_sort = (double*)malloc(sizeof(double) * (m_k - 2 * m + 2));
    double* U_sort = (double*)malloc(sizeof(double) * (m_k - 2 * m + 2));
    double* LB = (double*)malloc(sizeof(double) * (m_k - 2 * m + 2));
    double* m_p = (double*)malloc(sizeof(double) * (p + 1));

    double* h_l = (double*)malloc(sizeof(double) * m);
    double* h_u = (double*)malloc(sizeof(double) * m);
    double* U = (double*)malloc(sizeof(double) * m_k);
    double* L = (double*)malloc(sizeof(double) * m_k);
    double* U_normal = (double*)malloc(sizeof(double) * m_k);
    double* L_normal = (double*)malloc(sizeof(double) * m_k);
    double* g_l = (double*)malloc(sizeof(double) * m_k);
    double* g_u = (double*)malloc(sizeof(double) * m_k);
    double* z_l = (double*)malloc(sizeof(double) * m_k);
    double* z_u = (double*)malloc(sizeof(double) * m_k);
    double* z = (double*)malloc(sizeof(double) * k_m_l);
    double* miu_t = (double*)malloc(sizeof(double) * k_m_l);
    double* si_t = (double*)malloc(sizeof(double) * k_m_l);
    double* miu_min = (double*)malloc(sizeof(double) * m_k);
    double* miu_max = (double*)malloc(sizeof(double) * m_k);
    double* si_min = (double*)malloc(sizeof(double) * m_k);
    double* si_max = (double*)malloc(sizeof(double) * m_k);
    bool* invalid_mask = (bool*)malloc(sizeof(bool) * m_k);
    double* t_temp = (double*)malloc(sizeof(double) * m);
    double *U_Q = (double*) malloc(sizeof(double) * m);
    double *L_Q = (double*) malloc(sizeof(double) * m);
    lower_upper_lemire(q_normal, m, w, L_Q, U_Q);

    for (long long ll = 0; ll < forsize; ll++) {
        if (forsize_divide_100) {
            if ((ll % (forsize_divide_100)) == 0) {
                std::cout << "\b\b\b\b";
                std::cout << "|" << setw(3) << i_ << "%" << std::endl;
                i_++;


                std::cout << "Threshold=" << std::setprecision(17) << threshold << std::endl;
            }
        }

        if (ll != 0) memcpy(t, t + k_m_l, sizeof(double) * m);
        else input.read((char*)t, (m - 1) * sizeof(double));
        input.read((char*)(t + m - 1), (k_m_l) * sizeof(double));


        memset(miu_min, 0, sizeof(double) * m_k);
        memset(miu_max, 0, sizeof(double) * m_k);
        memset(si_min, 0, sizeof(double) * m_k);
        memset(si_max, 0, sizeof(double) * m_k);

        mvmean(t, m_k, m, miu_t, si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, m_k, m, invalid_mask);
        min_max_window(miu_t, m_k, m, miu_min, miu_max);
        min_max_window(si_t, m_k, m, si_min, si_max);

        memset(h_u, 0, sizeof(double) * m); memset(U, 0, sizeof(double) * m_k); memset(g_u, 0, sizeof(double) * m_k); memset(U_normal, 0, sizeof(double) * m_k);
        memset(h_l, 0, sizeof(double) * m); memset(L, 0, sizeof(double) * m_k); memset(g_l, 0, sizeof(double) * m_k); memset(L_normal, 0, sizeof(double) * m_k);
        memset(z_l, 0, sizeof(double) * m_k); memset(z_u, 0, sizeof(double) * m_k); memset(z, 0, sizeof(double) * k_m_l);

        lower_upper_lemire(t, m_k, w, L, U);
        compute_normalized_arrays(U, L, U_normal, L_normal, m_k, miu_max, miu_min, si_max, si_min);

        if (strcmp(tag, "dtw") == 0) {
            t_u = find_tu_lp_dtw(U_normal, m_k, q_normal, m, p, LB, m_p, tag, q_sorted);
            t_l = find_tl_lp_dtw(L_normal, m_k, q_normal, m, p, LB, m_p, tag, q_sorted);
        } else if(strcmp(tag,"erp")==0) {
            t_u=find_tu_lp_erp(U_normal,m_k,q_normal,m,p,arg,LB,m_p);
            t_l=find_tl_lp_erp(L_normal,m_k,q_normal,m,p,arg,LB,m_p);
        }else {
            t_u = find_tu_lp(U_normal, m_k, q_normal, m, p, arg, LB, m_p, tag, q_sorted);
            t_l = find_tl_lp(L_normal, m_k, q_normal, m, p, arg, LB, m_p, tag, q_sorted);
        }

        q_to_h_l(q_normal, m, t_l, arg, h_l, m, tag); q_to_h_u(q_normal, m, t_u, arg, h_u, m, tag);
        s_to_g_l(L_normal, m_k, t_l, g_l, m_k, tag, arg); s_to_g_u(U_normal, m_k, t_u, g_u, m_k, tag, arg);

        inv_and_padding_y(m_k, m, h_l, YY, p_r2c); FFT_x(g_l, XX, p_r2c); cov_IFFT(XX, YY, ZZ, m_k, z_l, p_c2r);
        inv_and_padding_y(m_k, m, h_u, YY, p_r2c); FFT_x(g_u, XX, p_r2c); cov_IFFT(XX, YY, ZZ, m_k, z_u, p_c2r);

        trans(z_l,z_u,z,m,k_m_l,tag);
        for (int j = 0, i = m - 1; j < k_m_l; j++, i++) {
            
            if (invalid_mask[j]) {
                continue;
            }

            kk1++;
            bool heap_full = top_k_heap.size() >= static_cast<size_t>(top_k);
            double active_threshold = (filter_by_initial_epsilon || heap_full) ? threshold : DBL_MAX;
            if (z[j] > active_threshold + DIST_TOL) {
                lb_pruned++;
                continue;
            }

            for (int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            if (filter_by_initial_epsilon || heap_full) {
                double bsf = threshold;
                double bglb_dist = bglb_elastic_lazy(metric, q_normal, m, t_temp, m, arg, L_Q, U_Q, w, bsf);
                if (bglb_dist > threshold + DIST_TOL) {
                    bglb_pruned++;
                    continue;
                }
            }

            double actual_dist = 0.0;
            if (strcmp("msm", tag) == 0) actual_dist = msm(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("dtw", tag) == 0) actual_dist = dtw(q_normal, m, t_temp, m, w);
            else if (strcmp("erp", tag) == 0) actual_dist = erp(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("edr", tag) == 0) actual_dist = edr(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("lcss", tag) == 0) actual_dist = lcss(q_normal, m, t_temp, m, w, arg);

            exact_count++;
            if (!std::isfinite(actual_dist)) {
                continue;
            }
            if (filter_by_initial_epsilon && actual_dist > threshold + DIST_TOL) {
                continue;
            }
            long long current_index = ll * k_m_l + j;
            if (top_k_heap.size() < static_cast<size_t>(top_k)) {
                top_k_heap.push({actual_dist, current_index});
                if (top_k_heap.size() == static_cast<size_t>(top_k)) {
                    threshold = top_k_heap.top().first;
                }
            } else if (actual_dist < top_k_heap.top().first) {
                top_k_heap.pop();
                top_k_heap.push({actual_dist, current_index});
                threshold = top_k_heap.top().first;
            }
        }
    }

    
    if (!flag_of_tail) {
        if (siz_t > 2 * m - 2) {
            for (int f = 0; f < m - 1; f++) t[f] = t[k_m_l + f];
            input.read((char*)(t + m - 1), tail * sizeof(double));
        } else {
            input.read((char*)t, siz_t * sizeof(double));
        }


        mvmean(t, siz_t, m, miu_t, si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, siz_t, m, invalid_mask);
        min_max_window(miu_t, siz_t, m, miu_min, miu_max);
        min_max_window(si_t, siz_t, m, si_min, si_max);

        lower_upper_lemire(t, siz_t, w, L, U);
        compute_normalized_arrays(U, L, U_normal, L_normal, siz_t, miu_max, miu_min, si_max, si_min);

        if (strcmp(tag, "dtw") == 0) {
            t_u = find_tu_lp_dtw(U_normal, siz_t, q_normal, m, p, LB, m_p, tag, q_sorted);
            t_l = find_tl_lp_dtw(L_normal, siz_t, q_normal, m, p, LB, m_p, tag, q_sorted);
        } else if(strcmp(tag,"erp")==0) {
            t_u=find_tu_lp_erp(U_normal,siz_t,q_normal,m,p,arg,LB,m_p);
            t_l=find_tl_lp_erp(L_normal,siz_t,q_normal,m,p,arg,LB,m_p);
        }else {
            t_u = find_tu_lp(U_normal, siz_t, q_normal, m, p, arg, LB, m_p, tag, q_sorted);
            t_l = find_tl_lp(L_normal, siz_t, q_normal, m, p, arg, LB, m_p, tag, q_sorted);
        }

        memset(h_u, 0, sizeof(double) * m); memset(U, 0, sizeof(double) * siz_t); memset(g_u, 0, sizeof(double) * siz_t);
        memset(h_l, 0, sizeof(double) * m); memset(L, 0, sizeof(double) * siz_t); memset(g_l, 0, sizeof(double) * siz_t);
        memset(z_l, 0, sizeof(double) * siz_t); memset(z_u, 0, sizeof(double) * siz_t); memset(z, 0, sizeof(double) * tail);

        q_to_h_l(q_normal, m, t_l, arg, h_l, m, tag); q_to_h_u(q_normal, m, t_u, arg, h_u, m, tag);
        s_to_g_l(L_normal, siz_t, t_l, g_l, siz_t, tag, arg); s_to_g_u(U_normal, siz_t, t_u, g_u, siz_t, tag, arg);

        inv_and_padding_y(siz_t, m, h_l, YY_tail, p_r2c_tail); FFT_x(g_l, XX_tail, p_r2c_tail); cov_IFFT(XX_tail, YY_tail, ZZ_tail, siz_t, z_l, p_c2r_tail);
        inv_and_padding_y(siz_t, m, h_u, YY_tail, p_r2c_tail); FFT_x(g_u, XX_tail, p_r2c_tail); cov_IFFT(XX_tail, YY_tail, ZZ_tail, siz_t, z_u, p_c2r_tail);

        for (int j = 0, i = m - 1; j < tail; j++, i++) {
            if (invalid_mask[j]) {
                continue;
            }
            if (strcmp("lcss", tag) == 0) z[j] = std::max(0.0, (z_l[i] + z_u[i]) / m);
            else if (strcmp("dtw", tag) == 0 || strcmp("erp", tag) == 0) z[j] = std::sqrt(std::max(0.0, z_l[i] + z_u[i]));
            else z[j] = z_l[i] + z_u[i];

            kk1++;
            bool heap_full = top_k_heap.size() >= static_cast<size_t>(top_k);
            double active_threshold = (filter_by_initial_epsilon || heap_full) ? threshold : DBL_MAX;
            if (z[j] > active_threshold + DIST_TOL) {
                lb_pruned++;
                continue;
            }

            for (int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            if (filter_by_initial_epsilon || heap_full) {
                double bsf = threshold;
                double bglb_dist = bglb_elastic_lazy(metric, q_normal, m, t_temp, m, arg, L_Q, U_Q, w, bsf);
                if (bglb_dist > threshold + DIST_TOL) {
                    bglb_pruned++;
                    continue;
                }
            }

            double actual_dist = 0.0;
            if (strcmp("msm", tag) == 0) actual_dist = msm(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("dtw", tag) == 0) actual_dist = dtw(q_normal, m, t_temp, m, w);
            else if (strcmp("erp", tag) == 0) actual_dist = erp(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("edr", tag) == 0) actual_dist = edr(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("lcss", tag) == 0) actual_dist = lcss(q_normal, m, t_temp, m, w, arg);

            exact_count++;
            if (!std::isfinite(actual_dist)) {
                continue;
            }
            if (filter_by_initial_epsilon && actual_dist > threshold + DIST_TOL) {
                continue;
            }
            long long current_index = forsize * k_m_l + j;
            if (top_k_heap.size() < static_cast<size_t>(top_k)) {
                top_k_heap.push({actual_dist, current_index});
                if (top_k_heap.size() == static_cast<size_t>(top_k)) {
                    threshold = top_k_heap.top().first;
                }
            } else if (actual_dist < top_k_heap.top().first) {
                top_k_heap.pop();
                top_k_heap.push({actual_dist, current_index});
                threshold = top_k_heap.top().first;
            }
        }
    }

    std::cout << "\n\n" << std::string(60, '=') << std::endl;
    std::cout << " [FSMELD Top-K Threshold Results]" << std::endl;
    std::cout << "  - Candidates checked by LB:   " << kk1 << std::endl;
    std::cout << "  - Pruned by LBIP:             " << lb_pruned << std::endl;
    std::cout << "  - Pruned by BGLB:             " << bglb_pruned << std::endl;
    std::cout << "  - Exact distances computed:   " << exact_count << std::endl;
    std::cout << std::setprecision(17);
    std::cout << "  - Final Top-K Threshold:      " << threshold << std::endl;
    if (filter_by_initial_epsilon) {
        std::cout << "  - Final FSMELD epsilon:       " << threshold << std::endl;
    }
    std::cout << std::string(60, '=') << "\n" << std::endl;

    
    while(!top_k_heap.empty()) {
        A.push_back(top_k_heap.top().second);
        top_k_heap.pop();
    }
    std::reverse(A.begin(), A.end()); 
    for (size_t k = 0; k < A.size(); ++k) {
        std::cout << "Rank " << k+1 << ": Index A[" << k+1 << "]=" << A[k] << std::endl;
    }

    
    free(t_temp);
    free(invalid_mask);
    free(z); free(z_l); free(z_u); free(U); free(L); free(h_l); free(h_u); free(g_l); free(g_u);
    free(Z); free(Z_tail); free(m_p); free(LB); free(U_sort); free(L_sort);
    free(t); free(q); free(q_normal); free(q_sorted);
    free(miu_t); free(si_t); free(miu_min); free(miu_max); free(si_min); free(si_max);
    free(U_normal); free(L_normal); free(U_Q); free(L_Q);

    fftw_free(XX); fftw_free(YY); fftw_free(ZZ);
    fftw_free(XX_tail); fftw_free(YY_tail); fftw_free(ZZ_tail);

    fftw_destroy_plan(p_c2r); fftw_destroy_plan(p_c2r_tail);
    fftw_destroy_plan(p_r2c); fftw_destroy_plan(p_r2c_tail);
}

void find_threshold(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w, double arg, int top_k) {
    find_threshold_impl(bi_t, bi_q, m, threshold, tag, w, arg, top_k, false);
}

void FSMELD_epsilon(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w, double arg, int top_k) {
    find_threshold_impl(bi_t, bi_q, m, threshold, tag, w, arg, top_k, true);
}


void get_tlb(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w, double arg, double target_index) {
    vector<int> A;
    
    int i;
    int m_k;
    double d;
    int p;
    if(strcmp(tag,"erp")==0) {
        p=2;
    }else if (strcmp(tag,"dtw")==0) {
        p=2;
    }else if (strcmp(tag,"msm")==0) {
        p=1;
    }else if (strcmp(tag,"edr")==0) {
        p=0;
    }else if (strcmp(tag,"lcss")==0) {
        p=0;
    }
    double average_ratio = 0.0;
    double max_ratio = 0.0;
    unsigned long long ratio_count = 0;
    unsigned long long zero_distance_count = 0;
    unsigned long long nonfinite_distance_count = 0;
    unsigned long long bound_violation_count = 0;

    double ex_q=0,ex2_q=0,std=0,mean=0;

    unsigned long long kk1=0,kk2=0;

    
    m_k = 1024;
    while(m_k <= 4*m) {
        m_k = m_k*2;  
    }
    double *t = (double *) malloc(sizeof(double) * m_k);
    double *q = (double *) malloc(sizeof(double)*m);
    double *q_normal = (double *) malloc(sizeof(double)*m);
    
    i=0;
    FILE *qp;
    qp = fopen(bi_q,"r");  
    while(fscanf(qp,"%lf",&d) != EOF && i < m)  
    {
        ex_q += d;
        ex2_q += d*d;
        q[i] = d;  
        i++;       
    }
    fclose(qp);    

    mean=ex_q/m;
    std=ex2_q/m;
    std=sqrt(std-mean*mean);
    fflush(stdout);

    
    for( i = 0 ; i < m ; i++ )
        q_normal[i] = (q[i] - mean)/std;

    double* q_sorted = (double*)malloc(sizeof(double) * m);
    memcpy(q_sorted, q_normal, sizeof(double) * m);
    std::sort(q_sorted, q_sorted + m);

    streampos sizeofT;
    string t_name = bi_t;
    ifstream input(t_name.c_str(),ios::binary | ios::in);
    input.seekg(0, ios::end);
    sizeofT = input.tellg();
    input.seekg(0, ios::beg);

    
    long long k_m_l = m_k - m + 1;  
    long long ts_n = (sizeofT/(sizeof(double))) - m + 1;  

    std::cout << "ts_n=" << ts_n << std::endl;
    std::cout << "m=" << m << std::endl;
    std::cout << "k_m_l=" << k_m_l << std::endl;
    std::cout << "m_k=" << m_k << std::endl;
    long long forsize = (ts_n)/(k_m_l);  
    int flag_of_tail = 0;
    double t_u,t_l;

    size_t tail = ts_n % k_m_l;
    int siz_t = tail + m - 1;  

    printf("tail=%zu\n",tail);
    if(tail==0)
    {
        flag_of_tail = 1;
    }
    if(siz_t<=2*m-2) {
        tail=m;
        siz_t=2*m-1;
    }
    printf("forsize=%lld\n",forsize);
    int forsize_divide_100 = forsize/100;
    int i_ = 0;
    
    double *Z = (double*) malloc(sizeof(double)*m_k);

    fftw_complex *XX = fftw_alloc_complex(m_k/2 + 1);
    fftw_complex *YY = fftw_alloc_complex(m_k/2 + 1);
    fftw_complex *ZZ = fftw_alloc_complex(m_k/2 + 1);

    fftw_plan p_r2c = fftw_plan_dft_r2c_1d(m_k, Z, ZZ, FFTW_ESTIMATE);
    fftw_plan p_c2r = fftw_plan_dft_c2r_1d(m_k, ZZ, Z, FFTW_ESTIMATE);

    double *Z_tail = (double*) malloc(sizeof(double)*siz_t);

    fftw_complex *XX_tail = fftw_alloc_complex(siz_t/2 + 1);
    fftw_complex *YY_tail = fftw_alloc_complex(siz_t/2 + 1);
    fftw_complex *ZZ_tail = fftw_alloc_complex(siz_t/2 + 1);

    fftw_plan p_r2c_tail = fftw_plan_dft_r2c_1d(siz_t, Z_tail, ZZ_tail, FFTW_ESTIMATE);
    fftw_plan p_c2r_tail = fftw_plan_dft_c2r_1d(siz_t, ZZ_tail, Z_tail, FFTW_ESTIMATE);

    
    double* L_sort=(double*)malloc(sizeof(double)*(m_k-2*m+2));
    double* U_sort=(double*)malloc(sizeof(double)*(m_k-2*m+2));
    double* LB=(double*)malloc(sizeof(double)*(m_k-2*m+2));
    double* m_p=(double*)malloc(sizeof(double)*(p+1));

    double *h_l = (double*) malloc(sizeof(double) * m);
    double *h_u = (double*) malloc(sizeof(double) * m);

    double *U = (double*) malloc(sizeof(double) * m_k);
    double *L = (double*) malloc(sizeof(double) * m_k);
    double *U_normal = (double*) malloc(sizeof(double) * m_k);
    double *L_normal = (double*) malloc(sizeof(double) * m_k);

    double *g_l = (double*) malloc(sizeof(double) * m_k);
    double *g_u = (double*) malloc(sizeof(double)* m_k);

    double *z_l = (double*) malloc(sizeof(double)*m_k);
    double *z_u = (double*) malloc(sizeof(double)*m_k);
    double *z = (double*) malloc(sizeof(double)*k_m_l);

    double *miu_t = (double*) malloc(sizeof(double)*k_m_l);
    double *si_t = (double*) malloc(sizeof(double)*k_m_l);

    double *miu_min = (double*) malloc(sizeof(double)*m_k);
    double *miu_max = (double*) malloc(sizeof(double)*m_k);
    double *si_min = (double*) malloc(sizeof(double)*m_k);
    double *si_max = (double*) malloc(sizeof(double)*m_k);
    bool* invalid_mask = (bool*)malloc(sizeof(bool) * m_k);

    
    for(long long ll = 0; ll < forsize; ll++) {
        if(forsize_divide_100) {
            if((ll%(forsize_divide_100)) == 0) {
                std::cout << "\b\b\b\b";
                std::cout << "|" << setw(3) << i_ << "%";
                i_++;
                std::cout << "guolv=" << (kk1 > 0 ? double(kk1-kk2)/kk1 : 0.0) << std::endl;
            }
        }
        
        if(ll != 0) {
            memcpy(t, t+k_m_l, sizeof(double)*m);
        } else {
            input.read((char *)t, (m-1)*sizeof(double));
        }
        
        input.read((char *)(t + m -1), (k_m_l)*sizeof(double));

        mvmean(t,m_k,m,miu_t,si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, m_k, m, invalid_mask);

        min_max_window(miu_t,m_k,m,miu_min,miu_max);
        min_max_window(si_t,m_k,m,si_min,si_max);

        memset(h_u, 0, sizeof(double) * m);
        memset(U, 0, sizeof(double) * m_k);
        memset(g_u, 0, sizeof(double) * m_k);
        memset(U_normal, 0, sizeof(double) * m_k);

        memset(h_l, 0, sizeof(double) * m);
        memset(L, 0, sizeof(double) * m_k);
        memset(g_l, 0, sizeof(double) * m_k);
        memset(L_normal, 0, sizeof(double) * m_k);

        memset(z_l, 0, sizeof(double) * m_k);
        memset(z_u, 0, sizeof(double) * m_k);
        memset(z, 0, sizeof(double) * k_m_l);

        lower_upper_lemire(t,m_k,w,L,U);

        compute_normalized_arrays(U,  L, U_normal,  L_normal, m_k,miu_max,miu_min, si_max,si_min);

        if(strcmp(tag,"dtw")==0) {
            t_u=find_tu_lp_dtw(U_normal,m_k,q_normal,m,p,LB,m_p,tag,q_sorted);
            t_l=find_tl_lp_dtw(L_normal,m_k,q_normal,m,p,LB,m_p,tag,q_sorted);
        }else if(strcmp(tag,"erp")==0) {
            t_u=find_tu_lp_erp(U_normal,m_k,q_normal,m,p,arg,LB,m_p);
            t_l=find_tl_lp_erp(L_normal,m_k,q_normal,m,p,arg,LB,m_p);
        }else {
            t_u=find_tu_lp(U_normal,m_k,q_normal,m,p,arg,LB,m_p,tag,q_sorted);
            t_l=find_tl_lp(L_normal,m_k,q_normal,m,p,arg,LB,m_p,tag,q_sorted);
        }

        
        q_to_h_l(q_normal,m,t_l,arg,h_l,m,tag);
        q_to_h_u(q_normal,m,t_u,arg,h_u,m,tag);

        s_to_g_l(L_normal,m_k,t_l, g_l, m_k, tag, arg);
        s_to_g_u(U_normal,m_k,t_u, g_u, m_k, tag, arg);

        
        
        
        inv_and_padding_y(m_k, m, h_l, YY, p_r2c);
        
        FFT_x(g_l, XX, p_r2c);
        
        cov_IFFT(XX, YY, ZZ, m_k, z_l, p_c2r);

        inv_and_padding_y(m_k, m, h_u, YY, p_r2c);
        FFT_x(g_u, XX, p_r2c);
        cov_IFFT(XX, YY, ZZ, m_k, z_u, p_c2r);

        trans(z_l, z_u, z, m, static_cast<int>(k_m_l), tag);
        for (int j=0;j<k_m_l;j++) {
            if (invalid_mask[j]) continue;
            kk1++;

            
            
            
            
            
            if(ll*k_m_l+j==target_index) {
                printf("z[%d]=%lf\n",j,z[j]);
                double t_temp[m];
                for(int tt = 0;tt < m;tt++)
                {
                    t_temp[tt] = (t[tt+j] - miu_t[j]) / si_t[j];
                }
                if(strcmp("msm",tag)==0) {
                    printf("msm%lf\n",msm(q_normal,m,t_temp,m,arg,w));
                }else if(strcmp("erp",tag)==0) {
                    printf("erp%lf\n",erp(q_normal,m,t_temp,m,arg,w));
                }else if(strcmp("lcss",tag)==0) {
                    printf("lcss%lf\n",lcss(q_normal,m,t_temp,m,w,arg));
                }else if(strcmp("edr",tag)==0) {
                    printf("edr%lf\n",edr(q_normal,m,t_temp,m,arg,w));
                }else if(strcmp("dtw",tag)==0) {
                    printf("dtw%lf\n",dtw(q_normal,m,t_temp,m,w));
                }
            }
            double t_temp[m];
            for (int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            double actual_dist = 0.0;
            if (strcmp("msm", tag) == 0) actual_dist = msm(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("dtw", tag) == 0) actual_dist = dtw(q_normal, m, t_temp, m, w);
            else if (strcmp("erp", tag) == 0) actual_dist = erp(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("edr", tag) == 0) actual_dist = edr(q_normal, m, t_temp, m, arg, w);
            else if (strcmp("lcss", tag) == 0) actual_dist = lcss(q_normal, m, t_temp, m, w, arg);

            if (!std::isfinite(actual_dist)) {
                nonfinite_distance_count++;
            } else if (actual_dist <= DIST_TOL) {
                zero_distance_count++;
            } else {
                double current_ratio = z[j] / actual_dist;
                max_ratio = std::max(max_ratio, current_ratio);
                if (current_ratio > 1.0 + DIST_TOL) {
                    bound_violation_count++;
                }
                ratio_count++;
                average_ratio += (current_ratio - average_ratio) / ratio_count;
            }
            if (z[j] <= threshold + DIST_TOL) {
                kk2++;
            }
        }
    }
    
        if(!flag_of_tail) {

            if(siz_t>2*m-2) {
                
                for (int f = 0; f < m - 1; f++) {
                    t[f] = t[k_m_l + f];  
                }
                
                input.read((char *)(t + m - 1), tail * sizeof(double));
            }else {
                input.read((char *)t , siz_t * sizeof(double));
            }

            mvmean(t,siz_t,m,miu_t,si_t);
            memset(invalid_mask, 0, sizeof(bool) * m_k);
            precompute_invalid_mask_fast(t, si_t, siz_t, m, invalid_mask);

            min_max_window(miu_t,siz_t,m,miu_min,miu_max);
            min_max_window(si_t,siz_t,m,si_min,si_max);

            lower_upper_lemire(t,siz_t,w,L,U);
            compute_normalized_arrays(U,  L, U_normal,  L_normal, siz_t,miu_max,miu_min, si_max,si_min);

            
            if(strcmp(tag,"dtw")==0) {
                t_u=find_tu_lp_dtw(U_normal,siz_t,q_normal,m,p,LB,m_p,tag,q_sorted);
                t_l=find_tl_lp_dtw(L_normal,siz_t,q_normal,m,p,LB,m_p,tag,q_sorted);
            }else if(strcmp(tag,"erp")==0) {
                t_u=find_tu_lp_erp(U_normal,siz_t,q_normal,m,p,arg,LB,m_p);
                t_l=find_tl_lp_erp(L_normal,siz_t,q_normal,m,p,arg,LB,m_p);
            }else {
                t_u=find_tu_lp(U_normal,siz_t,q_normal,m,p,arg,LB,m_p,tag,q_sorted);
                t_l=find_tl_lp(L_normal,siz_t,q_normal,m,p,arg,LB,m_p,tag,q_sorted);
            }

            memset(h_u, 0, sizeof(double) * m);
            memset(U, 0, sizeof(double) * siz_t);
            memset(g_u, 0, sizeof(double) * siz_t);

            memset(h_l, 0, sizeof(double) * m);
            memset(L, 0, sizeof(double) * siz_t);
            memset(g_l, 0, sizeof(double) * siz_t);

            memset(z_l, 0, sizeof(double) * siz_t);
            memset(z_u, 0, sizeof(double) * siz_t);
            memset(z, 0, sizeof(double) * tail);


            
            q_to_h_l(q_normal,m,t_l,arg,h_l,m,tag);
            q_to_h_u(q_normal,m,t_u,arg,h_u,m,tag);

            s_to_g_l(L_normal,siz_t,t_l, g_l, siz_t, tag, arg);
            s_to_g_u(U_normal,siz_t,t_u, g_u, siz_t, tag, arg);

            inv_and_padding_y(siz_t, m, h_l, YY_tail, p_r2c_tail);
            FFT_x(g_l, XX_tail, p_r2c_tail);
            cov_IFFT(XX_tail, YY_tail, ZZ_tail, siz_t, z_l, p_c2r_tail);

            inv_and_padding_y(siz_t, m, h_u, YY_tail, p_r2c_tail);
            FFT_x(g_u, XX_tail, p_r2c_tail);
            cov_IFFT(XX_tail, YY_tail, ZZ_tail, siz_t, z_u, p_c2r_tail);

            trans(z_l, z_u, z, m, static_cast<int>(tail), tag);
            for (int j=0;j<tail;j++) {
                if (invalid_mask[j]) continue;
                kk1++;

                double t_temp[m];
                for (int tt = 0; tt < m; tt++) {
                    t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
                }

                double actual_dist = 0.0;
                if (strcmp("msm", tag) == 0) actual_dist = msm(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("dtw", tag) == 0) actual_dist = dtw(q_normal, m, t_temp, m, w);
                else if (strcmp("erp", tag) == 0) actual_dist = erp(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("edr", tag) == 0) actual_dist = edr(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("lcss", tag) == 0) actual_dist = lcss(q_normal, m, t_temp, m, w, arg);

                if (!std::isfinite(actual_dist)) {
                    nonfinite_distance_count++;
                } else if (actual_dist <= DIST_TOL) {
                    zero_distance_count++;
                } else {
                    double current_ratio = z[j] / actual_dist;
                    max_ratio = std::max(max_ratio, current_ratio);
                    if (current_ratio > 1.0 + DIST_TOL) {
                        bound_violation_count++;
                    }
                    ratio_count++;
                    average_ratio += (current_ratio - average_ratio) / ratio_count;
                }
                if (z[j] <= threshold + DIST_TOL) {
                    kk2++;
                }
            }

        }

        std::cout << "kk1=" << kk1 << std::endl;
        std::cout << "kk2=" << kk2 << std::endl;
        std::cout << "Pruning Power (guolv)=" << (kk1 > 0 ? double(kk1-kk2)/kk1 : 0.0) << std::endl;
        std::cout << "Final dynamic threshold: " << threshold << std::endl;
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << " [Tightness Analysis Results]" << std::endl;
        std::cout << "  - Method:                LP_IP/" << tag << std::endl;
        std::cout << "  - TLB Formula:           mean(LP_IP_LB / ELD)" << std::endl;
        std::cout << "  - Total Subsequences:    " << kk1 << std::endl;
        std::cout << "  - Valid Ratios:          " << ratio_count << std::endl;
        std::cout << "  - Zero ELD Windows:      " << zero_distance_count << std::endl;
        std::cout << "  - Non-finite ELD Windows:" << nonfinite_distance_count << std::endl;
        std::cout << "  - LB Violations:         " << bound_violation_count << std::endl;
        std::cout << "  - Maximum LB/ELD Ratio:  " << std::fixed << std::setprecision(12) << max_ratio << std::endl;
        if (ratio_count > 0) {
            std::cout << "  - Average LB/Dist Ratio: " << std::fixed << std::setprecision(12) << average_ratio << std::endl;
        }
    std::cout << std::string(60, '=') << "\n" << std::endl;


        free(z);
        free(z_l);
        free(z_u);
        free(U);
        free(L);
        free(h_l);
        free(h_u);
        free(g_l);
        free(g_u);


        free(Z);
        free(Z_tail);

        free(m_p);
        free(LB);
        free(U_sort);
        free(q_sorted);
        free(L_sort);
        free(invalid_mask);

        fftw_free(XX);
        fftw_free(YY);
        fftw_free(ZZ);
        fftw_free(XX_tail);
        fftw_free(YY_tail);
        fftw_free(ZZ_tail);

        fftw_destroy_plan(p_c2r);
        fftw_destroy_plan(p_c2r_tail);
        fftw_destroy_plan(p_r2c);
        fftw_destroy_plan(p_r2c_tail);

}



void fsmeld_search_impl(const char* bi_t, const char* bi_q, int m,double threshold,char *tag,int w,double arg) {
    vector<int> A;
    MetricType metric = metric_from_tag(tag);
    int i;
    int m_k;
    double d;
    int p;
    if(strcmp(tag,"erp")==0) {
        p=2;
    }else if (strcmp(tag,"dtw")==0) {
        p=2;
    }else if (strcmp(tag,"msm")==0) {
        p=1;
    }else if (strcmp(tag,"edr")==0) {
        p=0;
    }else if (strcmp(tag,"lcss")==0) {
        p=0;
    }
    double ex_q=0,ex2_q=0,std=0,mean=0;

    unsigned long long kk1=0,kk2=0,kk3=0;

    streampos sizeofT;
    const string t_name = bi_t;
    ifstream input(t_name.c_str(),ios::binary | ios::in);
    input.seekg(0, ios::end);
    sizeofT = input.tellg();
    input.seekg(0, ios::beg);
    const long long available_points = sizeofT / sizeof(double);

    
    m_k = 1024;
    while(m_k <= 4*m) {
        m_k = m_k*2;  
    }
    if (available_points < m_k) {
        if (available_points < 2LL * m - 1) {
            std::cerr << "[Error] Target is too short for FSMELD: " << available_points << " points." << std::endl;
            return;
        }
        m_k = static_cast<int>(available_points);
    }
    double *t = (double *) malloc(sizeof(double) * m_k);
    double *q = (double *) malloc(sizeof(double)*m);
    double *q_normal = (double *) malloc(sizeof(double)*m);
    
    i=0;
    FILE *qp;
    qp = fopen(bi_q,"r");  
    while(fscanf(qp,"%lf",&d) != EOF && i < m)  
    {
        ex_q += d;
        ex2_q += d*d;
        q[i] = d;  
        i++;       
    }
    fclose(qp);    

    mean=ex_q/m;
    std=ex2_q/m;
    std=sqrt(std-mean*mean);
    fflush(stdout);

    
    for( i = 0 ; i < m ; i++ )
        q_normal[i] = (q[i] - mean)/std;

    double* q_sorted = (double*)malloc(sizeof(double) * m);
    memcpy(q_sorted, q_normal, sizeof(double) * m);
    std::sort(q_sorted, q_sorted + m);

    
    long long k_m_l = m_k - m + 1;  
    long long ts_n = (sizeofT/(sizeof(double))) - m + 1;  

    std::cout << "ts_n=" << ts_n << std::endl;
    std::cout << "m=" << m << std::endl;
    std::cout << "k_m_l=" << k_m_l << std::endl;
    std::cout << "m_k=" << m_k << std::endl;
    long long forsize = (ts_n)/(k_m_l);  
    int flag_of_tail = 0;
    double t_u,t_l;

    size_t tail = ts_n % k_m_l;
    int siz_t = tail + m - 1;  

    printf("tail=%zu\n",tail);
    if(tail==0)
    {
        flag_of_tail = 1;
    }
    if(siz_t<=2*m-2) {
        tail=m;
        siz_t=2*m-1;
    }
    printf("forsize=%lld\n",forsize);
    int forsize_divide_100 = forsize/100;
    int i_ = 0;
    
    double *Z = (double*) malloc(sizeof(double)*m_k);

    fftw_complex *XX = fftw_alloc_complex(m_k/2 + 1);
    fftw_complex *YY = fftw_alloc_complex(m_k/2 + 1);
    fftw_complex *ZZ = fftw_alloc_complex(m_k/2 + 1);

    fftw_plan p_r2c = fftw_plan_dft_r2c_1d(m_k, Z, ZZ, FFTW_ESTIMATE);
    fftw_plan p_c2r = fftw_plan_dft_c2r_1d(m_k, ZZ, Z, FFTW_ESTIMATE);

    double *Z_tail = (double*) malloc(sizeof(double)*siz_t);

    fftw_complex *XX_tail = fftw_alloc_complex(siz_t/2 + 1);
    fftw_complex *YY_tail = fftw_alloc_complex(siz_t/2 + 1);
    fftw_complex *ZZ_tail = fftw_alloc_complex(siz_t/2 + 1);

    fftw_plan p_r2c_tail = fftw_plan_dft_r2c_1d(siz_t, Z_tail, ZZ_tail, FFTW_ESTIMATE);
    fftw_plan p_c2r_tail = fftw_plan_dft_c2r_1d(siz_t, ZZ_tail, Z_tail, FFTW_ESTIMATE);

    
    double* L_sort=(double*)malloc(sizeof(double)*(m_k-2*m+2));
    double* U_sort=(double*)malloc(sizeof(double)*(m_k-2*m+2));
    double* LB=(double*)malloc(sizeof(double)*(m_k-2*m+2));
    double* m_p=(double*)malloc(sizeof(double)*(p+1));

    double *h_l = (double*) malloc(sizeof(double) * m);
    double *h_u = (double*) malloc(sizeof(double) * m);

    double *U = (double*) malloc(sizeof(double) * m_k);
    double *L = (double*) malloc(sizeof(double) * m_k);
    double *U_normal = (double*) malloc(sizeof(double) * m_k);
    double *L_normal = (double*) malloc(sizeof(double) * m_k);

    double *g_l = (double*) malloc(sizeof(double) * m_k);
    double *g_u = (double*) malloc(sizeof(double)* m_k);

    double *z_l = (double*) malloc(sizeof(double)*m_k);
    double *z_u = (double*) malloc(sizeof(double)*m_k);
    double *z = (double*) malloc(sizeof(double)*k_m_l);

    double *miu_t = (double*) malloc(sizeof(double)*k_m_l);
    double *si_t = (double*) malloc(sizeof(double)*k_m_l);

    double *miu_min = (double*) malloc(sizeof(double)*m_k);
    double *miu_max = (double*) malloc(sizeof(double)*m_k);
    double *si_min = (double*) malloc(sizeof(double)*m_k);
    double *si_max = (double*) malloc(sizeof(double)*m_k);
    bool* invalid_mask = (bool*)malloc(sizeof(bool) * m_k);
    double *U_Q = (double*) malloc(sizeof(double) * m);
    double *L_Q = (double*) malloc(sizeof(double) * m);
    lower_upper_lemire(q_normal, m, w, L_Q, U_Q);
    long long lp_ip_us_total = 0;

    
    
    for(long long ll = 0; ll < forsize; ll++) {
        if(forsize_divide_100) {
            if((ll%(forsize_divide_100)) == 0) {
                std::cout << "\b\b\b\b";
                std::cout << "|" << setw(3) << i_ << "%"<< std::endl;
                i_++;
                std::cout << "guolv=" << double(kk1-kk2)/kk1 << std::endl;
                std::cout << "guolv2=" << double(kk1-kk3)/kk1 << std::endl;
            }
        }

        if(ll != 0) {
            memcpy(t, t+k_m_l, sizeof(double)*m);
        } else {
            input.read((char *)t, (m-1)*sizeof(double));
        }
        
        input.read((char *)(t + m -1), (k_m_l)*sizeof(double));
        
        
        
        auto lp_ip_start = std::chrono::high_resolution_clock::now();


        mvmean(t,m_k,m,miu_t,si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, m_k, m, invalid_mask);

        min_max_window(miu_t,m_k,m,miu_min,miu_max);
        min_max_window(si_t,m_k,m,si_min,si_max);

        memset(h_u, 0, sizeof(double) * m);
        memset(U, 0, sizeof(double) * m_k);
        memset(g_u, 0, sizeof(double) * m_k);
        memset(U_normal, 0, sizeof(double) * m_k);

        memset(h_l, 0, sizeof(double) * m);
        memset(L, 0, sizeof(double) * m_k);
        memset(g_l, 0, sizeof(double) * m_k);
        memset(L_normal, 0, sizeof(double) * m_k);

        memset(z_l, 0, sizeof(double) * m_k);
        memset(z_u, 0, sizeof(double) * m_k);
        memset(z, 0, sizeof(double) * k_m_l);

        lower_upper_lemire(t,m_k,w,L,U);

        compute_normalized_arrays(U,  L, U_normal,  L_normal, m_k,miu_max,miu_min, si_max,si_min);

        if(strcmp(tag,"dtw")==0) {
            t_u=find_tu_lp_dtw(U_normal,m_k,q_normal,m,p,LB,m_p,tag,q_sorted);
            t_l=find_tl_lp_dtw(L_normal,m_k,q_normal,m,p,LB,m_p,tag,q_sorted);
        }else if(strcmp(tag,"erp")==0) {
            t_u=find_tu_lp_erp(U_normal,m_k,q_normal,m,p,arg,LB,m_p);
            t_l=find_tl_lp_erp(L_normal,m_k,q_normal,m,p,arg,LB,m_p);
        }else {
            t_u=find_tu_lp(U_normal,m_k,q_normal,m,p,arg,LB,m_p,tag,q_sorted);
            t_l=find_tl_lp(L_normal,m_k,q_normal,m,p,arg,LB,m_p,tag,q_sorted);
        }

        
        q_to_h_l(q_normal,m,t_l,arg,h_l,m,tag);
        q_to_h_u(q_normal,m,t_u,arg,h_u,m,tag);

        s_to_g_l(L_normal,m_k,t_l, g_l, m_k, tag, arg);
        s_to_g_u(U_normal,m_k,t_u, g_u, m_k, tag, arg);


        inv_and_padding_y(m_k, m, h_l, YY, p_r2c);

        FFT_x(g_l, XX, p_r2c);

        cov_IFFT(XX, YY, ZZ, m_k, z_l, p_c2r);

        inv_and_padding_y(m_k, m, h_u, YY, p_r2c);
        FFT_x(g_u, XX, p_r2c);
        cov_IFFT(XX, YY, ZZ, m_k, z_u, p_c2r);

        trans(z_l,z_u,z,m,k_m_l,tag);
        auto lp_ip_end = std::chrono::high_resolution_clock::now();
        lp_ip_us_total += std::chrono::duration_cast<std::chrono::microseconds>(lp_ip_end - lp_ip_start).count();
        for (int j=0,i=m-1;j<k_m_l;j++,i++) {
            if (invalid_mask[j]) {
                continue;
            }

            kk1++;
            if (z[j] <= threshold + DIST_TOL) {
                kk2++;
                double t_temp[m];
                for(int tt = 0;tt < m;tt++)
                {
                    t_temp[tt] = (t[tt+j] - miu_t[j]) / si_t[j];
                }
                double bsf = threshold;
                
                
                
                

















                double bglb_dist = bglb_elastic_lazy(metric, q_normal, m, t_temp, m, arg, L_Q, U_Q, w, bsf);

                
                if (bglb_dist <= threshold + DIST_TOL) {
                    kk3++;
                    double exact_dist = exact_elastic_distance(metric, q_normal, t_temp, m, w, arg);
                    if (!std::isfinite(exact_dist)) {
                        continue;
                    }
                    
                    
                    
                    if (exact_dist <= threshold + DIST_TOL) {
                        
                        printf("GETGET%s: %lf\n", tag, exact_dist);
                        printf("z[%d]=%lf\n", j, z[j]);
                        A.push_back(ll * k_m_l + j);
                    }
            }
            }
        }
    }
        if(!flag_of_tail) {

            if(siz_t>2*m-2) {
                
                for (int f = 0; f < m - 1; f++) {
                    t[f] = t[k_m_l + f];  
                }
                
                input.read((char *)(t + m - 1), tail * sizeof(double));
            }else {
                input.read((char *)t , siz_t * sizeof(double));
            }

            auto lp_ip_start = std::chrono::high_resolution_clock::now();
            mvmean(t,siz_t,m,miu_t,si_t);
            memset(invalid_mask, 0, sizeof(bool) * m_k);
            precompute_invalid_mask_fast(t, si_t, siz_t, m, invalid_mask);

            min_max_window(miu_t,siz_t,m,miu_min,miu_max);
            min_max_window(si_t,siz_t,m,si_min,si_max);

            lower_upper_lemire(t,siz_t,w,L,U);
            compute_normalized_arrays(U,  L, U_normal,  L_normal, siz_t,miu_max,miu_min, si_max,si_min);

            
            if(strcmp(tag,"dtw")==0) {
                t_u=find_tu_lp_dtw(U_normal,siz_t,q_normal,m,p,LB,m_p,tag,q_sorted);
                t_l=find_tl_lp_dtw(L_normal,siz_t,q_normal,m,p,LB,m_p,tag,q_sorted);
            }else if(strcmp(tag,"erp")==0) {
                t_u=find_tu_lp_erp(U_normal,siz_t,q_normal,m,p,arg,LB,m_p);
                t_l=find_tl_lp_erp(L_normal,siz_t,q_normal,m,p,arg,LB,m_p);
            }else {
                t_u=find_tu_lp(U_normal,siz_t,q_normal,m,p,arg,LB,m_p,tag,q_sorted);
                t_l=find_tl_lp(L_normal,siz_t,q_normal,m,p,arg,LB,m_p,tag,q_sorted);
            }

            memset(h_u, 0, sizeof(double) * m);
            memset(U, 0, sizeof(double) * siz_t);
            memset(g_u, 0, sizeof(double) * siz_t);

            memset(h_l, 0, sizeof(double) * m);
            memset(L, 0, sizeof(double) * siz_t);
            memset(g_l, 0, sizeof(double) * siz_t);

            memset(z_l, 0, sizeof(double) * siz_t);
            memset(z_u, 0, sizeof(double) * siz_t);
            memset(z, 0, sizeof(double) * tail);


            
            q_to_h_l(q_normal,m,t_l,arg,h_l,m,tag);
            q_to_h_u(q_normal,m,t_u,arg,h_u,m,tag);

            s_to_g_l(L_normal,siz_t,t_l, g_l, siz_t, tag, arg);
            s_to_g_u(U_normal,siz_t,t_u, g_u, siz_t, tag, arg);

            inv_and_padding_y(siz_t, m, h_l, YY_tail, p_r2c_tail);
            FFT_x(g_l, XX_tail, p_r2c_tail);
            cov_IFFT(XX_tail, YY_tail, ZZ_tail, siz_t, z_l, p_c2r_tail);

            inv_and_padding_y(siz_t, m, h_u, YY_tail, p_r2c_tail);
            FFT_x(g_u, XX_tail, p_r2c_tail);
            cov_IFFT(XX_tail, YY_tail, ZZ_tail, siz_t, z_u, p_c2r_tail);

            for (int j=0,i=m-1;j<tail;j++,i++) {
                if(strcmp(tag,"lcss")==0) {
                    z[j] = std::max(0.0, (z_l[i] + z_u[i]) / m);
                }else if(strcmp(tag,"dtw")==0 || strcmp(tag,"erp")==0) {
                    z[j]=std::sqrt(std::max(0.0, z_l[i]+z_u[i]));
                }else {
                    z[j]=z_l[i]+z_u[i];
                }
            }
            auto lp_ip_end = std::chrono::high_resolution_clock::now();
            lp_ip_us_total += std::chrono::duration_cast<std::chrono::microseconds>(lp_ip_end - lp_ip_start).count();

            for (int j=0,i=m-1;j<tail;j++,i++) {
                if (invalid_mask[j]) {
                    continue;
                }
                kk1++;
                if (z[j] <= threshold + DIST_TOL) {
                    kk2++;

                    double t_temp[m];
                    for(int tt = 0;tt < m;tt++)
                {
                    t_temp[tt] = (t[tt+j] - miu_t[j]) / si_t[j];
                }
                
                    double bsf = threshold;
                
                
                

















                double bglb_dist = bglb_elastic_lazy(metric, q_normal, m, t_temp, m, arg, L_Q, U_Q, w, bsf);

                
                if (bglb_dist <= threshold + DIST_TOL) {
                    kk3++;
                    double exact_dist = exact_elastic_distance(metric, q_normal, t_temp, m, w, arg);
                    if (!std::isfinite(exact_dist)) {
                        continue;
                    }

                    
                    if (exact_dist <= threshold + DIST_TOL) {
                        
                        printf("%s: %lf\n", tag, exact_dist);
                        printf("z[%d]=%lf\n", j, z[j]);
                        A.push_back(forsize * k_m_l + j);
                    }
                }
                }
            }

        }

        std::cout << "kk1=" << kk1 << std::endl;
        std::cout << "kk2=" << kk2 << std::endl;
        std::cout << "kk2=" << kk2 << std::endl;
        std::cout << "Pruning Power (guolv)=" << double(kk1-kk2)/kk1 << std::endl;
        std::cout << "Pruning Power (guolv2)=" << double(kk1-kk3)/kk1 << std::endl;
        std::cout << "LP_IP Lower Bound Time taken: " << (lp_ip_us_total / 1000.0) << " ms" << std::endl;
        std::cout << "LP_IP Lower Bound Time taken: " << (lp_ip_us_total / 1000000.0) << " s" << std::endl;
        std::cout << "Final dynamic threshold: " << threshold << std::endl;
        std::cout << "Ablation Total Valid Candidates: " << kk1 << std::endl;
        std::cout << "Ablation Candidates after LP/IP: " << kk2 << std::endl;
        std::cout << "Ablation Candidates after BGLB: " << kk3 << std::endl;
        std::cout << "Ablation Final Matches: " << A.size() << std::endl;

        for (size_t i = 0; i < A.size(); ++i) {
            printf("A[%d]=%d\n",i,A[i]);
        }
        free(invalid_mask);
        free(z);
        free(z_l);
        free(z_u);
        free(U);
        free(L);
        free(h_l);
        free(h_u);
        free(g_l);
        free(g_u);


        free(Z);
        free(Z_tail);

        free(U_Q);
        free(L_Q);

        free(m_p);
        free(LB);
        free(U_sort);
        free(L_sort);
        free(q_sorted);

        fftw_free(XX);
        fftw_free(YY);
        fftw_free(ZZ);
        fftw_free(XX_tail);
        fftw_free(YY_tail);
        fftw_free(ZZ_tail);

        fftw_destroy_plan(p_c2r);
        fftw_destroy_plan(p_c2r_tail);
        fftw_destroy_plan(p_r2c);
        fftw_destroy_plan(p_r2c_tail);

}

static void fsmeld_direct_impl(const char* bi_t, const char* bi_q, int m, double threshold,
                              char *tag, int w, double arg, bool use_bglb, int exact_top_k) {
    vector<long long> A;
    vector<pair<double, long long>> exact_results;
    priority_queue<double> exact_topk_heap;
    const bool derive_exact_topk = !use_bglb && exact_top_k > 0;
    MetricType metric = metric_from_tag(tag);
    int i;
    double d;

    double ex_q = 0, ex2_q = 0, std_val = 0, mean = 0;
    unsigned long long kk1 = 0; 
    unsigned long long kk2 = 0; 

    auto record_exact_result = [&](double distance, long long index) {
        if (derive_exact_topk) {
            exact_results.emplace_back(distance, index);
            if (exact_topk_heap.size() < static_cast<size_t>(exact_top_k)) {
                exact_topk_heap.push(distance);
            } else if (distance < exact_topk_heap.top()) {
                exact_topk_heap.pop();
                exact_topk_heap.push(distance);
            }
        }
    };

    streampos sizeofT;
    const string t_name = bi_t;
    ifstream input(t_name.c_str(), ios::binary | ios::in);
    input.seekg(0, ios::end);
    sizeofT = input.tellg();
    input.seekg(0, ios::beg);
    const long long available_points = sizeofT / sizeof(double);

    
    int m_k = 1024;
    while(m_k <= 4 * m) {
        m_k = m_k * 2;
    }

    if (available_points < m_k) {
        if (available_points < 2LL * m - 1) {
            std::cerr << "[Error] Target is too short for direct ELD scan: " << available_points << " points." << std::endl;
            return;
        }
        m_k = static_cast<int>(available_points);
    }

    double *t = (double *) malloc(sizeof(double) * m_k);
    double *q = (double *) malloc(sizeof(double) * m);
    double *q_normal = (double *) malloc(sizeof(double) * m);

    i = 0;
    FILE *qp = fopen(bi_q, "r");
    while(fscanf(qp, "%lf", &d) != EOF && i < m) {
        ex_q += d;
        ex2_q += d * d;
        q[i] = d;
        i++;
    }
    fclose(qp);

    mean = ex_q / m;
    std_val = ex2_q / m;
    std_val = sqrt(std_val - mean * mean);
    fflush(stdout);

    for(i = 0; i < m; i++) {
        q_normal[i] = (q[i] - mean) / std_val;
    }


    double *U_Q = (double*) malloc(sizeof(double) * m);
    double *L_Q = (double*) malloc(sizeof(double) * m);
    if (use_bglb) {
        lower_upper_lemire(q_normal, m, w, L_Q, U_Q);
    }

    long long k_m_l = m_k - m + 1;  
    long long ts_n = (sizeofT / (sizeof(double))) - m + 1;  

    std::cout << "ts_n=" << ts_n << std::endl;
    std::cout << "m=" << m << std::endl;
    std::cout << "k_m_l=" << k_m_l << std::endl;
    std::cout << "m_k=" << m_k << std::endl;

    long long forsize = (ts_n) / (k_m_l);
    int flag_of_tail = 0;

    size_t tail = ts_n % k_m_l;
    int siz_t = tail + m - 1;

    printf("tail=%zu\n", tail);
    if(tail == 0) flag_of_tail = 1;
    if(siz_t <= 2 * m - 2) {
        tail = m;
        siz_t = 2 * m - 1;
    }

    printf("forsize=%lld\n", forsize);
    int forsize_divide_100 = forsize / 100;
    int i_ = 0;

    double *miu_t = (double*) malloc(sizeof(double) * k_m_l);
    double *si_t = (double*) malloc(sizeof(double) * k_m_l);
    bool* invalid_mask = (bool*)malloc(sizeof(bool) * m_k);
    for(long long ll = 0; ll < forsize; ll++) {
        if(forsize_divide_100 && (ll % forsize_divide_100) == 0) {
            std::cout << "\b\b\b\b|" << setw(3) << i_ << "% ";
            i_++;
            if(kk1 > 0 && use_bglb) {
                std::cout << " BGLB Pruning Rate=" << double(kk1 - kk2) / kk1 << "\r" << std::flush;
            }
        }

        if(ll != 0) {
            memcpy(t, t + k_m_l, sizeof(double) * m);
        } else {
            input.read((char *)t, (m - 1) * sizeof(double));
        }
        input.read((char *)(t + m - 1), (k_m_l) * sizeof(double));




        mvmean(t, m_k, m, miu_t, si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, m_k, m, invalid_mask);

        for (int j = 0; j < k_m_l; j++) {
            if (invalid_mask[j]) {
                continue;
            }
            kk1++;

            
            double t_temp[m];
            for(int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            bool passed_bglb = true;
            if (use_bglb) {
                double bsf = threshold;
                double bglb_dist = bglb_elastic_lazy(
                    metric, q_normal, m, t_temp, m, arg, L_Q, U_Q, w, bsf);
                passed_bglb = bglb_dist <= threshold + DIST_TOL;
            }

            if (passed_bglb) {
                kk2++;
                double exact_dist = 0.0;

                
                if (strcmp("msm", tag) == 0) exact_dist = msm(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("erp", tag) == 0) exact_dist = erp(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("lcss", tag) == 0) exact_dist = lcss(q_normal, m, t_temp, m, w, arg);
                else if (strcmp("edr", tag) == 0) exact_dist = edr(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("dtw", tag) == 0) exact_dist = dtw(q_normal, m, t_temp, m, w);

                if (!std::isfinite(exact_dist)) {
                    continue;
                }
                record_exact_result(exact_dist, ll * k_m_l + j);
                if (exact_dist <= threshold + DIST_TOL) {
                    printf("\nGET %s: %lf (index: %lld)\n", tag, exact_dist, ll * k_m_l + j);
                    A.push_back(ll * k_m_l + j);
                }
            }
        }
    }
    
    if(!flag_of_tail) {
        if(siz_t > 2 * m - 2) {
            for (int f = 0; f < m - 1; f++) {
                t[f] = t[k_m_l + f];
            }
            input.read((char *)(t + m - 1), tail * sizeof(double));
        } else {
            input.read((char *)t, siz_t * sizeof(double));
        }

        mvmean(t, siz_t, m, miu_t, si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, siz_t, m, invalid_mask);

        for (int j = 0; j < tail; j++) {
            if (invalid_mask[j]) {
                continue;
            }
            kk1++;
            double t_temp[m];
            for(int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            bool passed_bglb = true;
            if (use_bglb) {
                double bsf = threshold;
                double bglb_dist = bglb_elastic_lazy(
                    metric, q_normal, m, t_temp, m, arg, L_Q, U_Q, w, bsf);
                passed_bglb = bglb_dist <= threshold + DIST_TOL;
            }

            if (passed_bglb) {
                kk2++;
                double exact_dist = 0.0;

                if (strcmp("msm", tag) == 0) exact_dist = msm(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("erp", tag) == 0) exact_dist = erp(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("lcss", tag) == 0) exact_dist = lcss(q_normal, m, t_temp, m, w, arg);
                else if (strcmp("edr", tag) == 0) exact_dist = edr(q_normal, m, t_temp, m, arg, w);
                else if (strcmp("dtw", tag) == 0) exact_dist = dtw(q_normal, m, t_temp, m, w);

                if (!std::isfinite(exact_dist)) {
                    continue;
                }
                record_exact_result(exact_dist, forsize * k_m_l + j);
                if (exact_dist <= threshold + DIST_TOL) {
                    printf("\nGET %s: %lf (index: %lld)\n", tag, exact_dist, forsize * k_m_l + j);
                    A.push_back(forsize * k_m_l + j);
                }
            }
        }
    }

    if (derive_exact_topk) {
        A.clear();
        if (exact_topk_heap.size() == static_cast<size_t>(exact_top_k)) {
            threshold = exact_topk_heap.top();
            for (const auto& result : exact_results) {
                if (result.first <= threshold + DIST_TOL) {
                    A.push_back(result.second);
                }
            }
            std::sort(A.begin(), A.end());
            std::cout << std::setprecision(17)
                      << "Ablation Exact Top-K Epsilon: " << threshold << std::endl;
        } else {
            std::cout << "Ablation Exact Top-K Epsilon: nan" << std::endl;
            std::cerr << "[Error] Only " << exact_topk_heap.size()
                      << " finite ELD values; cannot compute Top-" << exact_top_k << "." << std::endl;
        }
    }

    std::cout << "\n\n--- Statistics ---" << std::endl;
    std::cout << "Total comparisons (kk1) = " << kk1 << std::endl;
    if (use_bglb) {
        std::cout << "Candidates passed BGLB (kk2) = " << kk2 << std::endl;
        std::cout << "BGLB Pruning Rate = " << (kk1 == 0 ? 0.0 : double(kk1 - kk2) / kk1) << std::endl;
    } else {
        std::cout << "Exact ELD evaluations = " << kk2 << std::endl;
    }
    std::cout << "Final matches found: " << A.size() << std::endl;
    std::cout << "Ablation Total Valid Candidates: " << kk1 << std::endl;
    if (use_bglb) {
        std::cout << "Ablation Candidates after BGLB: " << kk2 << std::endl;
    } else {
        std::cout << "Ablation Exact ELD Evaluations: " << kk2 << std::endl;
    }
    std::cout << "Ablation Final Matches: " << A.size() << std::endl;

    for (size_t k = 0; k < A.size(); ++k) {
        printf("A[%zu] = %lld\n", k, A[k]);
    }

    
    free(t);
    free(q);
    free(q_normal);
    free(miu_t);
    free(si_t);
    free(U_Q);
    free(L_Q);
    free(invalid_mask);
}

void fsmeld_raw_impl(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w, double arg) {
    fsmeld_direct_impl(bi_t, bi_q, m, threshold, tag, w, arg, true, 0);
}

void fsmeld_exact_impl(const char* bi_t, const char* bi_q, int m, double threshold,
                  char *tag, int w, double arg, int top_k) {
    fsmeld_direct_impl(bi_t, bi_q, m, threshold, tag, w, arg, false, top_k);
}

static void print_candidate_exclusion_summary(const ScopedCandidateExclusion& exclusion,
                                              long long exclude_start, long long exclude_end) {
    if (exclude_start < 0 || exclude_end <= exclude_start) return;
    std::cout << "Candidate Exclusion Start: " << exclude_start << std::endl;
    std::cout << "Candidate Exclusion End Exclusive: " << exclude_end << std::endl;
    std::cout << "Candidate Exclusion Windows Visited: " << exclusion.visited_windows() << std::endl;
    std::cout << "Candidate Exclusion Windows Skipped: " << exclusion.excluded_windows() << std::endl;
}

void FSMELD(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w, double arg,
            long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    fsmeld_search_impl(bi_t, bi_q, m, threshold, tag, w, arg);
    print_candidate_exclusion_summary(exclusion, exclude_start, exclude_end);
}

void FSMELD_epsilon(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w, double arg,
                    int top_k, long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    FSMELD_epsilon(bi_t, bi_q, m, threshold, tag, w, arg, top_k);
    print_candidate_exclusion_summary(exclusion, exclude_start, exclude_end);
}

void FSMELD_raw(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w, double arg,
                long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    fsmeld_raw_impl(bi_t, bi_q, m, threshold, tag, w, arg);
    print_candidate_exclusion_summary(exclusion, exclude_start, exclude_end);
}

void FSMELD_exact(const char* bi_t, const char* bi_q, int m, double threshold,
                  char *tag, int w, double arg, int top_k,
                  long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    fsmeld_exact_impl(bi_t, bi_q, m, threshold, tag, w, arg, top_k);
    print_candidate_exclusion_summary(exclusion, exclude_start, exclude_end);
}

void FSMELD_find_threshold(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w,
                           double arg, int top_k, long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    find_threshold(bi_t, bi_q, m, threshold, tag, w, arg, top_k);
    print_candidate_exclusion_summary(exclusion, exclude_start, exclude_end);
}


void FSMELD_get_tlb(const char* bi_t, const char* bi_q, int m, double threshold, char *tag, int w,
                    double arg, double target_index,
                    long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    get_tlb(bi_t, bi_q, m, threshold, tag, w, arg, target_index);
    print_candidate_exclusion_summary(exclusion, exclude_start, exclude_end);
}

