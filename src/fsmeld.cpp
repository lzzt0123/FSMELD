#include <stdexcept>
//
// Created by LZT on 2025-08-09.
//
#include <stdlib.h>
#include <string.h>
#include "fstream"
#include<iostream>
#include "def.h"
#include "funcs.h"
#include "search_result.h"
#include "fsmeld_internal.h"
#include "fftw3.h"
#include <cmath>
#include <algorithm>
#include <cstring>   
#include <cfloat> 
#include <iomanip> 
#include <chrono>
#include <cerrno>
#include <climits>
#include <set>
#include <sstream>
#include <memory>

using namespace std;

static bool resolve_fsmeld_macro_length(int m, long long available_points,
                                        int* macro_length) {
    long long selected = 1;
    while (selected < 8LL * m) selected *= 2;
    if (available_points < 2LL * m - 1) return false;
    if (available_points < selected) selected = available_points;
    if (selected > INT_MAX) throw std::runtime_error("Macro segment exceeds the supported FFT length.");
    *macro_length = static_cast<int>(selected);
    return true;
}

struct CandidateExclusionState {
    bool active = false;
    long long begin = -1;
    long long end = -1;
    long long next_window_start = 0;
    long long excluded_windows = 0;
};

static thread_local CandidateExclusionState candidate_exclusion;

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
    }else if (strcmp(tag,"edr")==0 || strcmp(tag,"swale")==0) {
        for(int i=0,j=m-1;i<k_m_l;i++,j++) {
            z[i]=z_l[j]+z_u[j];
        }
    }else if (strcmp(tag,"lcss")==0) {
        for(int i=0,j=m-1;i<k_m_l;i++,j++) {
            z[i] = std::max(0.0, (z_l[j] + z_u[j]) / m);
        }
    }
}

// Mark windows that cannot be z-normalized safely: NaN/Inf input or near-zero sigma.

void precompute_invalid_mask_fast(const double* t, const double* si_t,
                                         int segment_size, int m, bool* invalid_mask,
                                         int first_processed_window,
                                         int processed_window_count) {
    const int window_count = segment_size - m + 1;
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

    if (processed_window_count < 0) {
        processed_window_count = window_count - first_processed_window;
    }
    if (first_processed_window < 0 || processed_window_count < 0 ||
        first_processed_window + processed_window_count > window_count) {
        return;
    }

    if (candidate_exclusion.active) {
        const long long block_start = candidate_exclusion.next_window_start;
        const long long overlap_begin = std::max(
            block_start, candidate_exclusion.begin);
        const long long overlap_end = std::min(
            block_start + processed_window_count, candidate_exclusion.end);
        for (long long global_index = overlap_begin; global_index < overlap_end; ++global_index) {
            invalid_mask[first_processed_window + global_index - block_start] = true;
        }
        if (overlap_end > overlap_begin) {
            candidate_exclusion.excluded_windows += overlap_end - overlap_begin;
        }
        candidate_exclusion.next_window_start += processed_window_count;
    }
}

// Load an EOF-aligned segment. Tail handling reads one complete macro segment
// and evaluates only its final, previously unprocessed windows. This preserves
// the exact n-m+1 candidate set and keeps threshold optimization well-defined
// even when the remainder contains fewer than m windows.
bool read_tail_segment(std::ifstream& input, double* t, int segment_size) {
    if (segment_size <= 0) return false;
    const std::streamoff byte_count =
        static_cast<std::streamoff>(segment_size) * sizeof(double);
    input.clear();
    input.seekg(-byte_count, std::ios::end);
    if (!input) return false;
    input.read(reinterpret_cast<char*>(t), byte_count);
    return input.gcount() == byte_count;
}

double exact_elastic_distance(MetricType metric, double* q, double* t,
                                     int m, int w, double arg) {
    switch (metric) {
        case METRIC_MSM: return msm(q, m, t, m, arg, w);
        case METRIC_ERP: return erp(q, m, t, m, arg, w);
        case METRIC_LCSS: return lcss(q, m, t, m, w, arg);
        case METRIC_EDR: return edr(q, m, t, m, arg, w);
        case METRIC_DTW: return dtw(q, m, t, m, w);
        case METRIC_SWALE:
            return swale(q, m, t, m, SWALE_MISMATCH_COST,
                         SWALE_MATCH_COST, arg, w);
        default: return NAN;
    }
}

static void find_threshold_impl(const char* bi_t, const char* bi_q, int m,
                                double threshold, char *tag, int w,
                                double arg, int top_k,
                                bool filter_by_initial_epsilon) {
    
    priority_queue<pair<double, long long>> top_k_heap;
    vector<long long> A; 
    MetricType metric = metric_from_tag(tag);

    int i;
    int m_k;
    double d;
    int p;
    const double swale_r = SWALE_MATCH_COST;
    const double swale_p = SWALE_MISMATCH_COST;

    if (strcmp(tag, "erp") == 0) p = 2;
    else if (strcmp(tag, "dtw") == 0) p = 2;
    else if (strcmp(tag, "msm") == 0) p = 1;
    else if (strcmp(tag, "edr") == 0) p = 0;
    else if (strcmp(tag, "lcss") == 0) p = 0;
    else if (strcmp(tag, "swale") == 0) p = 0;

    double ex_q = 0, ex2_q = 0, std_val = 0, mean = 0;
    unsigned long long kk1 = 0; 
    unsigned long long exact_count = 0;
    unsigned long long lb_pruned = 0;
    unsigned long long glb_pruned = 0;

    std::ifstream sizing_input(bi_t, std::ios::binary | std::ios::ate);
    const long long available_points = static_cast<long long>(sizing_input.tellg()) / sizeof(double);
    if (!resolve_fsmeld_macro_length(m, available_points, &m_k))
        throw std::runtime_error("Target is too short for cutoff search.");
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
    for (i = 0; i < m; i++) q_normal[i] = (q[i] - mean) / std_val;

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
    int siz_t = m_k;
    const int tail_first_window = static_cast<int>(k_m_l - tail);
    const int tail_end_window = tail_first_window + static_cast<int>(tail);

    if (tail == 0) flag_of_tail = 1;

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
//                std::cout << "guolv=" << (kk1 == 0 ? 0.0 : double(kk1 - lb_pruned) / kk1) << std::endl;
//                std::cout << "guolv2=" << (kk1 == 0 ? 0.0 : double(kk1 - lb_pruned - glb_pruned) / kk1) << std::endl;
                std::cout << "Threshold=" << std::setprecision(17) << threshold << std::endl;
            }
        }

        if (ll != 0) memcpy(t, t + k_m_l, sizeof(double) * (m - 1));
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
        compute_normalized_arrays(U, L, U_normal, L_normal, m_k,
                                  miu_max, miu_min, si_max, si_min);

        if (strcmp(tag, "dtw") == 0) {
            t_u = find_tu_lp_dtw(U_normal, m_k, q_normal, m, p, LB, m_p, tag, q_sorted);
            t_l = find_tl_lp_dtw(L_normal, m_k, q_normal, m, p, LB, m_p, tag, q_sorted);
        } else if (strcmp(tag, "swale") == 0) {
            t_u = find_tu_lp_swale(U_normal, m_k, q_normal, m, p, arg, LB, m_p, tag, swale_r, swale_p, q_sorted);
            t_l = find_tl_lp_swale(L_normal, m_k, q_normal, m, p, arg, LB, m_p, tag, swale_r, swale_p, q_sorted);
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
            double active_threshold =
                ((filter_by_initial_epsilon || heap_full))
                    ? threshold : DBL_MAX;
            if (fsmeld_lbip_exceeds_cutoff(z[j], active_threshold, metric)) {
                lb_pruned++;
                continue;
            }

            for (int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            if ((filter_by_initial_epsilon || heap_full)) {
                double bsf = threshold;
                double glb_dist = glb_elastic(metric, q_normal, m, t_temp, m, arg, L_Q, U_Q, w, bsf);
                if (glb_dist > threshold) {
                    glb_pruned++;
                    continue;
                }
            }

            const double actual_dist = exact_elastic_distance(
                metric, q_normal, t_temp, m, w, arg);

            exact_count++;
            if (!std::isfinite(actual_dist)) {
                continue;
            }
            if (filter_by_initial_epsilon && actual_dist > threshold) {
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
        if (!read_tail_segment(input, t, siz_t)) {
            std::cerr << "[Error] Failed to read the exact tail segment." << std::endl;
            return;
        }

        mvmean(t, siz_t, m, miu_t, si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, siz_t, m, invalid_mask,
                                     tail_first_window, static_cast<int>(tail));
        min_max_window(miu_t, siz_t, m, miu_min, miu_max);
        min_max_window(si_t, siz_t, m, si_min, si_max);

        lower_upper_lemire(t, siz_t, w, L, U);
        compute_normalized_arrays(U, L, U_normal, L_normal, siz_t,
                                  miu_max, miu_min, si_max, si_min);

        if (strcmp(tag, "dtw") == 0) {
            t_u = find_tu_lp_dtw(U_normal, siz_t, q_normal, m, p, LB, m_p, tag, q_sorted);
            t_l = find_tl_lp_dtw(L_normal, siz_t, q_normal, m, p, LB, m_p, tag, q_sorted);
        } else if (strcmp(tag, "swale") == 0) {
            t_u = find_tu_lp_swale(U_normal, siz_t, q_normal, m, p, arg, LB, m_p, tag, swale_r, swale_p, q_sorted);
            t_l = find_tl_lp_swale(L_normal, siz_t, q_normal, m, p, arg, LB, m_p, tag, swale_r, swale_p, q_sorted);
        } else if(strcmp(tag,"erp")==0) {
            t_u=find_tu_lp_erp(U_normal,siz_t,q_normal,m,p,arg,LB,m_p);
            t_l=find_tl_lp_erp(L_normal,siz_t,q_normal,m,p,arg,LB,m_p);
        }else {
            t_u = find_tu_lp(U_normal, siz_t, q_normal, m, p, arg, LB, m_p, tag, q_sorted);
            t_l = find_tl_lp(L_normal, siz_t, q_normal, m, p, arg, LB, m_p, tag, q_sorted);
        }

        memset(h_u, 0, sizeof(double) * m); memset(U, 0, sizeof(double) * siz_t); memset(g_u, 0, sizeof(double) * siz_t);
        memset(h_l, 0, sizeof(double) * m); memset(L, 0, sizeof(double) * siz_t); memset(g_l, 0, sizeof(double) * siz_t);
        memset(z_l, 0, sizeof(double) * siz_t); memset(z_u, 0, sizeof(double) * siz_t); memset(z, 0, sizeof(double) * k_m_l);

        q_to_h_l(q_normal, m, t_l, arg, h_l, m, tag); q_to_h_u(q_normal, m, t_u, arg, h_u, m, tag);
        s_to_g_l(L_normal, siz_t, t_l, g_l, siz_t, tag, arg); s_to_g_u(U_normal, siz_t, t_u, g_u, siz_t, tag, arg);

        inv_and_padding_y(siz_t, m, h_l, YY_tail, p_r2c_tail); FFT_x(g_l, XX_tail, p_r2c_tail); cov_IFFT(XX_tail, YY_tail, ZZ_tail, siz_t, z_l, p_c2r_tail);
        inv_and_padding_y(siz_t, m, h_u, YY_tail, p_r2c_tail); FFT_x(g_u, XX_tail, p_r2c_tail); cov_IFFT(XX_tail, YY_tail, ZZ_tail, siz_t, z_u, p_c2r_tail);

        for (int j = tail_first_window, i = m - 1 + tail_first_window;
             j < tail_end_window; j++, i++) {
            if (invalid_mask[j]) {
                continue;
            }
            if (strcmp("lcss", tag) == 0) z[j] = std::max(0.0, (z_l[i] + z_u[i]) / m);
            else if (strcmp("dtw", tag) == 0 || strcmp("erp", tag) == 0) z[j] = std::sqrt(std::max(0.0, z_l[i] + z_u[i]));
            else z[j] = z_l[i] + z_u[i];

            kk1++;
            bool heap_full = top_k_heap.size() >= static_cast<size_t>(top_k);
            double active_threshold =
                ((filter_by_initial_epsilon || heap_full))
                    ? threshold : DBL_MAX;
            if (fsmeld_lbip_exceeds_cutoff(z[j], active_threshold, metric)) {
                lb_pruned++;
                continue;
            }

            for (int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            if ((filter_by_initial_epsilon || heap_full)) {
                double bsf = threshold;
                double glb_dist = glb_elastic(metric, q_normal, m, t_temp, m, arg, L_Q, U_Q, w, bsf);
                if (glb_dist > threshold) {
                    glb_pruned++;
                    continue;
                }
            }

            const double actual_dist = exact_elastic_distance(
                metric, q_normal, t_temp, m, w, arg);

            exact_count++;
            if (!std::isfinite(actual_dist)) {
                continue;
            }
            if (filter_by_initial_epsilon && actual_dist > threshold) {
                continue;
            }
            long long current_index = forsize * k_m_l + (j - tail_first_window);
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

    if (top_k_heap.size() < static_cast<size_t>(top_k))
        throw std::runtime_error("Fewer valid, non-excluded candidates than requested k.");
    std::cout << "\n\n" << std::string(60, '=') << std::endl;
    std::cout << " [FSMELD Top-K Threshold Results]" << std::endl;
    std::cout << "  - Candidates checked by LB:   " << kk1 << std::endl;
    std::cout << "  - Pruned by LBIP:             " << lb_pruned << std::endl;
    std::cout << "  - Pruned by GLB:              " << glb_pruned << std::endl;
    std::cout << "  - Exact distances computed:   " << exact_count << std::endl;
    std::cout << std::setprecision(17);
    std::cout << "  - Final Top-K Threshold:      " << threshold << std::endl;
    if (filter_by_initial_epsilon) {
        std::cout << "  - Final FSMELD epsilon:       " << threshold << std::endl;
        const bool topk_complete =
            top_k_heap.size() >= static_cast<size_t>(top_k);
        std::cout << "TOP_K_THRESHOLD_RESULT status="
                  << (topk_complete ? "OK" : "FAILED")
                  << " retained=" << top_k_heap.size()
                  << " requested=" << top_k
                  << " initial_filter=enabled"
                  << " final_epsilon=" << threshold
                  << " exact_evaluations=" << exact_count
                  << std::endl;
    }
    std::cout << std::string(60, '=') << "\n" << std::endl;

    {
        
        while(!top_k_heap.empty()) {
            A.push_back(top_k_heap.top().second);
            top_k_heap.pop();
        }
        std::reverse(A.begin(), A.end()); 
        for (size_t k = 0; k < A.size(); ++k) {
            std::cout << "Rank " << k+1 << ": Index A[" << k+1 << "]=" << A[k] << std::endl;
        }
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
static void FSMDTW_impl(const char* bi_t, const char* bi_q, int m,
                        double threshold, char *tag, int w, double arg) {
    using ExperimentClock = std::chrono::steady_clock;
    const auto total_start = std::chrono::steady_clock::now();
    vector<long long> A;
    struct FSMELDSurvivor {
        long long index;
        std::vector<double> normalized;
    };
    MetricType metric = metric_from_tag(tag);
    int i;
    int m_k;
    double d;
    int p;
    const double swale_r = SWALE_MATCH_COST;
    const double swale_p = SWALE_MISMATCH_COST;
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
    }else if (strcmp(tag,"swale")==0) {
        p=0;
    }

    double ex_q=0,ex2_q=0,std=0,mean=0;

    unsigned long long kk1=0,kk2=0,kk3=0;
    const double evaluation_threshold = threshold;

    streampos sizeofT;
    const string t_name = bi_t;
    ifstream input(t_name.c_str(),ios::binary | ios::in);
    input.seekg(0, ios::end);
    sizeofT = input.tellg();
    input.seekg(0, ios::beg);
    const long long available_points = sizeofT / sizeof(double);

    // Resolve the macro-segment/FFT length.  For sensitivity experiments set
    // FSMELD_MACRO_LENGTH; otherwise the original automatic policy is used.
    if (!resolve_fsmeld_macro_length(m, available_points, &m_k)) {
        return;
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

    const auto query_normalization_start = ExperimentClock::now();
    mean=ex_q/m;
    std=ex2_q/m;
    std=sqrt(std-mean*mean);
    fflush(stdout);
    for(i = 0; i < m; i++) q_normal[i] = (q[i] - mean)/std;

    double* q_sorted = (double*)malloc(sizeof(double) * m);
    memcpy(q_sorted, q_normal, sizeof(double) * m);
    std::sort(q_sorted, q_sorted + m);
    long long lb_ns_total = std::chrono::duration_cast<std::chrono::nanoseconds>(
        ExperimentClock::now() - query_normalization_start).count();

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
    int siz_t = m_k;
    const int tail_first_window = static_cast<int>(k_m_l - tail);
    const int tail_end_window = tail_first_window + static_cast<int>(tail);

    printf("tail=%zu\n",tail);
    if(tail==0)
    {
        flag_of_tail = 1;
    }
    printf("forsize=%lld\n",forsize);
    int forsize_divide_100 = forsize/100;
    int i_ = 0;
    //
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
    const auto query_envelope_start = ExperimentClock::now();
    lower_upper_lemire(q_normal, m, w, L_Q, U_Q);
    lb_ns_total += std::chrono::duration_cast<std::chrono::nanoseconds>(
        ExperimentClock::now() - query_envelope_start).count();
    long long lp_ip_us_total = 0;
    unsigned long long invalid_candidates = 0;

    //int sample_rate =1;
    //for(long long ll = 0; ll < 10000; ll++) {
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
            memcpy(t, t+k_m_l, sizeof(double)*(m - 1));
        } else {
            input.read((char *)t, (m-1)*sizeof(double));
        }
        
        input.read((char *)(t + m -1), (k_m_l)*sizeof(double));
        // if (ll % sample_rate != 0) {
        //     continue;
        // }
        const auto stage1_start = ExperimentClock::now();
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

        compute_normalized_arrays(U, L, U_normal, L_normal, m_k,
                                  miu_max, miu_min, si_max, si_min);

        if(strcmp(tag,"dtw")==0) {
            t_u=find_tu_lp_dtw(U_normal,m_k,q_normal,m,p,LB,m_p,tag,q_sorted);
            t_l=find_tl_lp_dtw(L_normal,m_k,q_normal,m,p,LB,m_p,tag,q_sorted);
        }else if(strcmp(tag,"erp")==0) {
            t_u=find_tu_lp_erp(U_normal,m_k,q_normal,m,p,arg,LB,m_p);
            t_l=find_tl_lp_erp(L_normal,m_k,q_normal,m,p,arg,LB,m_p);
        }else if(strcmp(tag,"swale")==0) {
            t_u=find_tu_lp_swale(U_normal,m_k,q_normal,m,p,arg,LB,m_p,tag,swale_r,swale_p,q_sorted);
            t_l=find_tl_lp_swale(L_normal,m_k,q_normal,m,p,arg,LB,m_p,tag,swale_r,swale_p,q_sorted);
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
        std::vector<FSMELDSurvivor> block_survivors;
        for (int j=0,i=m-1;j<k_m_l;j++,i++) {
            if (invalid_mask[j]) {
                ++invalid_candidates;
                continue;
            }

            kk1++;
            
            if (!fsmeld_lbip_exceeds_cutoff(z[j], evaluation_threshold, metric)) {
                kk2++;
                std::vector<double> t_temp(m);
                for(int tt = 0;tt < m;tt++)
                {
                    t_temp[tt] = (t[tt+j] - miu_t[j]) / si_t[j];
                }
                double bsf = evaluation_threshold;
                // GLB builds the candidate envelope only if its forward bound survives.
                double glb_dist = glb_elastic(metric, q_normal, m, t_temp.data(), m, arg, L_Q, U_Q, w, bsf);

                if (glb_dist <= evaluation_threshold) {
                    kk3++;
                    {
                        block_survivors.push_back({
                            ll * k_m_l + j, std::move(t_temp)});
                    }
                }
            }
        }
        lb_ns_total += std::chrono::duration_cast<std::chrono::nanoseconds>(
            ExperimentClock::now() - stage1_start).count();
        {
            for (FSMELDSurvivor& survivor : block_survivors) {
                const double exact_dist = exact_elastic_distance(
                    metric, q_normal, survivor.normalized.data(), m, w, arg);
                if (std::isfinite(exact_dist) &&
                    exact_dist <= evaluation_threshold) {
                    A.push_back(survivor.index);
                }
            }
        }
    }
        if(!flag_of_tail) {
            if (!read_tail_segment(input, t, siz_t)) {
                std::cerr << "[Error] Failed to read the exact tail segment." << std::endl;
                return;
            }

            const auto tail_stage1_start = ExperimentClock::now();
            auto lp_ip_start = std::chrono::high_resolution_clock::now();
            mvmean(t,siz_t,m,miu_t,si_t);
            memset(invalid_mask, 0, sizeof(bool) * m_k);
            precompute_invalid_mask_fast(t, si_t, siz_t, m, invalid_mask,
                                         tail_first_window, static_cast<int>(tail));
            min_max_window(miu_t,siz_t,m,miu_min,miu_max);
            min_max_window(si_t,siz_t,m,si_min,si_max);

            lower_upper_lemire(t,siz_t,w,L,U);
            compute_normalized_arrays(U, L, U_normal, L_normal, siz_t,
                                      miu_max, miu_min, si_max, si_min);

            if(strcmp(tag,"dtw")==0) {
                t_u=find_tu_lp_dtw(U_normal,siz_t,q_normal,m,p,LB,m_p,tag,q_sorted);
                t_l=find_tl_lp_dtw(L_normal,siz_t,q_normal,m,p,LB,m_p,tag,q_sorted);
            }else if(strcmp(tag,"swale")==0) {
                t_u=find_tu_lp_swale(U_normal,siz_t,q_normal,m,p,arg,LB,m_p,tag,swale_r,swale_p,q_sorted);
                t_l=find_tl_lp_swale(L_normal,siz_t,q_normal,m,p,arg,LB,m_p,tag,swale_r,swale_p,q_sorted);
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
            memset(z, 0, sizeof(double) * k_m_l);

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

            for (int j = tail_first_window, i = m - 1 + tail_first_window;
                 j < tail_end_window; j++, i++) {
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
            std::vector<FSMELDSurvivor> tail_survivors;

            for (int j = tail_first_window, i = m - 1 + tail_first_window;
                 j < tail_end_window; j++, i++) {
                if (invalid_mask[j]) {
                    ++invalid_candidates;
                    continue;
                }
                kk1++;
                
                if (!fsmeld_lbip_exceeds_cutoff(z[j], evaluation_threshold, metric)) {
                    kk2++;

                    std::vector<double> t_temp(m);
                    for(int tt = 0;tt < m;tt++)
                {
                    t_temp[tt] = (t[tt+j] - miu_t[j]) / si_t[j];
                }
                    double bsf = evaluation_threshold;
                // Use the same GLB cascade for tail candidates as for full blocks.
                double glb_dist = glb_elastic(metric, q_normal, m, t_temp.data(), m, arg, L_Q, U_Q, w, bsf);

                if (glb_dist <= evaluation_threshold) {
                    kk3++;
                    {
                        tail_survivors.push_back({
                            forsize * k_m_l + (j - tail_first_window),
                            std::move(t_temp)});
                    }
                }
                }
            }
            lb_ns_total += std::chrono::duration_cast<std::chrono::nanoseconds>(
                ExperimentClock::now() - tail_stage1_start).count();
            {
                for (FSMELDSurvivor& survivor : tail_survivors) {
                    const double exact_dist = exact_elastic_distance(
                        metric, q_normal, survivor.normalized.data(), m, w,
                        arg);
                    if (std::isfinite(exact_dist) &&
                        exact_dist <= evaluation_threshold) {
                        A.push_back(survivor.index);
                    }
                }
            }

        }

        const double total_seconds = std::chrono::duration<double>(
            ExperimentClock::now() - total_start).count();
        const double lb_seconds = lb_ns_total / 1000000000.0;
        const double lbip_pruning_rate = kk1 == 0
            ? 0.0 : static_cast<double>(kk1 - kk2) / kk1;
        const double pruning_rate = kk1 == 0
            ? 0.0 : static_cast<double>(kk1 - kk3) / kk1;

        {
        std::cout << "kk1=" << kk1 << std::endl;
        std::cout << "kk2=" << kk2 << std::endl;
        std::cout << "kk2=" << kk2 << std::endl;
        std::cout << "Pruning Power (guolv)=" << lbip_pruning_rate << std::endl;
        std::cout << "Pruning Power (guolv2)=" << pruning_rate << std::endl;
        std::cout << "LP_IP Lower Bound Time taken: " << (lp_ip_us_total / 1000.0) << " ms" << std::endl;
        std::cout << "LP_IP Lower Bound Time taken: " << (lp_ip_us_total / 1000000.0) << " s" << std::endl;
        std::cout << "T_LB: " << std::setprecision(17) << lb_seconds << " s" << std::endl;
        std::cout << "T_total: " << total_seconds << " s" << std::endl;
        std::cout << "Invalid candidates: " << invalid_candidates << std::endl;
        std::cout << "Final dynamic threshold: " << threshold << std::endl;
        std::cout << "Total Valid Candidates: " << kk1 << std::endl;
        std::cout << "Candidates after LP/IP: " << kk2 << std::endl;
        std::cout << "Cascade: LB_IP+GLB+ELD" << std::endl;
        std::cout << "Candidates after GLB: " << kk3 << std::endl;
        std::cout << "Final Matches: " << A.size() << std::endl;

        const char* quiet_matches_env = std::getenv("FSMELD_QUIET_MATCHES");
        const bool quiet_matches = quiet_matches_env != nullptr &&
            std::string(quiet_matches_env) != "0" &&
            std::string(quiet_matches_env) != "false";
        if (!quiet_matches) {
            for (size_t i = 0; i < A.size(); ++i) {
                printf("A[%zu]=%lld\n", i, A[i]);
            }
        }

        ExperimentResult result;
        result.dataset = infer_dataset_name(bi_t);
        result.normalization = "zscore";
        result.metric = tag;
        result.m = m;
        result.w = w;
        result.xi = threshold;
        result.method = "FSMELD";
        result.lb_time_seconds = lb_seconds;
        result.total_time_seconds = total_seconds;
        result.valid_candidates = kk1;
        result.invalid_candidates = invalid_candidates;
        result.pruned_candidates = kk1 - kk3;
        result.survivors = kk3;
        result.matches = A.size();
        result.has_lbip = true;
        result.lbip_time_seconds = lp_ip_us_total / 1000000.0;
        result.lbip_pruning_rate = lbip_pruning_rate;
        emit_experiment_result(result);
        }
        free(invalid_mask);
        free(t);
        free(q);
        free(q_normal);
        free(z);
        free(z_l);
        free(z_u);
        free(U);
        free(L);
        free(U_normal);
        free(L_normal);
        free(h_l);
        free(h_u);
        free(g_l);
        free(g_u);
        free(miu_t);
        free(si_t);
        free(miu_min);
        free(miu_max);
        free(si_min);
        free(si_max);

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
static void FSMDTW_direct_impl(const char* bi_t, const char* bi_q, int m, double threshold,
                              char *tag, int w, double arg,
                              int exact_top_k) {
    const auto total_start = std::chrono::steady_clock::now();
    vector<long long> A;
    vector<pair<double, long long>> exact_results;
    priority_queue<double> exact_topk_heap;
    
    const bool derive_exact_topk = exact_top_k > 0;
    const char* quiet_matches_env = std::getenv("FSMELD_QUIET_MATCHES");
    const bool quiet_matches = quiet_matches_env != nullptr &&
        std::string(quiet_matches_env) != "0" &&
        std::string(quiet_matches_env) != "false";
    MetricType metric = metric_from_tag(tag);
    int i;
    double d;

    double ex_q = 0, ex2_q = 0, std_val = 0, mean = 0;
    unsigned long long kk1 = 0; //
    unsigned long long kk2 = 0; //
    unsigned long long invalid_candidates = 0;

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

    //
    int m_k;
    if (!resolve_fsmeld_macro_length(m, available_points, &m_k))
        throw std::runtime_error("Target is too short for direct ELD scan.");

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
    for(i = 0; i < m; i++) q_normal[i] = (q[i] - mean) / std_val;

    double *U_Q = (double*) malloc(sizeof(double) * m);
    double *L_Q = (double*) malloc(sizeof(double) * m);

    long long k_m_l = m_k - m + 1;  
    long long ts_n = (sizeofT / (sizeof(double))) - m + 1;  

    std::cout << "ts_n=" << ts_n << std::endl;
    std::cout << "m=" << m << std::endl;
    std::cout << "k_m_l=" << k_m_l << std::endl;
    std::cout << "m_k=" << m_k << std::endl;

    long long forsize = (ts_n) / (k_m_l);
    int flag_of_tail = 0;

    size_t tail = ts_n % k_m_l;
    int siz_t = m_k;
    const int tail_first_window = static_cast<int>(k_m_l - tail);
    const int tail_end_window = tail_first_window + static_cast<int>(tail);

    printf("tail=%zu\n", tail);
    if(tail == 0) flag_of_tail = 1;

    printf("forsize=%lld\n", forsize);
    int forsize_divide_100 = forsize / 100;
    int i_ = 0;

    double *miu_t = (double*) malloc(sizeof(double) * k_m_l);
    double *si_t = (double*) malloc(sizeof(double) * k_m_l);
    bool* invalid_mask = (bool*)malloc(sizeof(bool) * m_k);
    int sample_rate=1000;

    for(long long ll = 0; ll < forsize; ll++) {
        if(forsize_divide_100 && (ll % forsize_divide_100) == 0) {
            std::cout << "\b\b\b\b|" << setw(3) << i_ << "% ";
            i_++;
            
        }

        if(ll != 0) {
            memcpy(t, t + k_m_l, sizeof(double) * (m - 1));
        } else {
            input.read((char *)t, (m - 1) * sizeof(double));
        }
        input.read((char *)(t + m - 1), (k_m_l) * sizeof(double));
//        if (ll % sample_rate != 0) {
//            continue;
//        }

        mvmean(t, m_k, m, miu_t, si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, m_k, m, invalid_mask);

        for (int j = 0; j < k_m_l; j++) {
            if (invalid_mask[j]) {
                ++invalid_candidates;
                continue;
            }
            kk1++;

            double t_temp[m];
            for(int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            bool passed_lower_bound = true;

            if (passed_lower_bound) {
                kk2++;
                const double exact_dist = exact_elastic_distance(
                    metric, q_normal, t_temp, m, w, arg);

                if (!std::isfinite(exact_dist)) {
                    continue;
                }
                record_exact_result(exact_dist, ll * k_m_l + j);
                if (exact_dist <= threshold) {
                    if (!quiet_matches) {
                        printf("\nGET %s: %lf (index: %lld)\n", tag, exact_dist,
                               ll * k_m_l + j);
                    }
                    A.push_back(ll * k_m_l + j);
                }
            }
        }
    }
    
    if(!flag_of_tail) {
        if (!read_tail_segment(input, t, siz_t)) {
            std::cerr << "[Error] Failed to read the exact tail segment." << std::endl;
            return;
        }

        mvmean(t, siz_t, m, miu_t, si_t);
        memset(invalid_mask, 0, sizeof(bool) * m_k);
        precompute_invalid_mask_fast(t, si_t, siz_t, m, invalid_mask,
                                     tail_first_window, static_cast<int>(tail));

        for (int j = tail_first_window; j < tail_end_window; j++) {
            if (invalid_mask[j]) {
                ++invalid_candidates;
                continue;
            }
            kk1++;
            double t_temp[m];
            for(int tt = 0; tt < m; tt++) {
                t_temp[tt] = (t[tt + j] - miu_t[j]) / si_t[j];
            }

            bool passed_lower_bound = true;

            if (passed_lower_bound) {
                kk2++;
                const double exact_dist = exact_elastic_distance(
                    metric, q_normal, t_temp, m, w, arg);

                if (!std::isfinite(exact_dist)) {
                    continue;
                }
                const long long global_index = forsize * k_m_l + (j - tail_first_window);
                record_exact_result(exact_dist, global_index);
                if (exact_dist <= threshold) {
                    if (!quiet_matches) {
                        printf("\nGET %s: %lf (index: %lld)\n", tag, exact_dist,
                               global_index);
                    }
                    A.push_back(global_index);
                }
            }
        }
    }

    if (derive_exact_topk) {
        A.clear();
        if (exact_topk_heap.size() == static_cast<size_t>(exact_top_k)) {
            threshold = exact_topk_heap.top();
            for (const auto& result : exact_results) {
                if (result.first <= threshold) {
                    A.push_back(result.second);
                }
            }
            std::sort(A.begin(), A.end());
            std::cout << std::setprecision(17)
                      << "Exact Top-K Epsilon: " << threshold << std::endl;
        } else {
            std::cout << "Exact Top-K Epsilon: nan" << std::endl;
            std::cerr << "[Error] Only " << exact_topk_heap.size()
                      << " finite ELD values; cannot compute Top-" << exact_top_k << "." << std::endl;
        }
    }

    const double total_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - total_start).count();

    std::cout << "\n\n--- Statistics ---" << std::endl;
    std::cout << "Total comparisons (kk1) = " << kk1 << std::endl;
    {
        std::cout << "Exact ELD evaluations = " << kk2 << std::endl;
    }
    std::cout << "Final matches found: " << A.size() << std::endl;
    std::cout << "Total Valid Candidates: " << kk1 << std::endl;
    std::cout << "Invalid candidates: " << invalid_candidates << std::endl;
    {
        std::cout << "Exact ELD Evaluations: " << kk2 << std::endl;
    }
    std::cout << "Final Matches: " << A.size() << std::endl;

    if (!quiet_matches) {
        for (size_t k = 0; k < A.size(); ++k) {
            printf("A[%zu] = %lld\n", k, A[k]);
        }
    }

    {
        std::cout << "T_total: " << std::setprecision(17)
                  << total_seconds << " s" << std::endl;

        ExperimentResult result;
        result.dataset = infer_dataset_name(bi_t);
        result.normalization = "zscore";
        result.metric = tag;
        result.m = m;
        result.w = w;
        result.xi = threshold;
        result.method = "ELD";
        result.lb_time_seconds = 0.0;
        result.total_time_seconds = total_seconds;
        result.valid_candidates = kk1;
        result.invalid_candidates = invalid_candidates;
        result.pruned_candidates = 0;
        result.survivors = kk1;
        result.matches = A.size();
        emit_experiment_result(result);
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

void FSMELD(const char* bi_t, const char* bi_q, int m, double threshold,
            char* tag, int w, double arg, long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    FSMDTW_impl(bi_t, bi_q, m, threshold, tag, w, arg);
}
void FSMELD_find_threshold(const char* bi_t, const char* bi_q, int m, double threshold,
                          char* tag, int w, double arg, int top_k,
                          long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    find_threshold_impl(bi_t, bi_q, m, threshold, tag, w, arg, top_k, false);
}
void FSMELD_exact(const char* bi_t, const char* bi_q, int m, double threshold,
                  char* tag, int w, double arg, int top_k,
                  long long exclude_start, long long exclude_end) {
    ScopedCandidateExclusion exclusion(exclude_start, exclude_end);
    FSMDTW_direct_impl(bi_t, bi_q, m, threshold, tag, w, arg, top_k);
}
