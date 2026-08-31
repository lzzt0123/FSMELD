#ifndef FSMELD_FUNCS_H
#define FSMELD_FUNCS_H

#include "fftw3.h"

enum MetricType {
    METRIC_UNKNOWN = 0,
    METRIC_DTW,
    METRIC_EDR,
    METRIC_ERP,
    METRIC_LCSS,
    METRIC_MSM
};

MetricType metric_from_tag(const char* tag);

void inv_and_padding_y(int n, int m, double* y, fftw_complex* YY, fftw_plan p);
void FFT_x(double* x, fftw_complex* XX, fftw_plan p);
void cov_IFFT(fftw_complex* XX, fftw_complex* YY, fftw_complex* ZZ, int n, double* z, fftw_plan p);
void SIP(double* h, int m, double* g, int n, double* z);

void compute_normalized_arrays(const double* U, const double* L, double* U_normal, double* L_normal,
                               int length, double* max_mean, double* min_mean,
                               double* max_std, double* min_std);
void mvmean(double* a, int len_a, int l, double* miu, double* si);

void q_to_h_u(const double* q, int q_size, double t_u, double C, double* h_u, int h_size, const char* tag);
void q_to_h_l(const double* q, int q_size, double t_l, double C, double* h_l, int h_size, const char* tag);
void s_to_g_l(const double* L, int L_size, double t_l, double* g_l, int g_size);
void s_to_g_u(const double* U, int U_size, double t_u, double* g_u, int g_size);
void s_to_g_l(const double* L, int L_size, double t_l, double* g_l, int g_size, const char* tag, double C);
void s_to_g_u(const double* U, int U_size, double t_u, double* g_u, int g_size, const char* tag, double C);

void lower_upper_lemire(double a[], int n, int r, double l[], double u[]);
void min_max_window(double a[], int n, int r, double min_out[], double max_out[]);

double find_tl_lp(double* L, int l, double* q, int m, int p, double C_, double* LB, double* m_p,
                  const char* tag, const double* q_sorted = nullptr);
double find_tl_lp_dtw(double* L, int l, double* q, int m, int p, double* LB, double* m_p,
                      const char* tag, const double* q_sorted = nullptr);
double find_tl_lp_erp(double* L, int l, double* q, int m, int p, double C_, double* LB, double* m_p);
double find_tu_lp(double* U, int l, double* q, int m, int p, double C_, double* LB, double* m_p,
                  const char* tag, const double* q_sorted = nullptr);
double find_tu_lp_dtw(double* U, int l, double* q, int m, int p, double* LB, double* m_p,
                      const char* tag, const double* q_sorted = nullptr);
double find_tu_lp_erp(double* U, int l, double* q, int m, int p, double C_, double* LB, double* m_p);

void FSMELD(const char* bi_t, const char* bi_q, int m, double threshold, char* tag, int w, double arg,
            long long exclude_start = -1, long long exclude_end = -1);
void FSMELD_raw(const char* bi_t, const char* bi_q, int m, double threshold, char* tag, int w, double arg,
                long long exclude_start = -1, long long exclude_end = -1);
void FSMELD_exact(const char* bi_t, const char* bi_q, int m, double threshold, char* tag, int w,
                  double arg, int top_k, long long exclude_start = -1, long long exclude_end = -1);
void FSMELD_epsilon(const char* bi_t, const char* bi_q, int m, double threshold, char* tag, int w,
                    double arg, int top_k);
void FSMELD_epsilon(const char* bi_t, const char* bi_q, int m, double threshold, char* tag, int w,
                    double arg, int top_k, long long exclude_start, long long exclude_end);
void FSMELD_find_threshold(const char* bi_t, const char* bi_q, int m, double threshold, char* tag, int w,
                           double arg, int top_k, long long exclude_start = -1, long long exclude_end = -1);
void FSMELD_get_tlb(const char* bi_t, const char* bi_q, int m, double threshold, char* tag, int w,
                    double arg, double target_index, long long exclude_start = -1, long long exclude_end = -1);

double dtw(double* A, int m, double* B, int n, int r);
double lcss(double* Q, int n, double* A, int m, int w, double epsilon);
double erp(double* x, int lenx, double* y, int leny, double g, int w);
double edr(double* x, int n, double* y, int m, double epsilon, int w);
double msm(double* x, int xlen, double* y, int ylen, double c, int w);

double bglb_dtw(double* Q, int lenQ, double* T, int lenT, double* L_Q, double* U_Q,
                double* L_T, double* U_T, int w, double& bsf);
double bglb_edr(double* Q, int lenQ, double* T, int lenT, double epsilon, double* U_Q, double* L_Q,
                double* U_T, double* L_T, int w, double& bsf);
double bglb_lcss(double* Q, int lenQ, double* T, int lenT, double epsilon, double* U_Q, double* L_Q,
                 double* U_T, double* L_T, int w, double& bsf);
double bglb_erp(double* Q, int lenQ, double* T, int lenT, double g, double* U_Q, double* L_Q,
                double* U_T, double* L_T, int w, double& bsf);
double bglb_msm(double* Q, int lenQ, double* T, int lenT, double c, double* U_Q, double* L_Q,
                double* U_T, double* L_T, int w, double& bsf);
double bglb_elastic_lazy(MetricType metric, double* Q, int lenQ, double* T, int lenT, double arg,
                         double* L_Q, double* U_Q, int w, double& bsf);

#endif
