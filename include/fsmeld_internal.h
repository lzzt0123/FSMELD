#ifndef FSMELD_INTERNAL_H
#define FSMELD_INTERNAL_H

#include "funcs.h"

#include <fstream>
#include <cmath>

inline bool fsmeld_lbip_exceeds_cutoff(double bound, double cutoff, MetricType metric) {
    constexpr double tolerance = 1e-9;
    // FFT roundoff is in squared-cost units for DTW/ERP. Near zero, taking
    // sqrt first amplifies it; uncertain candidates must reach GLB/exact ELD.
    if (!std::isfinite(bound) || !(bound > cutoff + tolerance)) return false;
    if (metric == METRIC_DTW || metric == METRIC_ERP) {
        const double limit = cutoff + tolerance;
        return bound > std::sqrt(limit * limit + tolerance);
    }
    return true;
}

void trans(double z_l[], double z_u[], double z[], int m, int window_count, char* tag);
void precompute_invalid_mask_fast(const double* t, const double* si_t,
                                  int segment_size, int m, bool* invalid_mask,
                                  int first_processed_window = 0,
                                  int processed_window_count = -1);
bool read_tail_segment(std::ifstream& input, double* t, int segment_size);
double exact_elastic_distance(MetricType metric, double* q, double* t,
                              int m, int w, double arg);

#endif
