#include <algorithm>
#include <cmath>
#include <cstring>
#include "def.h"
#include "funcs.h"
using namespace std;

double edr_subcost(double a, double b, double epsilon) {
    return (abs(a - b) <= epsilon) ? 0.0 : 1.0;
}
static double swale_subcost(double a, double b, double epsilon,
                            double p, double r_cost) {
    return (abs(a - b) <= epsilon) ? r_cost : p;
}

// Match exact LCSS distance comparisons at inclusive tolerance boundaries.
static inline bool glb_envelope_match(double value, double lower, double upper,
                                      double epsilon) {
    if (value < lower) return lower - value <= epsilon;
    if (value > upper) return value - upper <= epsilon;
    return true;
}

static inline bool glb_exceeds_bsf(double partial_cost, double bsf) {
    return partial_cost > bsf;
}

static inline bool glb_exceeds_squared_bsf(double partial_cost, double bsf) {
    const double cutoff = bsf;
    return partial_cost > cutoff * cutoff;
}
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------

// GLB implementations (see THIRD_PARTY_NOTICES.txt), adapted to the
// distance scales used by this project. DTW and ERP return sqrt(cost) here,
// while the reference repository keeps their squared costs.
double glb_dtw(double* Q, int lenQ, double* T, int lenT,
               double* L_Q, double* U_Q, double* L_T, double* U_T,
               double& bsf) {
    const double boundary = dist(T[0], Q[0]) + dist(T[lenT - 1], Q[lenQ - 1]);
    double t_to_q = 0.0;
    for (int i = 1; i < lenT - 1; ++i) {
        if (T[i] < L_Q[i]) t_to_q += dist(T[i], L_Q[i]);
        else if (T[i] > U_Q[i]) t_to_q += dist(T[i], U_Q[i]);
        if (glb_exceeds_squared_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
    }

    double q_to_t = 0.0;
    for (int i = 1; i < lenQ - 1; ++i) {
        if (Q[i] < L_T[i]) q_to_t += dist(Q[i], L_T[i]);
        else if (Q[i] > U_T[i]) q_to_t += dist(Q[i], U_T[i]);
        if (glb_exceeds_squared_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
    }
    return sqrt(max(0.0, boundary + MAX(t_to_q, q_to_t)));
}

double glb_lcss(double* Q, int lenQ, double* T, int lenT, double epsilon,
                 double* L_Q, double* U_Q, double* L_T, double* U_T) {
    // ElasticMeasures::lcss returns the loss 1 - LCSS/min(lenQ,lenT) and
    // does not force either endpoint into the common subsequence.  Hence an
    // admissible similarity upper bound must count every position that could
    // match inside the opposite envelope; adding fixed endpoint rewards would
    // double-count/force endpoints and can make the resulting loss exceed the
    // exact LCSS distance.
    double t_reward = 0.0;
    for (int i = 0; i < lenT; ++i) {
        if (glb_envelope_match(T[i], L_Q[i], U_Q[i], epsilon)) {
            t_reward += 1.0;
        }
    }

    double q_reward = 0.0;
    for (int i = 0; i < lenQ; ++i) {
        if (glb_envelope_match(Q[i], L_T[i], U_T[i], epsilon)) {
            q_reward += 1.0;
        }
    }

    const double min_len = MIN(lenQ, lenT);
    const double similarity_upper = MIN(t_reward, q_reward);
    return max(0.0, 1.0 - similarity_upper / min_len);
}

double glb_erp(double* Q, int lenQ, double* T, int lenT, double g_val,
               double* L_Q, double* U_Q, double* L_T, double* U_T,
               double& bsf) {
    const double boundary = min({
        dist(Q[lenQ - 1], T[lenT - 1]),
        dist(Q[lenQ - 1], g_val),
        dist(T[lenT - 1], g_val)
    });

    double t_to_q = 0.0;
    for (int i = 1; i < lenT - 1; ++i) {
        if (T[i] > U_Q[i]) {
            t_to_q += MIN(dist(T[i], U_Q[i]), dist(T[i], g_val));
        } else if (T[i] < L_Q[i]) {
            t_to_q += MIN(dist(T[i], L_Q[i]), dist(T[i], g_val));
        }
        if (glb_exceeds_squared_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
    }

    double q_to_t = 0.0;
    for (int i = 1; i < lenQ - 1; ++i) {
        if (Q[i] > U_T[i]) {
            q_to_t += MIN(dist(Q[i], U_T[i]), dist(Q[i], g_val));
        } else if (Q[i] < L_T[i]) {
            q_to_t += MIN(dist(Q[i], L_T[i]), dist(Q[i], g_val));
        }
        if (glb_exceeds_squared_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
    }
    return sqrt(max(0.0, boundary + MAX(t_to_q, q_to_t)));
}

double glb_edr(double* Q, int lenQ, double* T, int lenT, double epsilon,
               double* L_Q, double* U_Q, double* L_T, double* U_T,
               double& bsf) {
    const double boundary = MIN(
        edr_subcost(Q[lenQ - 1], T[lenT - 1], epsilon), 1.0);
    double t_to_q = 0.0;
    for (int i = 1; i < lenT - 1; ++i) {
        if (T[i] > U_Q[i]) t_to_q += edr_subcost(T[i], U_Q[i], epsilon);
        else if (T[i] < L_Q[i]) t_to_q += edr_subcost(T[i], L_Q[i], epsilon);
        if (glb_exceeds_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
    }

    double q_to_t = 0.0;
    for (int i = 1; i < lenQ - 1; ++i) {
        if (Q[i] > U_T[i]) q_to_t += edr_subcost(Q[i], U_T[i], epsilon);
        else if (Q[i] < L_T[i]) q_to_t += edr_subcost(Q[i], L_T[i], epsilon);
        if (glb_exceeds_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
    }
    return boundary + MAX(t_to_q, q_to_t);
}

double glb_swale(double* Q, int lenQ, double* T, int lenT, double epsilon,
                 double p, double r_cost, double* L_Q, double* U_Q,
                 double* L_T, double* U_T, double& bsf) {
    const double boundary = MIN(
        swale_subcost(Q[lenQ - 1], T[lenT - 1], epsilon, p, r_cost), p);
    double t_to_q = 0.0;
    for (int i = 1; i < lenT - 1; ++i) {
        if (T[i] > U_Q[i]) {
            t_to_q += swale_subcost(T[i], U_Q[i], epsilon, p, r_cost);
        } else if (T[i] < L_Q[i]) {
            t_to_q += swale_subcost(T[i], L_Q[i], epsilon, p, r_cost);
        }
        if (glb_exceeds_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
    }

    double q_to_t = 0.0;
    for (int i = 1; i < lenQ - 1; ++i) {
        if (Q[i] > U_T[i]) {
            q_to_t += swale_subcost(Q[i], U_T[i], epsilon, p, r_cost);
        } else if (Q[i] < L_T[i]) {
            q_to_t += swale_subcost(Q[i], L_T[i], epsilon, p, r_cost);
        }
        if (glb_exceeds_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
    }
    return boundary + MAX(t_to_q, q_to_t);
}

double glb_msm(double* Q, int lenQ, double* T, int lenT, double c,
               double* L_Q, double* U_Q, double* L_T, double* U_T,
               double& bsf) {
    double boundary;
    // The source implementation uses chained comparisons, which do not express the
    // intended interval tests in C++. Expand them explicitly.
    if ((T[lenT - 2] >= T[lenT - 1] && T[lenT - 1] >= Q[lenQ - 1]) ||
        (T[lenT - 2] <= T[lenT - 1] && T[lenT - 1] <= Q[lenQ - 1]) ||
        (Q[lenQ - 2] <= Q[lenQ - 1] && Q[lenQ - 1] <= T[lenT - 1]) ||
        (Q[lenQ - 2] >= Q[lenQ - 1] && Q[lenQ - 1] >= T[lenT - 1])) {
        boundary = abs(Q[0] - T[0]) + MIN(abs(Q[lenQ - 1] - T[lenT - 1]), c);
    } else {
        boundary = abs(Q[0] - T[0]) + min({
            abs(Q[lenQ - 1] - T[lenT - 1]),
            c + abs(T[lenT - 1] - T[lenT - 2]),
            c + abs(Q[lenQ - 1] - Q[lenQ - 2])
        });
    }

    double t_to_q = 0.0;
    for (int i = 1; i < lenT - 1; ++i) {
        if (T[i] > U_Q[i]) t_to_q += MIN(abs(T[i] - U_Q[i]), c);
        else if (T[i] < L_Q[i]) t_to_q += MIN(abs(T[i] - L_Q[i]), c);
        if (glb_exceeds_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
    }

    double q_to_t = 0.0;
    for (int i = 1; i < lenQ - 1; ++i) {
        if (Q[i] > U_T[i]) q_to_t += MIN(abs(Q[i] - U_T[i]), c);
        else if (Q[i] < L_T[i]) q_to_t += MIN(abs(Q[i] - L_T[i]), c);
        if (glb_exceeds_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
    }
    return boundary + MAX(t_to_q, q_to_t);
}

double glb_elastic(MetricType metric, double* Q, int lenQ, double* T, int lenT,
                    double arg, double* L_Q, double* U_Q, int w, double& bsf) {
    // Endpoint formulas need two positions. Defer one-point cases to ELD;
    // zero remains admissible and avoids duplicate endpoints/MSM index -1.
    if (lenQ < 2 || lenT < 2) return 0.0;
    // Lazy candidate-envelope construction.  The query envelope is reusable
    // and has already been computed by the caller.  First evaluate the full
    // T->Q directional GLB; only candidates surviving that direction pay for
    // construction of their own envelope.  The directional contribution is
    // retained, so unlike a precheck-plus-full-call implementation this pass
    // is not recomputed.
    double boundary = 0.0;
    double t_to_q = 0.0;
    double lcss_t_reward = 0.0;

    switch (metric) {
        case METRIC_DTW: {
            boundary = dist(T[0], Q[0]) + dist(T[lenT - 1], Q[lenQ - 1]);
            for (int i = 1; i < lenT - 1; ++i) {
                if (T[i] < L_Q[i]) t_to_q += dist(T[i], L_Q[i]);
                else if (T[i] > U_Q[i]) t_to_q += dist(T[i], U_Q[i]);
                if (glb_exceeds_squared_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
            }
            break;
        }
        case METRIC_ERP: {
            boundary = min({
                dist(Q[lenQ - 1], T[lenT - 1]),
                dist(Q[lenQ - 1], arg),
                dist(T[lenT - 1], arg)
            });
            for (int i = 1; i < lenT - 1; ++i) {
                if (T[i] > U_Q[i]) {
                    t_to_q += MIN(dist(T[i], U_Q[i]), dist(T[i], arg));
                } else if (T[i] < L_Q[i]) {
                    t_to_q += MIN(dist(T[i], L_Q[i]), dist(T[i], arg));
                }
                if (glb_exceeds_squared_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
            }
            break;
        }
        case METRIC_EDR:
            boundary = MIN(edr_subcost(Q[lenQ - 1], T[lenT - 1], arg), 1.0);
            for (int i = 1; i < lenT - 1; ++i) {
                if (T[i] > U_Q[i]) t_to_q += edr_subcost(T[i], U_Q[i], arg);
                else if (T[i] < L_Q[i]) t_to_q += edr_subcost(T[i], L_Q[i], arg);
                if (glb_exceeds_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
            }
            break;
        case METRIC_LCSS: {
            const double min_len = MIN(lenQ, lenT);
            for (int i = 0; i < lenT; ++i) {
                if (glb_envelope_match(T[i], L_Q[i], U_Q[i], arg)) {
                    lcss_t_reward += 1.0;
                }
            }
            // The second direction can only reduce the min(reward1,reward2),
            // so this one-direction loss is already a valid early-stop bound.
            if (glb_exceeds_bsf(max(0.0, 1.0 - lcss_t_reward / min_len), bsf)) {
                return MAX_VAL;
            }
            break;
        }
        case METRIC_MSM:
            if ((T[lenT - 2] >= T[lenT - 1] && T[lenT - 1] >= Q[lenQ - 1]) ||
                (T[lenT - 2] <= T[lenT - 1] && T[lenT - 1] <= Q[lenQ - 1]) ||
                (Q[lenQ - 2] <= Q[lenQ - 1] && Q[lenQ - 1] <= T[lenT - 1]) ||
                (Q[lenQ - 2] >= Q[lenQ - 1] && Q[lenQ - 1] >= T[lenT - 1])) {
                boundary = abs(Q[0] - T[0]) + MIN(abs(Q[lenQ - 1] - T[lenT - 1]), arg);
            } else {
                boundary = abs(Q[0] - T[0]) + min({
                    abs(Q[lenQ - 1] - T[lenT - 1]),
                    arg + abs(T[lenT - 1] - T[lenT - 2]),
                    arg + abs(Q[lenQ - 1] - Q[lenQ - 2])
                });
            }
            for (int i = 1; i < lenT - 1; ++i) {
                if (T[i] > U_Q[i]) t_to_q += MIN(abs(T[i] - U_Q[i]), arg);
                else if (T[i] < L_Q[i]) t_to_q += MIN(abs(T[i] - L_Q[i]), arg);
                if (glb_exceeds_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
            }
            break;
        case METRIC_SWALE:
            boundary = MIN(
                swale_subcost(Q[lenQ - 1], T[lenT - 1], arg,
                              SWALE_MISMATCH_COST, SWALE_MATCH_COST),
                SWALE_MISMATCH_COST);
            for (int i = 1; i < lenT - 1; ++i) {
                if (T[i] > U_Q[i]) {
                    t_to_q += swale_subcost(T[i], U_Q[i], arg,
                                            SWALE_MISMATCH_COST, SWALE_MATCH_COST);
                } else if (T[i] < L_Q[i]) {
                    t_to_q += swale_subcost(T[i], L_Q[i], arg,
                                            SWALE_MISMATCH_COST, SWALE_MATCH_COST);
                }
                if (glb_exceeds_bsf(boundary + t_to_q, bsf)) return MAX_VAL;
            }
            break;
        default:
            return MAX_VAL;
    }

    double L_T[lenT];
    double U_T[lenT];
    lower_upper_lemire(T, lenT, w, L_T, U_T);

    switch (metric) {
        case METRIC_DTW: {
            double q_to_t = 0.0;
            for (int i = 1; i < lenQ - 1; ++i) {
                if (Q[i] < L_T[i]) q_to_t += dist(Q[i], L_T[i]);
                else if (Q[i] > U_T[i]) q_to_t += dist(Q[i], U_T[i]);
                if (glb_exceeds_squared_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
            }
            return sqrt(max(0.0, boundary + MAX(t_to_q, q_to_t)));
        }
        case METRIC_ERP: {
            double q_to_t = 0.0;
            for (int i = 1; i < lenQ - 1; ++i) {
                if (Q[i] > U_T[i]) {
                    q_to_t += MIN(dist(Q[i], U_T[i]), dist(Q[i], arg));
                } else if (Q[i] < L_T[i]) {
                    q_to_t += MIN(dist(Q[i], L_T[i]), dist(Q[i], arg));
                }
                if (glb_exceeds_squared_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
            }
            return sqrt(max(0.0, boundary + MAX(t_to_q, q_to_t)));
        }
        case METRIC_EDR: {
            double q_to_t = 0.0;
            for (int i = 1; i < lenQ - 1; ++i) {
                if (Q[i] > U_T[i]) q_to_t += edr_subcost(Q[i], U_T[i], arg);
                else if (Q[i] < L_T[i]) q_to_t += edr_subcost(Q[i], L_T[i], arg);
                if (glb_exceeds_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
            }
            return boundary + MAX(t_to_q, q_to_t);
        }
        case METRIC_LCSS: {
            double q_reward = 0.0;
            for (int i = 0; i < lenQ; ++i) {
                if (glb_envelope_match(Q[i], L_T[i], U_T[i], arg)) {
                    q_reward += 1.0;
                }
            }
            const double min_len = MIN(lenQ, lenT);
            return max(0.0, 1.0 - MIN(lcss_t_reward, q_reward) / min_len);
        }
        case METRIC_MSM: {
            double q_to_t = 0.0;
            for (int i = 1; i < lenQ - 1; ++i) {
                if (Q[i] > U_T[i]) q_to_t += MIN(abs(Q[i] - U_T[i]), arg);
                else if (Q[i] < L_T[i]) q_to_t += MIN(abs(Q[i] - L_T[i]), arg);
                if (glb_exceeds_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
            }
            return boundary + MAX(t_to_q, q_to_t);
        }
        case METRIC_SWALE: {
            double q_to_t = 0.0;
            for (int i = 1; i < lenQ - 1; ++i) {
                if (Q[i] > U_T[i]) {
                    q_to_t += swale_subcost(Q[i], U_T[i], arg,
                                            SWALE_MISMATCH_COST, SWALE_MATCH_COST);
                } else if (Q[i] < L_T[i]) {
                    q_to_t += swale_subcost(Q[i], L_T[i], arg,
                                            SWALE_MISMATCH_COST, SWALE_MATCH_COST);
                }
                if (glb_exceeds_bsf(boundary + q_to_t, bsf)) return MAX_VAL;
            }
            return boundary + MAX(t_to_q, q_to_t);
        }
        default:
            return MAX_VAL;
    }
}

// Metric name conversion.

MetricType metric_from_tag(const char* tag) {
    if (strcmp("dtw", tag) == 0) return METRIC_DTW;
    if (strcmp("edr", tag) == 0) return METRIC_EDR;
    if (strcmp("erp", tag) == 0) return METRIC_ERP;
    if (strcmp("lcss", tag) == 0) return METRIC_LCSS;
    if (strcmp("msm", tag) == 0) return METRIC_MSM;
    if (strcmp("swale", tag) == 0) return METRIC_SWALE;
    return METRIC_UNKNOWN;
}
