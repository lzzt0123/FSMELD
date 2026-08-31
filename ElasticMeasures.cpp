#include <algorithm>
#include <cmath>
#include <cstddef>
#include <math.h>

#define MAX_VAL 1e18
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

inline double dist(double a, double b) {
    return (a - b) * (a - b);
}

static double* ensure_buffer(double*& buffer, size_t& capacity, size_t needed) {
    if (capacity < needed) {
        delete[] buffer;
        buffer = new double[needed];
        capacity = needed;
    }
    return buffer;
}

double dtw(double* A, int m, double* B, int n, int r) {
    int row_len = 2 * r + 1;
    static double* cost_buffer = nullptr;
    static double* cost_prev_buffer = nullptr;
    static size_t cost_capacity = 0;
    static size_t cost_prev_capacity = 0;
    double* cost = ensure_buffer(cost_buffer, cost_capacity, row_len);
    double* cost_prev = ensure_buffer(cost_prev_buffer, cost_prev_capacity, row_len);

    for (int i = 0; i < row_len; i++) {
        cost[i] = MAX_VAL;
        cost_prev[i] = MAX_VAL;
    }

    for (int i = 0; i < m; i++) {
        int k = MAX(0, r - i);
        for (int j = std::max(0, i - r); j <= std::min(n - 1, i + r); j++, k++) {
            if ((i == 0) && (j == 0)) {
                cost[k] = dist(A[0], B[0]);
                continue;
            }

            double x, y, z;
            y = ((j - 1 < 0) || (k - 1 < 0)) ? MAX_VAL : cost[k - 1];
            x = ((i - 1 < 0) || (k + 1 > 2 * r)) ? MAX_VAL : cost_prev[k + 1];
            z = ((i - 1 < 0) || (j - 1 < 0)) ? MAX_VAL : cost_prev[k];

            cost[k] = std::min({x, y, z}) + dist(A[i], B[j]);
        }
        double* tmp = cost;
        cost = cost_prev;
        cost_prev = tmp;
    }

    double final_dtw = cost_prev[r + (n - 1) - (m - 1)];
    if (!std::isfinite(final_dtw)) return MAX_VAL;
    return std::sqrt(std::max(0.0, final_dtw));
}

double lcss(double* Q, int n, double* A, int m, int w, double epsilon) {
    static double* arr_buffer = nullptr;
    static size_t arr_capacity = 0;
    double* arr = ensure_buffer(arr_buffer, arr_capacity, static_cast<size_t>(n) * m);
    for (int i = 0; i < n * m; i++) arr[i] = 0.0;

    for (int i = 0; i < n; ++i) {
        int wmin = std::max(0, i - w);
        int wmax = std::min(m, i + w + 1);
        for (int j = wmin; j < wmax; ++j) {
            int idx = i * m + j;
            if (i + j == 0) {
                if (std::abs(Q[i] - A[j]) <= epsilon) arr[idx] = 1;
            } else if (i == 0) {
                if (std::abs(Q[i] - A[j]) <= epsilon) arr[idx] = 1;
                else arr[idx] = (j > wmin) ? arr[i * m + (j - 1)] : 0;
            } else if (j == 0) {
                if (std::abs(Q[i] - A[j]) <= epsilon) arr[idx] = 1;
                else arr[idx] = arr[(i - 1) * m + j];
            } else {
                double val1 = arr[(i - 1) * m + (j - 1)];
                double val2 = arr[(i - 1) * m + j];
                double val3 = (j > wmin) ? arr[idx - 1] : 0;
                if (std::abs(Q[i] - A[j]) <= epsilon) arr[idx] = std::max({val1 + 1, val2, val3});
                else arr[idx] = std::max({val1, val2, val3});
            }
        }
    }
    double result = arr[(n - 1) * m + (m - 1)];
    return 1 - result / MIN(n, m);
}

double erp(double* x, int lenx, double* y, int leny, double g, int w) {
    static double* cur_buffer = nullptr;
    static double* prev_buffer = nullptr;
    static size_t cur_capacity = 0;
    static size_t prev_capacity = 0;
    const size_t row_len = static_cast<size_t>(leny) + 1;
    double* cur = ensure_buffer(cur_buffer, cur_capacity, row_len);
    double* prev = ensure_buffer(prev_buffer, prev_capacity, row_len);

    std::fill(prev, prev + row_len, MAX_VAL);
    prev[0] = 0.0;
    for (int j = 1; j <= std::min(leny, w); ++j) {
        prev[j] = prev[j - 1] + dist(y[j - 1], g);
    }

    for (int i = 1; i <= lenx; ++i) {
        std::fill(cur, cur + row_len, MAX_VAL);
        if (i <= w) {
            cur[0] = prev[0] + dist(x[i - 1], g);
        }

        const int minw = std::max(1, i - w);
        const int maxw = std::min(leny, i + w);
        for (int j = minw; j <= maxw; ++j) {
            const double match = prev[j - 1] + dist(x[i - 1], y[j - 1]);
            const double insert = cur[j - 1] + dist(y[j - 1], g);
            const double delete_c = prev[j] + dist(x[i - 1], g);
            cur[j] = std::min({match, insert, delete_c});
        }
        std::swap(cur, prev);
    }

    double res = prev[leny];
    if (!std::isfinite(res)) return MAX_VAL;
    return std::sqrt(std::max(0.0, res));
}

double edr(double* x, int n, double* y, int m_size, double epsilon, int w) {
    static double* cur_buffer = nullptr;
    static double* prev_buffer = nullptr;
    static size_t cur_capacity = 0;
    static size_t prev_capacity = 0;
    const size_t row_len = static_cast<size_t>(m_size) + 1;
    double* cur = ensure_buffer(cur_buffer, cur_capacity, row_len);
    double* prev = ensure_buffer(prev_buffer, prev_capacity, row_len);

    std::fill(prev, prev + row_len, MAX_VAL);
    prev[0] = 0.0;
    for (int j = 1; j <= std::min(m_size, w); ++j) {
        prev[j] = static_cast<double>(j);
    }

    for (int i = 1; i <= n; ++i) {
        std::fill(cur, cur + row_len, MAX_VAL);
        if (i <= w) {
            cur[0] = static_cast<double>(i);
        }

        const int minw = std::max(1, i - w);
        const int maxw = std::min(m_size, i + w);
        for (int j = minw; j <= maxw; ++j) {
            const double substitution =
                (std::abs(x[i - 1] - y[j - 1]) <= epsilon) ? 0.0 : 1.0;
            cur[j] = std::min({prev[j] + 1.0,
                               cur[j - 1] + 1.0,
                               prev[j - 1] + substitution});
        }
        std::swap(cur, prev);
    }
    return prev[m_size];
}


double msm(double* x, int xlen, double* y, int ylen, double c, int w) {
    static double* cost_buffer = nullptr;
    static size_t cost_capacity = 0;
    double* cost = ensure_buffer(cost_buffer, cost_capacity, static_cast<size_t>(xlen) * ylen);
    for (int i = 0; i < xlen * ylen; i++) cost[i] = MAX_VAL;

    auto msm_dist_local = [](double np, double x, double y, double c) {
        if ((x <= np && np <= y) || (y <= np && np <= x)) return c;
        return c + std::min(std::abs(np - x), std::abs(np - y));
    };

    cost[0] = std::abs(x[0] - y[0]);
    for (int i = 1; i < xlen && i <= w; ++i) {
        cost[i * ylen] = cost[(i - 1) * ylen] + msm_dist_local(x[i], x[i - 1], y[0], c);
    }
    for (int j = 1; j < ylen && j <= w; ++j) {
        cost[j] = cost[j - 1] + msm_dist_local(y[j], x[0], y[j - 1], c);
    }

    for (int i = 1; i < xlen; ++i) {
        int j_start = MAX(1, i - w);
        int j_end = MIN(ylen, i + w + 1);
        for (int j = j_start; j < j_end; ++j) {
            double c1 = cost[(i - 1) * ylen + (j - 1)] + std::abs(x[i] - y[j]);
            double c2 = cost[(i - 1) * ylen + j] + msm_dist_local(x[i], x[i - 1], y[j], c);
            double c3 = cost[i * ylen + (j - 1)] + msm_dist_local(y[j], x[i], y[j - 1], c);
            cost[i * ylen + j] = std::min({c1, c2, c3});
        }
    }
    double res = cost[xlen * ylen - 1];
    return res;
}

