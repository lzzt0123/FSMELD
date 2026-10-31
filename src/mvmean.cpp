#include <cstdio>
#include "funcs.h"
#include "cmath"
#include <limits>
using namespace std;
#include <fstream>
/// calculates the moving average (mean) and standard deviation of a sliding window over an array
/// also returns the maximum and minimum values of the moving averages and standard deviations
void mvmean(double *a, int len_a, int l, double *miu, double *si)
{
    
    double sum1 = 0;
    double sum2 = 0;
    int nonfinite_count = 0;
    const double invalid = std::numeric_limits<double>::quiet_NaN();
    for(int i = 0; i < l; i++)
    {
        if (std::isfinite(a[i])) {
            sum1 += a[i];
            sum2 += a[i] * a[i];
        } else {
            ++nonfinite_count;
        }
    }

    int ll = len_a - l + 1;

    miu[0] = nonfinite_count ? invalid : sum1 / l;
    si[0] = nonfinite_count ? invalid : sqrt(sum2 / l - miu[0] * miu[0]);

    for(int i = 1; i < ll; i++)
    {
        // Keep nonfinite values out of the sums so later finite windows recover.
        const bool outgoing_finite = std::isfinite(a[i-1]);
        const bool incoming_finite = std::isfinite(a[i+l-1]);
        const double outgoing = outgoing_finite ? a[i-1] : 0.0;
        const double incoming = incoming_finite ? a[i+l-1] : 0.0;
        nonfinite_count += !incoming_finite - !outgoing_finite;
        sum1 = sum1 - outgoing + incoming;
        sum2 = sum2 - outgoing * outgoing + incoming * incoming;
        
        miu[i] = nonfinite_count ? invalid : sum1 / l;
        si[i] = nonfinite_count ? invalid : sqrt(sum2 / l - miu[i] * miu[i]);

    }
}
