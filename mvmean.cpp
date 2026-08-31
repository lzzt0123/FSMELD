#include <cstdio>
#include "funcs.h"
#include "cmath"
using namespace std;
#include <fstream>


void mvmean(double *a, int len_a, int l, double *miu, double *si)
{
    
    double sum1 = 0;
    double sum2 = 0;
    for(int i = 0; i < l; i++)
    {
        sum1 += a[i];
        sum2 += a[i] * a[i];
    }

    int ll = len_a - l + 1;

    
    miu[0] = sum1 / l;
    si[0] = sqrt(sum2 / l - miu[0] * miu[0]);

    
    for(int i = 1; i < ll; i++)
    {
        sum1 = sum1 - a[i-1] + a[i+l-1];
        sum2 = sum2 - a[i-1] * a[i-1] + a[i+l-1] * a[i+l-1];
        
        miu[i] = sum1 / l;
        si[i] = sqrt(sum2 / l - miu[i] * miu[i]);

    }
}

