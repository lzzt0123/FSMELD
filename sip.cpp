



#include "fftw3.h"
#include "funcs.h"
#include<stdio.h>
#include "fstream"
#include<iostream>
#include <random>

using namespace std;

static fftw_plan p_r2c = nullptr;
static fftw_plan p_c2r = nullptr;

void SIP(double* h,int m,double* g,int n,double *z) {



    fftw_complex *XX = fftw_alloc_complex(n/2 + 1);
    fftw_complex *YY = fftw_alloc_complex(n/2 + 1);
    fftw_complex *ZZ = fftw_alloc_complex(n/2 + 1);

    if (p_r2c == nullptr) {
        p_r2c = fftw_plan_dft_r2c_1d(n, z, ZZ, FFTW_ESTIMATE);
        p_c2r = fftw_plan_dft_c2r_1d(n, ZZ, z, FFTW_ESTIMATE);
    }


    
    inv_and_padding_y(n, m, h, YY, p_r2c);
    
    FFT_x(g, XX, p_r2c);

    
    cov_IFFT(XX, YY, ZZ, n, z, p_c2r);

    fftw_free(XX); fftw_free(YY); fftw_free(ZZ);

}


void  inv_and_padding_y(int n,int m,double *y, fftw_complex * YY,fftw_plan p)

{
    double y_temp[n];
    for(int i = 0 ; i < n ; i++ )
    {
        if(i < m )
            y_temp[i] = y[m-i-1];
        else
            y_temp[i] = 0;
    }
    fftw_execute_dft_r2c(p,y_temp,YY);

}

void FFT_x(double *x, fftw_complex (*XX), fftw_plan p)
{
    fftw_execute_dft_r2c(p,x,XX);
}

void cov_IFFT(fftw_complex (*XX), fftw_complex (*YY), fftw_complex (*ZZ), int n, double *z, fftw_plan p)
{
    
    int i_ceil = n/2;
    double inv_n = 1.0/n;

    
    for(int i = 0 ; i <= i_ceil; i++)
    {
        
        ZZ[i][0] = (XX[i][0]*YY[i][0] - XX[i][1]*YY[i][1])*inv_n;
        
        ZZ[i][1] = (XX[i][1]*YY[i][0] + XX[i][0]*YY[i][1])*inv_n;
    }

    fftw_execute_dft_c2r(p,ZZ,z);
}

