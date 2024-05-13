#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

double sinc(double x) {
    if (x == 0.0) {
        return 1.0;
    } else {
        return (sin(x) / x);
    }
}

int main() {
    int n = 512; // Number of sample points
    float x_min = -50.0, x_max = 50.0, d = 0.0;
    float *k_vals, *x_vals; // Declare x_vals, k_vals instead of x_arr, k_arr
    fftw_complex *in_data, *out_data, *ft_factors, *product;
    FILE *ft_file;
    fftw_plan fft_plan;

    x_vals = calloc(n, sizeof(float));
    k_vals = calloc(n, sizeof(float));
    in_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    ft_factors = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    product = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fft_plan = fftw_plan_dft_1d(n, in_data, out_data, FFTW_FORWARD, FFTW_ESTIMATE);

    ft_file = fopen("C:/Users/sandesh/Desktop/mayuresh/Tifr sem2/computational/Assignment 3/fftw_1_data.txt", "w");

    d = (x_max - x_min) / (n - 1);
    for (int i = 0; i < n; i++) {
        x_vals[i] = x_min + i * d;
        if (i < n / 2) {
            k_vals[i] = 2 * M_PI * (i / (n * d));
        } else {
            k_vals[i] = 2 * M_PI * ((i - n) / (n * d));
        }

        in_data[i][0] = sinc(x_vals[i]);
        in_data[i][1] = 0.0;
        ft_factors[i][0] = cos(k_vals[i] * x_min);
        ft_factors[i][1] = -sin(k_vals[i] * x_min);
    }

    fftw_execute(fft_plan);

    printf("DFT Printing\n");
    for (int i = 0; i < n; i++) {
        printf("%f  %f\n", out_data[i][0], out_data[i][1]);
    }

    // Normalizing
    for (int i = 0; i < n; i++) {
        out_data[i][0] = (1.0 / sqrt(n)) * out_data[i][0];
        out_data[i][1] = (1.0 / sqrt(n)) * out_data[i][1];
    }

    // Complex multiplication by factors
    for (int i = 0; i < n; i++) {
        product[i][0] = ft_factors[i][0] * out_data[i][0] - ft_factors[i][1] * out_data[i][1];
        product[i][1] = ft_factors[i][0] * out_data[i][1] + ft_factors[i][1] * out_data[i][0];
    }

    // Constructing FT
    for (int i = 0; i < n; i++) {
        out_data[i][0] = sqrt(n / (2 * M_PI)) * d * product[i][0];
        out_data[i][1] = sqrt(n / (2 * M_PI)) * d * product[i][1];
    }

    fprintf(ft_file, "#X_val f(x)  k_val   FT(f(x))\n");
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            printf("Creating file\n");
        }
        fprintf(ft_file, "%f %f %f  %f\n", x_vals[i], in_data[i][0], k_vals[i], out_data[i][0]);
    }

    fclose(ft_file);

    return 0;
}
