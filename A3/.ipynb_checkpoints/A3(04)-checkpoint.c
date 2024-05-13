#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

double gaussian(double x) {
    return exp(-(x * x));
}

int main() {
    int n = 512; // Number of sample points
    float x_min = -5.0, x_max = 5.0, delta = 0.0, *k_values, *x_values; // declaring x_min, x_max, and delta
    fftw_complex *input, *output, *fft_factors, *product;
    FILE *ft_data;
    fftw_plan plan;

    x_values = calloc(n, sizeof(float));
    k_values = calloc(n, sizeof(float));
    input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    fft_factors = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    product = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    plan = fftw_plan_dft_1d(n, input, output, FFTW_FORWARD, FFTW_ESTIMATE);

    // Change file path name here
    ft_data = fopen("fft_data.txt", "w");

    delta = (x_max - x_min) / (n - 1);
    for (int i = 0; i < n; i++) {
        x_values[i] = x_min + i * delta;
        if (i < n / 2)
            k_values[i] = 2 * M_PI * (i / (n * delta));
        else
            k_values[i] = 2 * M_PI * ((i - n) / (n * delta));

        input[i][0] = gaussian(x_values[i]);
        input[i][1] = 0.0;
        fft_factors[i][0] = cos(k_values[i] * x_min);
        fft_factors[i][1] = -sin(k_values[i] * x_min);
    }

    fftw_execute(plan);
    printf("DFT Printing\n");
    for (int i = 0; i < n; i++) {
        printf("%f\t%f\n", output[i][0], output[i][1]);
    }

    // Normalizing
    for (int i = 0; i < n; i++) {
        output[i][0] = (1.0 / sqrt(n)) * output[i][0];
        output[i][1] = (1.0 / sqrt(n)) * output[i][1];
    }

    // Complex multiplication by factors
    for (int i = 0; i < n; i++) {
        product[i][0] = fft_factors[i][0] * output[i][0] - fft_factors[i][1] * output[i][1];
        product[i][1] = fft_factors[i][0] * output[i][1] + fft_factors[i][1] * output[i][0];
    }

    // Constructing FT
    for (int i = 0; i < n; i++) {
        output[i][0] = sqrt(n / (2 * M_PI)) * delta * product[i][0];
        output[i][1] = sqrt(n / (2 * M_PI)) * delta * product[i][1];
    }

    fprintf(ft_data, "#X_val f(x)  k_val   FT(f(x))\n");
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            printf("Creating file\n");
        }
        fprintf(ft_data, "%f %f %f  %f\n", x_values[i], input[i][0], k_values[i], output[i][0]);
    }
    
    fclose(ft_data);
    fftw_destroy_plan(plan);
    fftw_free(input);
    fftw_free(output);
    fftw_free(fft_factors);
    fftw_free(product);
    free(x_values);
    free(k_values);

    return 0;
}

