#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double exponential_inverse_cdf(double u, double mu) {
    return -mu * log(1.0 - u);
}

int main() {
    const int num_samples = 10000;
    const double mu = 0.5;
    double U, X;

    for (int i = 0; i < num_samples; ++i) {
        U = (double)rand() / RAND_MAX; // Generate U in [0, 1]
        X = exponential_inverse_cdf(U, mu);
        printf("%lf\n", X); // Print the generated random number
    }

    return 0;
}

