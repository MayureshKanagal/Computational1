#include <stdio.h>
#include <gsl/gsl_linalg.h>

void lu_decomposition(double A[], int size, double L[], double U[]) {
    gsl_matrix_view gsl_A = gsl_matrix_view_array(A, size, size);
    gsl_permutation *p = gsl_permutation_alloc(size);
    int signum;

    gsl_linalg_LU_decomp(&gsl_A.matrix, p, &signum);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i > j)
                L[i * size + j] = gsl_matrix_get(&gsl_A.matrix, i, j);
            else if (i == j)
                L[i * size + j] = 1.0;
            else
                L[i * size + j] = 0.0;

            if (i <= j)
                U[i * size + j] = gsl_matrix_get(&gsl_A.matrix, i, j);
            else
                U[i * size + j] = 0.0;
        }
    }

    gsl_permutation_free(p);
}

int main() {
    // Example matrix of size 3x3
    double A[] = {4, 1, 1, 0, 1, -1, -3, 1, 1, 0, 2, 1, 5, -1, -1, -1, -1 -1, 4, 0, 0, 2, -1, 1, 4};
    int size = 5;

    // Allocate space for L and U matrices
    double L[size * size];
    double U[size * size];

    // Perform LU decomposition
    lu_decomposition(A, size, L, U);

    // Print the results
    printf("Matrix L:\n");
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%8.7f ", L[i * size + j]);
        }
        printf("\n");
    }

    printf("Matrix U:\n");
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%8.7f ", U[i * size + j]);
        }
        printf("\n");
    }

    return 0;
}