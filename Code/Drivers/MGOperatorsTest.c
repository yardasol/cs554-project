#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../Helpers/mg.c"
#include "../Helpers/poisson_matrix_creaters.c"

int main(int argc, char *argv[]) {
    unsigned int n_fine = 3;
    unsigned int n_coarse;
    unsigned int n_rows, n_cols;

    struct A_csr* P = mg_2d_interpolation(n_fine, &n_coarse, &n_rows, &n_cols);
    printf("Interpolation from %u points to %u points, 2D\n", n_coarse, n_fine);
    float** P_dense = csr_to_dense(P, n_rows, n_cols);
    print_mat(n_rows, n_cols, P_dense);

    struct A_csr* R = spmatTransposeCSR(P, n_rows, n_cols);
    printf("Restriction from %u points to %u points, 2D\n", n_fine, n_coarse);
    float** R_dense = csr_to_dense(R, n_cols, n_rows);
    print_mat(n_cols, n_rows, R_dense);

    printf("2D Poisson system\n");
    float** A_dense = create2dPoissonMat(n_fine);
    struct A_csr A = create2dPoissonMatCSR(A_dense, n_fine);
    print_mat(n_fine*n_fine, n_fine*n_fine, A_dense);

    printf("2D coarse Poisson system\n");
    struct A_csr* A_coarse = mg_galerkin_product(R, &A, n_cols, n_rows);
    float** A_coarse_dense = csr_to_dense(A_coarse, n_cols, n_cols);
    print_mat(n_cols, n_cols, A_coarse_dense);

    return 0;
}
