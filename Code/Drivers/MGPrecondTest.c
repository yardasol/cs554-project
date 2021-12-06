#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../Helpers/mg.c"
#include "../Helpers/poisson_matrix_creaters.c"
#include "../Helpers/preconditioner_solve.c"

int main(int argc, char *argv[]) {
    struct PCdata pcdata;

    int n = 11;
    float** A_dense = create2dPoissonMat(n);
    struct A_csr A = create2dPoissonMatCSR(A_dense, n);

    setup_mg_pc(&pcdata, n * n, &A, 2, 1);
    float* r = create1dZeroVec(n*n);
    PC_Solve(n*n, r, A, "MG2D", pcdata);

    print_vec(n*n, r);
}
