#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../Helpers/matrix_arithmetic.c"
#include "../Helpers/mpi_helpers.c"
#include "../Helpers/mg.c"
#include "../Helpers/poisson_matrix_creaters.c"
#include "../Helpers/preconditioner_solve.c"

int main(int argc, char *argv[]) {
    struct PCdata pcdata;

    int numprocs, rank, rc;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    int n = 11;
    float** A_dense = create1dPoissonMat(n);
    struct A_csr A = create1dPoissonMatCSR(n);

    struct par_multdat pmdat = parmult_create(n);

    setup_mg_pc(&pcdata, n, &A, -1, 0);
    srand(0);
    float* r = create1dZeroVec(n);
    for (unsigned int i = 0; i < n; i++) {
        r[i] = rand() % 10;
    }

    struct PCret pc = PC_Solve(pmdat, n, r, A, A, A, "MG1D", pcdata);

    if (rank == 0) {
        printf("RHS\n");
        print_vec(n, r);

        printf("Soln\n");
        print_vec(n, pc.sol);
    }

    MPI_Finalize();
    return 0;
}
