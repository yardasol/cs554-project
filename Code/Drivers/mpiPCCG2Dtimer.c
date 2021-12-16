
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../Helpers/poisson_matrix_creaters.c"
#include "../Helpers/matrix_arithmetic.c"
#include "../Helpers/helpers.c"
#include "../Helpers/mpi_helpers.c"
#include "../Helpers/cg_solvers.c"
#include "../Helpers/mpi_cg_solvers.c"
#include "../Helpers/preconditioner_solve.c"
#include "../Helpers/mpi_pccg_solvers.c"


int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    int ntimerf, ntimerc, ntimerd;
    ntimerf = atoi(argv[2]); // No. of times CG-Full solve done- time will be averaged
    ntimerc = atoi(argv[3]); // No. of times CG-CSR solve done- time will be averaged
    ntimerd = atoi(argv[4]); // No. of times CG-DIA solve done- time will be averaged
    int N2D = N*N;
    int dim = 2; // Poisson Matrix type 2D
    int n_diag = 5;
    int offset = 2*(N+1);
    int ntimer = 10;
    float tol = 0.000001;

    struct CGret cgret, cgretcsr, cgretdia;
    struct CGret cgret_p, cgretcsr_p, cgretdia_p;
    clock_t beg, end;
    float t_tot, err; int iters;
    int numprocs, rank, rc;

    struct A_csr A_CSR= create2dPoissonMatCSR(N);
    struct A_csr ILU_CSR = createPoissonILUCSR(N2D, n_diag, offset, &A_CSR);
    struct A_csr L_CSR = getLfromPoissonILUCSR(N2D, n_diag, offset, &ILU_CSR);
    struct A_csr U_CSR = getUfromPoissonILUCSR(N2D, n_diag, offset, &ILU_CSR);
    float* x = create1dRandRHS(N2D);
    float* b = create1dZeroVec(N2D);
    float* b_p = create1dZeroVec(N2D);
    float* b_csr = create1dZeroVec(N2D);
    float* b_p_csr = create1dZeroVec(N2D);
    float* b_dia = create1dZeroVec(N2D);
    float* b_p_dia = create1dZeroVec(N2D);
    float* xsol = create1dZeroVec(N2D);
    float* xguess = create1dZeroVec(N2D);

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    if (numprocs < 2 ) {
      fprintf(stderr, "Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }

    int nwrks, rows;
    nwrks = numprocs-1; // no. of workers

    int* offsvec =  row_load_allot(N2D,numprocs); // vector of offsets for each worker
    if (rank>0) { rows = offsvec[rank]-offsvec[rank-1]; } //Calc. no. of rows of A this worker rank deals with

    // make sure all workers receive the same random vec x created by master
    communicate_xvec(N2D,rank,nwrks,x);

    // create struct with data for parallel mult routines
    struct par_multdat pmdat = parmult_struct_assign(offsvec,rows,rank,N2D,nwrks,n_diag);

    //b_p = mpiMatVecProduct1(pmdat, x, A);
    b_p_csr = mpiMatVecProductCSR1(pmdat, x, A_CSR);
    //b_p_dia = mpiMatVecProductDIA1(pmdat, x, A_DIA);

    /*cgret_p = mpiCGsolveFull(pmdat,A,b_p_csr,xguess,tol);
      cgretcsr_p = mpiCGsolveCSR(pmdat,A_CSR,b_p_csr,xguess,tol);
      cgretdia_p = mpiCGsolveDIA(pmdat,A_DIA,b_p_csr,xguess,tol);*/

    /* CG_FullMPI_timer_output(ntimerf,pmdat,x,xsol,b_p_csr,xguess,A,tol,dim); */
    CG_CSRMPI_timer_output(ntimerc,pmdat,x,xsol,b_p_csr,xguess,A_CSR,tol,dim);
    /* CG_DIAMPI_timer_output(ntimerd,pmdat,x,xsol,b_p_csr,xguess,A_DIA,tol,dim); */

    char pctype[10];
    //Multigrid
    strcpy(pctype, "MG2D"); //pctype = "Jacobi";

    //Jacobi PC
    //strcpy(pctype,"Jacobi"); //pctype = "Jacobi";

    //ILU PC
    //strcpy(pctype,"ILU");

    PCCG_CSRMPI_timer_output(ntimerc,pmdat,x,xsol,b_p_csr,xguess,A_CSR,L_CSR,U_CSR,tol,dim,pctype);

    MPI_Finalize();
    return 0;
}
