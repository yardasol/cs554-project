
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
    float** A = create1dPoissonMat(N);
    struct A_csr A_CSR= create1dPoissonMatCSR(N);
    struct A_dia A_DIA= create1dPoissonMatDIA(N);
    int n_diag = 3;
    float* x = create1dRandRHS(N);
    float* b = create1dZeroVec(N);  float* b_p = create1dZeroVec(N);         
    float* b_csr = create1dZeroVec(N);  float* b_p_csr = create1dZeroVec(N);
    float* b_dia = create1dZeroVec(N);  float* b_p_dia = create1dZeroVec(N);
    float* xsol = create1dZeroVec(N); 
    float* xguess = create1dZeroVec(N);   
    float tol = 0.000001;
    struct CGret cgret, cgretcsr, cgretdia; 
    struct CGret cgret_p, cgretcsr_p, cgretdia_p, pccgretcsr_p;     
    clock_t beg, end; 
    float t_tot, err; int iters;
    char pctype[10];
    
    int numprocs, rank, rc;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);


    if (numprocs < 2 ) {
        printf("Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }

    int nwrks, rows;
    nwrks = numprocs-1; // no. of workers

    int* offsvec =  row_load_allot(N,numprocs); // vector of offsets for each worker
    if (rank>0) { rows = offsvec[rank]-offsvec[rank-1]; } //Calc. no. of rows of A this worker rank deals with

    // make sure all workers receive the same random vec x created by master 
    communicate_xvec(N,rank,nwrks,x); 

    // create struct with data for parallel mult routines
    struct par_multdat pmdat = parmult_struct_assign(offsvec,rows,rank,N,nwrks,n_diag); 

    //b_p = mpiMatVecProduct1(pmdat, x, A);
    b_p_csr = mpiMatVecProductCSR1(pmdat, x, A_CSR);        
    //b_p_dia = mpiMatVecProductDIA1(pmdat, x, A_DIA);    
    
    cgret_p = mpiCGsolveFull(pmdat,A,b_p_csr,xguess,tol);
    cgretcsr_p = mpiCGsolveCSR(pmdat,A_CSR,b_p_csr,xguess,tol);
    cgretdia_p = mpiCGsolveDIA(pmdat,A_DIA,b_p_csr,xguess,tol);

    //Jacobi PC
    strcpy(pctype,"Jacobi"); //pctype = "Jacobi";
    pccgretcsr_p = mpiPCCG_solveCSR(pmdat,A_CSR,b_p_csr,xguess,tol,pctype);

    
    /*mult_Output_verify(N,n_diag,numprocs,rank,offsvec,
        A,A_CSR,A_DIA,x,b,b_csr,b_dia,b_p,b_p_csr,b_p_dia); // Print output on rank 0 to validate mult routines
    */
    
    //Output
    //cg_Output_verify(cgret_p,cgretcsr_p,cgretdia_p,x,b_p_csr,A,N,rank);
    cg_Output_verify_v2(cgret_p,cgretcsr_p,pccgretcsr_p,x,b_p_csr,A,N,rank);


    MPI_Finalize();
    return 0;
}




