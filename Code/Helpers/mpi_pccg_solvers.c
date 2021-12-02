#ifndef MPI_PCCG_SOLVERS_C
#define MPI_PCCG_SOLVERS_C

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "poisson_matrix_creaters.c"
#include "matrix_arithmetic.c"
#include "helpers.c"
#include "preconditioner_solve.c"
#include "mpi_helpers.c"

struct CGret mpiPCCG_solveCSR(struct par_multdat pmd, struct A_csr A, struct A_csr L, struct A_csr U, float* b, float* xg, float tol, char* pctype)
{
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    struct CGret cgret;
    struct PCret pcret;
    float rsold, rsnew, alpha, rnorm;
    float* Ap = create1dZeroVec(n);
    float* r = create1dZeroVec(n); 
    float* dum = create1dZeroVec(n);
    float* p = create1dZeroVec(n); 
    float* x = create1dZeroVec(n); 
    float* z = create1dZeroVec(n);

    x = VecAdd1(n,xg,r,0);
    dum = mpiMatVecProductCSR1(pmd,xg,A);
    r = VecAdd1(n,b,dum,-1); //r=b-A*x
    pcret = PC_Solve(pmd, n,r,A,L,U,pctype, -1); //Solve M*z0=r0
    z = pcret.sol;
    p = z;

    rsold = innerProd1(n,r,z); //rsold=r'*r
    int iters = 0;
    for (int i=0; i<n; i++){
        Ap = mpiMatVecProductCSR1(pmd,p,A); // Ap=A*p
        alpha = rsold / innerProd1(n,p,Ap); // alpha=rsold/(p'*Ap)
        x = VecAdd1(n,x,p,alpha); // x=x+alpha*p
        r = VecAdd1(n,r,Ap,-alpha); // r=r-alpha*Ap
        rnorm = innerProd1(n,r,r); // rnorm=r'*r;
        if (sqrt(rnorm) < tol) break; // Check convergence
        pcret = PC_Solve(pmd, n,r,A,L,U,pctype,iters); //Solve M*z=r
        z = pcret.sol;
        rsnew = innerProd1(n,r,z); // rsnew = r'*z
        p = VecAdd1(n,z,p,(rsnew / rsold)); // p=z+(rsnew/rsold)*p
        rsold = rsnew;
        iters = i;
    }
    free(r); free(p); free(Ap); free(dum); free(z);
    cgret.x = x;
    cgret.iter = iters;
    return cgret;
}


void PCCG_CSRMPI_timer_output(int ntimer, struct par_multdat pmd, float* x,
        float* xsol, float* b, float* xguess, struct A_csr A_CSR, struct A_csr L_CSR, struct A_csr U_CSR, float tol, int dim, char* pctype){

    int rank = pmd.rank_d;
    int N = pmd.n_d;
    int nwrks = pmd.nwrks_d;
    clock_t beg, end; 
    struct CGret cgret;
    // Start Timer
    beg = clock();

    for(int j=0; j<ntimer; j++){cgret = mpiPCCG_solveCSR(pmd,A_CSR,L_CSR,U_CSR,b,xguess,tol,pctype);}
    xsol = cgret.x;
    int iters = cgret.iter;
    // End timer
    end = clock();
    float t_tot = 1.0*(end-beg)/CLOCKS_PER_SEC;
    float err = VecErrNorm(N, xsol, x);
    //Output
    if(rank==0){  
        printf("CSR Matrix rep (PCCG):- \n");                
        printf("Size %dD: %d on %d worker ranks\n",dim,N,nwrks);
        printf("Error, tol, iters %3f %3f %d \n", err, tol, iters);// Error norm to check if soln converged    
        printf("Total time taken: %f \n", t_tot/ntimer);
        printf("-------------------------------------------------------\n\n");
    }
    return;
}


void cg_Output_verify_v2(struct CGret cgret, struct CGret cgretcsr, struct CGret pccgretcsr,
        float* x, float* b, float** A, int n, int rank)
{
    if (rank == 0)
    {   
        //Output                            
        //print_mat(n,n,A);
        printf("b:\n"); print_vec(n,b);        
        printf("Exact sol x:\n"); print_vec(n,x);
        printf("CGFull sol x:\n"); print_vec(n,cgret.x);
        printf("CGCSR sol x:\n"); print_vec(n,cgretcsr.x);
        printf("PCCGCSR sol x:\n"); print_vec(n,pccgretcsr.x);        
        float errnorm = VecErrNorm(n,cgret.x,x);
        float errnormcsr = VecErrNorm(n,cgretcsr.x,x);
        float errnormpccgcsr = VecErrNorm(n,pccgretcsr.x,x);
        printf("Within %d CGFull iters we converge to error norm: %f \n", cgret.iter, errnorm);
        printf("Within %d CGCSR iters we converge to error norm: %f \n", cgretcsr.iter, errnormcsr);
        printf("Within %d PCCGCSR iters we converge to error norm: %f \n \n", pccgretcsr.iter, errnormpccgcsr);
    }
    return;

}

#endif

