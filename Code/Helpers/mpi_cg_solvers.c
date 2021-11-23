#ifndef MPI_CG_SOLVERS_C
#define MPI_CG_SOLVERS_C

#include <time.h>
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
#include "cg_solvers.c"
#include "mpi_helpers.c"

struct CGret mpiCGsolveFull(struct par_multdat pmd, float** A, float* b, float* xg, float tol)
{
    int source, dest, rows, offset, i, j, ind, ioff;
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    int numworkers = pmd.nwrks_d;
    int* offsv = pmd.offsetv_d;  
    offset = pmd.offsetv_d[rank-1];
    rows = pmd.rows_d;

    struct CGret cgret;
    float rsold, rsnew, alpha;
    float* Ap = create1dZeroVec(n);
    float* r = create1dZeroVec(n); 
    float* dum = create1dZeroVec(n);
    float* p = create1dZeroVec(n); 
    float* x = create1dZeroVec(n); 

    x = VecAdd1(n,xg,r,0);
    dum = mpiMatVecProduct1(pmd,xg,A);
    r = VecAdd1(n,b,dum,-1); //r=b-A*x
    p = r;
    rsold = innerProd1(n,r,r); //rsold=r'*r
    int iters = 0;
    for (int i=0; i<n; i++){
        Ap = mpiMatVecProduct1(pmd,p,A); // Ap=A*p
        alpha = rsold / innerProd1(n,p,Ap); // alpha=rsold/(p'*Ap)
        x = VecAdd1(n,x,p,alpha); // x=x+alpha*p
        r = VecAdd1(n,r,Ap,-alpha); // r=r-alpha*Ap
        rsnew = innerProd1(n,r,r); // rsnew=r'*r;
        if (sqrt(rsnew) < tol) break; // Check convergence
        p = VecAdd1(n,r,p,(rsnew / rsold)); // p=r+(rsnew/rsold)*p
        rsold = rsnew;
        iters = i;
    }
    free(r); free(p); free(Ap); free(dum);
    cgret.x = x;
    cgret.iter = iters;
    return cgret;
}

struct CGret mpiCGsolveFull_allpar(struct par_multdat pmd, float** A, float* b, float* xg, float tol)
{
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    struct CGret cgret;
    float rsold, rsnew, alpha;
    float* Ap = create1dZeroVec(n);
    float* r = create1dZeroVec(n); 
    float* dum = create1dZeroVec(n);
    float* p = create1dZeroVec(n); 
    float* x = create1dZeroVec(n); 

    x = mpiVecAdd1(pmd,xg,r,0);
    dum = mpiMatVecProduct1(pmd,xg,A);
    r = mpiVecAdd1(pmd,b,dum,-1); //r=b-A*x
    p = r;
    rsold = mpiInnerProd1(pmd,r,r); //rsold=r'*r
    int iters = 0;
    for (int i=0; i<n; i++){
        Ap = mpiMatVecProduct1(pmd,p,A); // Ap=A*p
        alpha = rsold / mpiInnerProd1(pmd,p,Ap); // alpha=rsold/(p'*Ap)
        x = mpiVecAdd1(pmd,x,p,alpha); // x=x+alpha*p
        r = mpiVecAdd1(pmd,r,Ap,-alpha); // r=r-alpha*Ap
        rsnew = mpiInnerProd1(pmd,r,r); // rsnew=r'*r;
        if (sqrt(rsnew) < tol) break; // Check convergence
        p = mpiVecAdd1(pmd,r,p,(rsnew / rsold)); // p=r+(rsnew/rsold)*p
        rsold = rsnew;
        iters = i;
    }
    free(r); free(p); free(Ap); free(dum);
    cgret.x = x;
    cgret.iter = iters;
    return cgret;
}

struct CGret mpiCGsolveCSR(struct par_multdat pmd, struct A_csr A, float* b, float* xg, float tol)
{
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    struct CGret cgret;
    float rsold, rsnew, alpha;
    float* Ap = create1dZeroVec(n);
    float* r = create1dZeroVec(n); 
    float* dum = create1dZeroVec(n);
    float* p = create1dZeroVec(n); 
    float* x = create1dZeroVec(n); 

    x = VecAdd1(n,xg,r,0);
    dum = mpiMatVecProductCSR1(pmd,xg,A);
    r = VecAdd1(n,b,dum,-1); //r=b-A*x
    p = r;
    rsold = innerProd1(n,r,r); //rsold=r'*r
    int iters = 0;
    for (int i=0; i<n; i++){
        Ap = mpiMatVecProductCSR1(pmd,p,A); // Ap=A*p
        alpha = rsold / innerProd1(n,p,Ap); // alpha=rsold/(p'*Ap)
        x = VecAdd1(n,x,p,alpha); // x=x+alpha*p
        r = VecAdd1(n,r,Ap,-alpha); // r=r-alpha*Ap
        rsnew = innerProd1(n,r,r); // rsnew=r'*r;
        if (sqrt(rsnew) < tol) break; // Check convergence
        p = VecAdd1(n,r,p,(rsnew / rsold)); // p=r+(rsnew/rsold)*p
        rsold = rsnew;
        iters = i;
    }
    free(r); free(p); free(Ap); free(dum);
    cgret.x = x;
    cgret.iter = iters;
    return cgret;
}

// Extract diagonal from CSR matrix
float* getCSRdiagonal(int n, struct A_csr A){
    float* d = create1dZeroVec(n);
    int col=0;
    for (int i = 0; i<n; i++){   
        for (int j=A.row_ptr[i]; j<A.row_ptr[i+1]; j++){
            col = A.col_ind[j];
            if(col==i){
                d[col] = A.val[j];
                col++; 
                if(col>n){break;}   
                }
            }
        }     
    return d; 
}

/*Preconditioner solve*/
float *PC_Solve(int n, float* r, struct A_csr A, char* pctype){
    float* zero = create1dZeroVec(n);
    float* z = create1dZeroVec(n);
    float* diag = create1dZeroVec(n);
    if (strcmp(pctype,"Jacobi")==0){
        //z = VecAdd1(n,zero,r,2);
        diag = getCSRdiagonal(n,A);
        for(int i = 0; i < n; i++){
            z[i] = r[i]/diag[i];
        }

    }
    return z;
}

struct CGret mpiPCCG_solveCSR(struct par_multdat pmd, struct A_csr A, float* b, float* xg, float tol, char* pctype)
{
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    struct CGret cgret;
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
    z = PC_Solve(n,r,A,pctype); //Solve M*z0=r0
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
        z = PC_Solve(n,r,A,pctype); //Solve M*z=r
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

struct CGret mpiCGsolveDIA(struct par_multdat pmd, struct A_dia A, float* b, float* xg, float tol)
{
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    struct CGret cgret;
    float rsold, rsnew, alpha;
    float* Ap = create1dZeroVec(n);
    float* r = create1dZeroVec(n); 
    float* dum = create1dZeroVec(n);
    float* p = create1dZeroVec(n); 
    float* x = create1dZeroVec(n); 

    x = VecAdd1(n,xg,r,0);
    dum = mpiMatVecProductDIA1(pmd,xg,A);
    r = VecAdd1(n,b,dum,-1); //r=b-A*x
    p = r;
    rsold = innerProd1(n,r,r); //rsold=r'*r
    int iters = 0;
    for (int i=0; i<n; i++){
        Ap = mpiMatVecProductDIA1(pmd,p,A); // Ap=A*p
        alpha = rsold / innerProd1(n,p,Ap); // alpha=rsold/(p'*Ap)
        x = VecAdd1(n,x,p,alpha); // x=x+alpha*p
        r = VecAdd1(n,r,Ap,-alpha); // r=r-alpha*Ap
        rsnew = innerProd1(n,r,r); // rsnew=r'*r;
        if (sqrt(rsnew) < tol) break; // Check convergence
        p = VecAdd1(n,r,p,(rsnew / rsold)); // p=r+(rsnew/rsold)*p
        rsold = rsnew;
        iters = i;
    }
    free(r); free(p); free(Ap); free(dum);
    cgret.x = x;
    cgret.iter = iters;
    return cgret;
}

void cg_Output_verify(struct CGret cgret, struct CGret cgretcsr, struct CGret cgretdia,
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
        printf("CGDIA sol x:\n"); print_vec(n,cgretdia.x);        
        float errnorm = VecErrNorm(n,cgret.x,x);
        float errnormcsr = VecErrNorm(n,cgretcsr.x,x);
        float errnormdia = VecErrNorm(n,cgretdia.x,x);
        printf("Within %d CGFull iters we converge to error norm: %f \n", cgret.iter, errnorm);
        printf("Within %d CGCSR iters we converge to error norm: %f \n", cgretcsr.iter, errnormcsr);
        printf("Within %d CGDIA iters we converge to error norm: %f \n \n", cgretdia.iter, errnormdia);
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

void CG_FullMPI_timer_output(int ntimer, struct par_multdat pmd, float* x,
 float* xsol, float* b, float* xguess, float** A, float tol, int dim){
    
    int rank = pmd.rank_d;
    int N = pmd.n_d;
    int nwrks = pmd.nwrks_d;
    clock_t beg, end; 
    struct CGret cgret;
    // Start Timer
    beg = clock();
        
    for(int j=0; j<ntimer; j++){cgret = mpiCGsolveFull(pmd,A,b,xguess,tol);}
    xsol = cgret.x;
    int iters = cgret.iter;
    // End timer
    end = clock();
    float t_tot = 1.0*(end-beg)/CLOCKS_PER_SEC;
    float err = VecErrNorm(N, xsol, x);
    //Output
    if(rank==0){  
        printf("Full Matrix rep:- \n");        
        printf("Size %dD: %d on %d worker ranks\n",dim,N,nwrks);
        printf("Error, tol, iters %3f %3f %d \n", err, tol, iters);// Error norm to check if soln converged    
        printf("Total time taken: %f \n", t_tot/ntimer);
        printf("-------------------------------------------------------\n\n");
    }
    return;
}

void CG_CSRMPI_timer_output(int ntimer, struct par_multdat pmd, float* x,
 float* xsol, float* b, float* xguess, struct A_csr A_CSR, float tol, int dim){
    
    int rank = pmd.rank_d;
    int N = pmd.n_d;
    int nwrks = pmd.nwrks_d;
    clock_t beg, end; 
    struct CGret cgret;
    // Start Timer
    beg = clock();
        
    for(int j=0; j<ntimer; j++){cgret = mpiCGsolveCSR(pmd,A_CSR,b,xguess,tol);}
    xsol = cgret.x;
    int iters = cgret.iter;
    // End timer
    end = clock();
    float t_tot = 1.0*(end-beg)/CLOCKS_PER_SEC;
    float err = VecErrNorm(N, xsol, x);
    //Output
    if(rank==0){  
        printf("CSR Matrix rep:- \n");                
        printf("Size %dD: %d on %d worker ranks\n",dim,N,nwrks);
        printf("Error, tol, iters %3f %3f %d \n", err, tol, iters);// Error norm to check if soln converged    
        printf("Total time taken: %f \n", t_tot/ntimer);
        printf("-------------------------------------------------------\n\n");
    }
    return;
}

void PCCG_CSRMPI_timer_output(int ntimer, struct par_multdat pmd, float* x,
 float* xsol, float* b, float* xguess, struct A_csr A_CSR, float tol, int dim, char* pctype){
    
    int rank = pmd.rank_d;
    int N = pmd.n_d;
    int nwrks = pmd.nwrks_d;
    clock_t beg, end; 
    struct CGret cgret;
    // Start Timer
    beg = clock();
        
    for(int j=0; j<ntimer; j++){cgret = mpiPCCG_solveCSR(pmd,A_CSR,b,xguess,tol,pctype);}
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

void CG_DIAMPI_timer_output(int ntimer, struct par_multdat pmd, float* x,
 float* xsol, float* b, float* xguess, struct A_dia A_DIA, float tol, int dim){
    
    int rank = pmd.rank_d;
    int N = pmd.n_d;
    int nwrks = pmd.nwrks_d;
    clock_t beg, end; 
    struct CGret cgret;
    // Start Timer
    beg = clock();
        
    for(int j=0; j<ntimer; j++){cgret = mpiCGsolveDIA(pmd,A_DIA,b,xguess,tol);}
    xsol = cgret.x;
    int iters = cgret.iter;
    // End timer
    end = clock();
    float t_tot = 1.0*(end-beg)/CLOCKS_PER_SEC;
    float err = VecErrNorm(N, xsol, x);
    //Output
    if(rank==0){  
        printf("DIA Matrix rep:- \n");                
        printf("Size %dD: %d on %d worker ranks\n",dim,N,nwrks);
        printf("Error, tol, iters %3f %3f %d \n", err, tol, iters);// Error norm to check if soln converged    
        printf("Total time taken: %f \n", t_tot/ntimer);
        printf("-------------------------------------------------------\n\n");
    }
    return;
}

#endif