#ifndef CG_SOLVERS_C
#define CG_SOLVERS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "poisson_matrix_creaters.c"
#include "matrix_arithmetic.c"

struct CGret {
    int iter;
    float* x;
};

struct CGret CGsolveFull(int n, float** A, float* b, float* xg, float tol)
{
    struct CGret cgret;
    float rsold, rsnew, alpha;
    float* Ap;
    float* r = create1dZeroVec(n); 
    float* dum = create1dZeroVec(n);
    float* p = create1dZeroVec(n); Ap = create1dZeroVec(n); 
    float* x;
    x = malloc(sizeof(float*) * n);
    x = VecAdd1(n,xg,r,0);
    dum = matVecProduct1(n,xg,A);
    r = VecAdd1(n,b,dum,-1); //r=b-A*x
    p = r;
    rsold = innerProd1(n,r,r); //rsold=r'*r
    int iters;
    for (int i=0; i<n; i++){
        Ap = matVecProduct1(n,p,A); // Ap=A*p
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

struct CGret CGsolveCSR(int n, struct A_csr A, float* b, float* xg, float tol)
{
    struct CGret cgret;
    float rsold, rsnew, alpha;
    float* Ap;
    float* r = create1dZeroVec(n); 
    float* dum = create1dZeroVec(n);
    float* p = create1dZeroVec(n); Ap = create1dZeroVec(n); 
    float* x;
    x = malloc(sizeof(float*) * n);
    x = VecAdd1(n,xg,r,0);
    dum = matVecProductCSR1(n,xg,&A);
    r = VecAdd1(n,b,dum,-1); //r=b-A*x
    p = r;
    rsold = innerProd1(n,r,r); //rsold=r'*r
    int iters;
    for (int i=0; i<n; i++){
        Ap = matVecProductCSR1(n,p,&A); // Ap=A*p
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

struct CGret CGsolveDIA_1D(int n, int ndiag, struct A_dia A, float* b, float* xg, float tol)
{
    struct CGret cgret;
    float rsold, rsnew, alpha;
    float* Ap;
    float* r = create1dZeroVec(n); 
    float* dum = create1dZeroVec(n);
    float* p = create1dZeroVec(n); Ap = create1dZeroVec(n); 
    float* x;
    x = malloc(sizeof(float*) * n);
    x = VecAdd1(n,xg,r,0);
    dum = matVecProductDIA1(n,ndiag,xg,A);
    r = VecAdd1(n,b,dum,-1); //r=b-A*x
    p = r;
    rsold = innerProd1(n,r,r); //rsold=r'*r
    int iters;
    for (int i=0; i<n; i++){
        Ap = matVecProductDIA1(n,ndiag,p,A); // Ap=A*p
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

#endif
