#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>

#include "poisson_matrix_creaters.c"
#include "helpers.c"
#include "mpi_helpers.c"

struct PCret {
    float* diag;
    float* sol;
};

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

struct PCret Jacobi_PC_Solve(int n, struct A_csr A, float* r){
    struct PCret pcret;
    float* z = create1dZeroVec(n);
    pcret.diag = getCSRdiagonal(n,A);
    for(int i = 0; i < n; i++){
        z[i] = r[i]/pcret.diag[i];
    }
    pcret.sol = z;
    return pcret;
}

/*Preconditioner solve*/
struct PCret PC_Solve(struct par_multdat pmd, int n, float* r, struct A_csr A, struct A_csr L, struct A_csr U, char* pctype, int iters){

    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A,r);
    }
	if (strcmp(pctype, "ILU")==0){
	pcret.sol = mpiTriangularSolveCSR1(pmd, r, L, U);
	}

    return pcret;
}


#endif
