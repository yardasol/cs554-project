#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>

#include "poisson_matrix_creaters.c"

struct PCret {
    float* diag;
    float* sol;
};

struct PCdata {
    float* diag;
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

// Setup data in PCdata struct for info to be extracted before solve 
struct PCdata setupPCdata(struct PCdata pcdata, char* pctype, int n, struct A_csr A)
{
    if (strcmp(pctype,"Jacobi")==0){
        pcdata.diag = getCSRdiagonal(n,A);
    }
    return pcdata;
}



struct PCret Jacobi_PC_Solve(int n, struct A_csr A, float* r, struct PCdata pcdata){
    struct PCret pcret;
    float* z = create1dZeroVec(n);
    //pcret.diag = getCSRdiagonal(n,A);
    for(int i = 0; i < n; i++){
        z[i] = r[i]/pcdata.diag[i];
    }
    pcret.sol = z;
    return pcret;
}


/*Preconditioner solve*/
struct PCret PC_Solve(int n, float* r, struct A_csr A, char* pctype, struct PCdata pcdata){

    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A,r, pcdata);
    }

    return pcret;
}


#endif
