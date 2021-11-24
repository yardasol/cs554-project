#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>

#include "poisson_matrix_creaters.c"

struct PCret {
    float* diag;
    float* sol;
};

// Get ILU in-place factorization
// hardcoded for 1D poisson matrix structure right now
struct A_csr* getILU(int n, struct A_csr A){
    int i, nnz, n_diag;
    struct A_csr L, U;
    n_diag = 2; //hardcoded
    nnz = n*n_diag - 1; //-1 hardcoded

    n_diag_A = 3;
    nnz_A = n*n_diag_A - 2;

    L.val = malloc(sizeof(float*) * (nnz));
    L.col_ind = malloc(sizeof(int*) * (nnz));
    L.row_ptr = malloc(sizeof(int*) * (nnz + 1));
    L.val[0] = 1; L.col_ind[0] = 0; L.row_ptr[0] = 0;

    U.val = malloc(sizeof(float*) * (nnz));
    U.col_ind = malloc(sizeof(int*) * (nnz));
    U.row_ptr = malloc(sizeof(int*) * (nnz + 1));
    for (int l = 0; l < n_diag; l++){
        U.val[l] = A.val[l];
	U.col_ind[l] = A.col_ind[l];
	U.row_ptr[l] = A.row_ptr[l];
    }

    int cnt_L = 1;
    int cnt_U = n_diag; 
    int cnt_A = n_diag;
    int i,j,k;
    // ijk version
    for (i = 1; i < n; i++){
	for (k = 0; k < i; k++){
      	    c = //current a index //Update lower diag of L
	    d = (i+1)*n_diag_A - k//index of A diagonal
	    L.val[cnt_L] = A.val[l]/A.val[d];
	    L.col_ind[cnt_L] = A.col_ind[l];
	    L.row_prt[cnt_L] = A.row_ptr[l];
	    for (j = k+1; k < n; k++){
	      // update U using L
	    }
	}
	//add 1 to main diag of L
	index_A = ...
	L.val[cnt_L] = 1;
	L.col_ind[cnt_L] = A.col_ind[index_A];
	L.row_prt[cnt_L] = A.row_ptr[index_A];
	cnt_L += 1;
    }

    //single index version
    for (l = 1; l < nnz_A; l++){
	//update L
	if (A.col_ind[l] > A.row_ptr[l]){
	  d = //index of diagonal
	  L.val[cnt_L] = A.val[l]/A.val[d];
	  L.col_ind[cnt_L] = A.col_ind[l];
	  L.row_prt[cnt_L] = A.row_ptr[l];
	}
	//update U
	else if (A.col_ind[l] <= A.row_ptr[l]){		
	  U.val[cnt_U] = A.val[l]
	}
    }
}

/* Create an ILU-preconditioned A */
//float **create1dILUPoissonMat(int n, float **A){
//    int i, k, j;

    /* need to add some lines that multiply the L and U components */

//    for(i = 1; i < n; i++){
//   		
//      for(k = 0; k < i-1; k++){
//        if (A[i][k] != 0){
//          L[i][k] = A[i][k] / A[k][k];
//            for(j = k + 1; j < n; j++){
//              if (A[i][j] != 0){
//                U[i][j] = A[i][j] - L[i][k] * U[k][j]
//              }
//            } 
//        }
//      }
//    }
//    return A;
//}

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

struct PCret ILU_PC_Solve(int n, struct A_csr A, float* r){
    struct PCret pcret;
    struct A_csr ILU = ...;
    float* z = create1dZeroVec(n);
    pcret.diag = getILU(n,A);
    for(int i = 0; i < n; i++){
        z[i] = r[i]/pcret.diag[i];
    }
    pcret.sol = z;
    return pcret;
}

/*Preconditioner solve*/
struct PCret PC_Solve(int n, float* r, struct A_csr A, char* pctype){

    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A,r);
    }

    return pcret;
}


#endif
