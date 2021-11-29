#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>

#include "poisson_matrix_creaters.c"
#include "helpers.c"

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

struct PCret ILU_PC_Solve(int n, struct A_csr L, struct A_csr U, float* r){
    struct PCret pcret;
	int l = 0; //CSR index
	
	// Ly = r
	int row, col;
	int row_ptr, next_row_ptr;
	int i;
	float *y = create1dZeroVec(n);
        for (row = 0; row < n; row++){
            row_ptr = L.row_ptr[row];
            next_row_ptr = L.row_ptr[row+1];
            for (i = row_ptr; i < next_row_ptr; i++){
                col = L.col_ind[i];
                y[row] -= L.val[i] * y[col];
            }
            y[row] += r[row];
            y[row] = y[row] / L.val[i-1];
        }

	// Uz = y;
	float *z = create1dZeroVec(n);
	for (row = n-1; row >= 0; row--){
		row_ptr = U.row_ptr[row+1];
		next_row_ptr = U.row_ptr[row];
		for (i = row_ptr-1; i >= next_row_ptr; i--){
			col = U.col_ind[i];
			z[row] -= U.val[i] * z[col];
		}
		z[row] += y[row];
		z[row] = z[row] / U.val[i+1];
	}
    pcret.sol = z;
    return pcret;
}

/*Preconditioner solve*/
struct PCret PC_Solve(int n, float* r, struct A_csr A, struct A_csr L, struct A_csr U, char* pctype){

    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A,r);
    }
	if (strcmp(pctype, "ILU")==0){
			pcret = ILU_PC_Solve(n, L, U, r);
	}

    return pcret;
}


#endif
