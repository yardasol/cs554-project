#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>

#include "poisson_matrix_creaters.c"

struct PCret {
    float* diag;
    float* sol;
};

// Extract diagonal from CSR matrix
float* getCSRdiagonal(int n, struct A_csr *A_ptr){
    float* d = create1dZeroVec(n);
    int col=0;
    for (int i = 0; i<n; i++){   
        for (int j=A_ptr->row_ptr[i]; j<A_ptr->row_ptr[i+1]; j++){
            col = A_ptr->col_ind[j];
            if(col==i){
                d[col] = A_ptr->val[j];
                col++; 
                if(col>n){break;}   
                }
            }
        }     
    return d; 
}

struct PCret Jacobi_PC_Solve(int n, struct A_csr *A_ptr, float* r){
    struct PCret pcret;
    float* z = create1dZeroVec(n);
    pcret.diag = getCSRdiagonal(n,A_ptr);
    for(int i = 0; i < n; i++){
        z[i] = r[i]/pcret.diag[i];
    }
    pcret.sol = z;
    return pcret;
}

struct PCret ILU_PC_Solve(int n, struct A_csr *L_ptr, struct A_csr *U_ptr, float* r){
    struct PCret pcret;
	int l = 0; //CSR index

	// Ly = r
	float *y = create1dZeroVec(n);
	for (int i = 0; i < n; i++){
		while (L_ptr->row_ptr[l+1] == i){
			y[i] -= L_ptr->val[l] * y[L_ptr->col_ind[l]];
			l++;
		}
		y[i] += r[i];
		y[i] = y[i] / L_ptr->val[l];
	}

	// Uz = y;
	float *z = create1dZeroVec(n);
	for (int i = n-1; i >= 0  ; i--){
		while (U_ptr->row_ptr[l-1] == i){
			z[i] -= U_ptr->val[l] * z[U_ptr->col_ind[l]];
			l--;
		}
		z[i] += y[i];
		z[i] = z[i] / U_ptr->val[l];
	}
    pcret.sol = z;
    return pcret;
}

/*Preconditioner solve*/
struct PCret PC_Solve(int n, float* r, struct A_csr *A_ptr, struct A_csr *L_ptr, struct A_csr *U_ptr, char* pctype){

    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A_ptr,r);
    }
	if (stcmp(pctype, "ILU")==0){
		pcret = ILU_PC_Solve(n, L_ptr, U_ptr, r);
	}

    return pcret;
}


#endif
