#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>

#include "poisson_matrix_creaters.c"

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

struct PCret ILU_PC_Solve(int n, int n_diag, int offset, struct A_csr A, float* r){
    struct PCret pcret;
    struct A_csr ILU, L, U;
	ILU = createPoissonIncompleteLUCSR(n, n_diag, offset, &A);
	L = getLfromPoissonILUCSR(n, n_diag, offset, &ILU);
	U = getUfromPoissonILUCSR(n, n_diag, offset, &ILU);

	int l = 0; //CSR index

	// Ly = r
	float *y = create1dZeroVec(n);
	for (int i = 0; i < n; i++){
		while (L.row_ptr[l+1] == i){
			y[i] -= L.val[l] * y[L.col_ind[l]];
			l++;
		}
		y[i] += r[i];
		y[i] = y[i] / L.val[l];
	}

	// Uz = y;
	float *z = create1dZeroVec(n);
	for (int i = n-1; i >= 0  ; i--){
		while (L.row_ptr[l-1] == i){
			z[i] -= L.val[l] * z[L.col_ind[l]];
			l--;
		}
		z[i] += y[i];
		z[i] = z[i] / L.val[l];
	}
    pcret.sol = z;
    return pcret;
}

/*Preconditioner solve*/
struct PCret PC_Solve(int n, int n_diag, int offset, float* r, struct A_csr A, char* pctype){

    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A,r);
    }
	if (stcmp(pctype, "ILU")==0){
		pcret = ILU_PC_Solve(n, n_diag, offset, A, r);
	}

    return pcret;
}


#endif
