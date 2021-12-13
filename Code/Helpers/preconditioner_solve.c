#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>
#include <stdbool.h>

#include "poisson_matrix_creaters.c"
#include "mg.c"
#include "helpers.c"
#include "mpi_helpers.c"

struct PCret {
    float* diag;
    float* sol;
};

struct PCdata {
    float* diag;
    struct MGData* mg;
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

void setup_mg_pc(struct PCdata* pcdata, int n, struct A_csr Aval, int num_levels, bool is_2d) {
    if (num_levels < 2) {
        if (is_2d) {
            num_levels = (int) log2f(sqrt(n));
        } else {
            num_levels = (int) log2f(n);
        }
    }

    pcdata->mg = calloc(1, sizeof(struct MGData));
    pcdata->mg->num_levels = num_levels;
    pcdata->mg->levels = calloc(num_levels, sizeof(struct MGLevelData));
    pcdata->mg->num_pre_relax = 4;
    pcdata->mg->num_post_relax = 4;
    pcdata->mg->num_coarsest_relax = 4;

    struct A_csr* A = malloc(sizeof(struct A_csr));
    memcpy(A, &Aval, sizeof(struct A_csr));

    struct MGLevelData* levels = pcdata->mg->levels;
    /* Set up first MG level */
    levels[0].A = A;
    levels[0].diag = getCSRdiagonal(n, *A);
    if (is_2d) {
        levels[0].n_pts_per_dim = (int) sqrt(n);
    } else {
        levels[0].n_pts_per_dim = n;
    }
    levels[0].n_pts = n;
    levels[0].pmdat = parmult_create(n);

    for (unsigned int i = 1; i < num_levels; i++) {
        unsigned int n_coarse_per_dim, n_coarse;
        struct A_csr* AH, *P, *R;

        /* Create interpolation from this level to higher (more fine) level */
        if (is_2d) {
            P = mg_2d_interpolation(levels[i-1].n_pts_per_dim, &n_coarse_per_dim, NULL, &n_coarse);
        } else {
            P = mg_1d_interpolation(levels[i-1].n_pts_per_dim, &n_coarse_per_dim);
            n_coarse = n_coarse_per_dim;
        }

        R = spmatTransposeCSR(P, levels[i-1].n_pts, n_coarse);
        AH = mg_galerkin_product(R, levels[i-1].A, n_coarse, levels[i-1].n_pts);

        levels[i].A = AH;
        levels[i].diag = getCSRdiagonal(n_coarse, *AH);
        levels[i].R = R;
        levels[i].P = P;
        levels[i].n_pts_per_dim = n_coarse_per_dim;
        levels[i].n_pts = n_coarse;
        levels[i].pmdat = parmult_create(n_coarse);
    }
}

// Setup data in PCdata struct for info to be extracted before solve
struct PCdata setupPCdata(struct PCdata pcdata, char* pctype, int n, struct A_csr A)
{
    if (strcmp(pctype,"Jacobi")==0){
        pcdata.diag = getCSRdiagonal(n,A);
    } else if (strcmp(pctype, "MG2D") == 0) {
        setup_mg_pc(&pcdata, n, A, -1, 1);
    } else if (strcmp(pctype, "MG1D") == 0) {
        setup_mg_pc(&pcdata, n, A, -1, 0);
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

void mg_solve_recursive(struct par_multdat pmd, struct MGData* mg, struct MGLevelData* levels, unsigned int level, float* b, float* x) {
    if (level >= mg->num_levels) {
        return;
    }
    int n = levels[level].n_pts;
    struct A_csr* A = levels[level].A;

    /* pre-compute omega * D^{-1} */
    const float omega = 2.f/3.f;
    float* omegaDinv = create1dZeroVec(n);
    for (unsigned int i = 0; i < n; i++) {
        omegaDinv[i] = (omega / levels[level].diag[i]);
    }

    /* Jacobi pre-relax */
    for (unsigned int i = 0; i < mg->num_pre_relax; i++) {
        float* DinvAx = mpiMatVecInvDProductCSR1(levels[level].pmdat, x, *levels[level].A);

        /* Perform a weighted Jacobi sweep */
        for (unsigned int j = 0; j < n; j++) {
            x[j] = omegaDinv[j] * b[j] + x[j] - omega * DinvAx[j];
        }

        free(DinvAx);
    }

    /* Coarsen and recurse */
    /* if (level < mg->num_levels - 1) { */
    if (0) {
        /* Coarsen */
        unsigned int n_coarse = levels[level+1].n_pts;

        float* r_h = mpiMatVecProductCSR1(levels[level].pmdat, x, *levels[level].A); /* r_h = A * x */
        for (unsigned int i = 0; i < n; i++) {
            r_h[i] = b[i] - r_h[i]; /* r_h <- b - Ax */
        }
        float* r_H = mpiMatVecProductCSR1(levels[level + 1].pmdat, r_h, *levels[level + 1].R);
        float* e_H = create1dZeroVec(levels[level + 1].n_pts);

        mg_solve_recursive(pmd, mg, levels, level + 1, r_H, e_H);

        /* Interpolate and add error */
        float* e_h = mpiMatVecProductCSR1(levels[level].pmdat, e_H, *levels[level+1].P);
        for (unsigned int i = 0; i < n; i++) {
            x[i] = x[i] + e_h[i];
        }

        free(e_h);
        free(e_H);
        free(r_H);
        free(r_h);
    } else {
        for (unsigned int i = 0; i < mg->num_coarsest_relax; i++) {
            float* DinvAx = mpiMatVecInvDProductCSR1(levels[level].pmdat, x, *levels[level].A);

            /* Perform a weighted Jacobi sweep */
            for (unsigned int j = 0; j < n; j++) {
                x[j] = omegaDinv[j] * b[j] + x[j] - omega * DinvAx[j];
            }

            free(DinvAx);
        }
    }

    /* Jacobi post-relax */
    for (unsigned int i = 0; i < mg->num_pre_relax; i++) {
        float* DinvAx = mpiMatVecInvDProductCSR1(levels[level].pmdat, x, *levels[level].A);

        /* Perform a weighted Jacobi sweep */
        for (unsigned int j = 0; j < n; j++) {
            x[j] = omegaDinv[j] * b[j] + x[j] - omega * DinvAx[j];
        }

        free(DinvAx);
    }

    free(omegaDinv);
}

struct PCret mg_pc_solve(struct par_multdat pmd, struct PCdata pcdata, float* r){
    struct PCret pcret;
    /* Approximately solve Az = r */
    struct MGData* mg = pcdata.mg;
    struct MGLevelData* levels = mg->levels;

    float* z = create1dZeroVec(levels[0].n_pts);
    mg_solve_recursive(pmd, mg, levels, 0, r, z);
    pcret.sol = z;

    return pcret;
}


/*Preconditioner solve*/
struct PCret PC_Solve(struct par_multdat pmd, int n, float* r, struct A_csr A, struct A_csr L, struct A_csr U, char* pctype, struct PCdata pcdata){
    /* printf("PC_Solve .%s.\n", pctype); */
    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A,r, pcdata);
    }
    else if (strcmp(pctype, "MG2D") == 0 ||
             strcmp(pctype, "MG1D") == 0) {
        pcret = mg_pc_solve(pmd, pcdata, r);
    }
    else if (strcmp(pctype, "ILU")==0){
        pcret.sol = mpiTriangularSolveCSR1(pmd, r, L, U);
    } else if (strcmp(pctype, "Identity") == 0) {
        pcret.sol = copy_vec(n, r);
    }

    return pcret;
}


#endif
