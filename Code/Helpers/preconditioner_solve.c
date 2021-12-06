#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>
#include <stdbool.h>

#include "poisson_matrix_creaters.c"
#include "mg.c"
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

void setup_mg_pc(struct PCdata* pcdata, int n, struct A_csr* A, unsigned int num_levels, bool is_2d) {
    assert(num_levels >= 2);

    pcdata->mg = calloc(1, sizeof(struct MGData));
    pcdata->mg->num_levels = num_levels;
    pcdata->mg->levels = calloc(num_levels, sizeof(struct MGLevelData));
    pcdata->mg->num_pre_relax = 2;
    pcdata->mg->num_post_relax = 2;

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

    /* for (unsigned int i = 0; i < num_levels; i++) { */
    /*     printf("Level %u\n", i); */
    /*     printf("Number of points %u\n", levels[i].n_pts); */
    /*     printf("Number of points per dim %u\n", levels[i].n_pts_per_dim); */
    /* } */
}

// Setup data in PCdata struct for info to be extracted before solve
struct PCdata setupPCdata(struct PCdata pcdata, char* pctype, int n, struct A_csr A)
{
    if (strcmp(pctype,"Jacobi")==0){
        pcdata.diag = getCSRdiagonal(n,A);
    } else if (strcmp(pctype, "MG2D") == 0) {
        setup_mg_pc(&pcdata, n, &A, 2, 1);
    } else if (strcmp(pctype, "MG1D") == 0) {
        setup_mg_pc(&pcdata, n, &A, 2, 0);
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

    /* pre-compute omega * D^{-1} b */
    const float omega = 2.f/3.f;
    float* omegaDinvB = create1dZeroVec(n);
    for (unsigned int i = 0; i < n; i++) {
        omegaDinvB[i] = (b[i] / levels[level].diag[i]) * omega;
    }

    /* Jacobi pre-relax */
    for (unsigned int i = 0; i < mg->num_pre_relax; i++) {
        float* Ax = mpiMatVecProductCSR1(levels[level].pmdat, x, *levels[level].A);

        /* Perform a weighted Jacobi sweep */
        for (unsigned int j = 0; j < n; j++) {
            x[i] += omegaDinvB[i] - omega * (Ax[i] / levels[level].diag[i]);
        }

        free(Ax);
    }

    /* Coarsen and recurse */
    if (level < mg->num_levels - 1) {
        /* Coarsen */
        unsigned int n_coarse = levels[level+1].n_pts;
        float* x_H = mpiMatVecProductCSR1(levels[level + 1].pmdat, x, *levels[level+1].R);
        float* b_H = mpiMatVecProductCSR1(levels[level + 1].pmdat, b, *levels[level+1].R);
        mg_solve_recursive(pmd, mg, levels, level + 1, b_H, x_H);

        /* Interpolate and copy solution */
        float* x_h = mpiMatVecProductCSR1(levels[level].pmdat, x_H, *levels[level+1].P);
        for (unsigned int i = 0; i < n; i++) {
            x[i] = x_h[i];
        }

        free(x_H);
        free(b_H);
        free(x_h);
    }

    /* Jacobi post-relax */
    for (unsigned int i = 0; i < mg->num_pre_relax; i++) {
        float* Ax = mpiMatVecProductCSR1(levels[level].pmdat, x, *levels[level].A);

        /* Perform a weighted Jacobi sweep */
        for (unsigned int j = 0; j < n; j++) {
            x[i] += omegaDinvB[i] - omega * (Ax[i] / levels[level].diag[i]);
        }

        free(Ax);
    }
}

struct PCret mg_pc_solve(struct par_multdat pmd, struct PCdata pcdata, float* r){
    struct PCret pcret;
    /* Approximately solve Ar = 0 */
    struct MGData* mg = pcdata.mg;
    struct MGLevelData* levels = mg->levels;

    float* z = create1dZeroVec(levels[0].n_pts);
    mg_solve_recursive(pmd, mg, levels, 0, r, z);
    pcret.sol = z;

    /* printf("%d\n", levels[0].n_pts); */
    /* print_vec(levels[1].n_pts, z); */

    return pcret;
}


/*Preconditioner solve*/
struct PCret PC_Solve(struct par_multdat pmd, int n, float* r, struct A_csr A, char* pctype, struct PCdata pcdata){
    /* printf("PC_Solve .%s.\n", pctype); */
    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A,r, pcdata);
    } else if (strcmp(pctype, "MG2D") == 0 ||
               strcmp(pctype, "MG1D") == 0) {
        pcret = mg_pc_solve(pmd, pcdata, r);
    }

    return pcret;
}


#endif
